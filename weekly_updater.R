# ==============================================================================
# B+LNZ: AUTOMATED WEEKLY RISK ENGINE (Annual July 1 Reset)
# ==============================================================================
library(ecmwfr)
library(terra)
library(dplyr)
library(stringr)
library(lubridate)
library(geosphere)
library(deSolve)
library(zoo)
library(geodata)

# ---- 0. Manual Testing Controls (Keep NULL / FALSE for normal automation) ----
MANUAL_END_DATE <- NULL   # Example for test Step 1: as.Date("2026-11-08")
FORCE_RESET     <- FALSE  # Set TRUE only when forcing a fresh July 1 baseline

# ---- 1. API Credentials ----
# [ACTION REQUIRED BY B+LNZ]:
# 1. Register a free account at Copernicus CDS: https://cds.climate.copernicus.eu/
# 2. Add your CDS API key as a GitHub Secret named ECMWF_API_KEY.
# 3. Replace the default email below (or set an ECMWF_USER secret).
cds_key <- Sys.getenv("ECMWF_API_KEY")
if (!nzchar(cds_key)) {
  stop("ECMWF_API_KEY environment variable is missing. Set it in GitHub Secrets.")
}
cds_user <- Sys.getenv("ECMWF_USER", unset = "twang@kreavet.com")
wf_set_key(user = cds_user, key = cds_key)
Sys.setenv(WF_USER = cds_user)

# ---- 2. Self-Healing Date Window (Overlap-Proof) ----
dir.create("processed_data", showWarnings = FALSE)
state_file     <- "processed_data/model_state.rds"
history_file   <- "processed_data/risk_history.rds"
reservoir_file <- "processed_data/reservoir_history.rds"

end_date <- if (!is.null(MANUAL_END_DATE)) as.Date(MANUAL_END_DATE) else (Sys.Date() - 8)

# Determine the most recent July 1st on or before end_date (Southern Hemisphere season)
end_year <- as.integer(format(end_date, "%Y"))
season_start <- as.Date(paste0(end_year, "-07-01"))
if (end_date < season_start) {
  season_start <- as.Date(paste0(end_year - 1, "-07-01"))
}

has_saved_state <- file.exists(state_file) && !FORCE_RESET
prev_state_peek <- if (has_saved_state) readRDS(state_file) else NULL
last_saved_date <- if (!is.null(prev_state_peek) && !is.null(prev_state_peek$last_date)) {
  as.Date(prev_state_peek$last_date)
} else {
  NULL
}

# Decide whether to reset to July 1 or continue seamlessly from last_saved_date + 1
if (!has_saved_state || is.null(last_saved_date) || last_saved_date < season_start) {
  is_reset_week <- TRUE
  start_date <- season_start
  print(paste("Reset / fresh season run from", start_date, "to", end_date))
} else {
  if (end_date <= last_saved_date) {
    print(paste("Already up to date through", last_saved_date, "- skipping run to prevent double-counting."))
    quit(save = "no", status = 0)
  }
  is_reset_week <- FALSE
  start_date <- last_saved_date + 1
  print(paste("Continuation run from", start_date, "to", end_date))
}

# ---- 3. Download ERA5-Land Climate Data (into temporary directory) ----
download_dir <- tempdir()
target_name  <- paste0("NZ_ERA5_Weekly_", format(end_date, "%Y_%m_%d"), ".nc")
target_path  <- file.path(download_dir, target_name)

if (!file.exists(target_path)) {
  request_weekly <- list(
    dataset_short_name = "reanalysis-era5-land",
    variable = c("2m_temperature", "total_precipitation"),
    date = paste(start_date, end_date, sep = "/"),
    time = sprintf("%02d:00", 0:23),
    area = c(-34, 166, -48, 179),
    format = "netcdf",
    target = target_name
  )

  wf_request(user = cds_user, request = request_weekly, transfer = TRUE, path = download_dir)

  zip_version <- sub("\\.nc$", ".zip", target_path)
  if (file.exists(zip_version)) {
    internal_file <- unzip(zip_version, list = TRUE)$Name[1]
    unzip(zip_version, exdir = download_dir, overwrite = TRUE)
    file.rename(file.path(download_dir, internal_file), target_path)
    file.remove(zip_version)
  }
}

# ---- 4. The Calibrated NZ GLOWORM-FL Engine & Risk Scaler ----
safe_approxfun <- function(v) {
  if (length(v) == 1) {
    force(v)
    return(function(t) rep(v[1], length(t)))
  }
  approxfun(x = v, method = "linear", rule = 2)
}

flsim <- function(init, start, end, input, temp, precip, species, lat,
                  is_continuation = FALSE, C_window = NULL, u5_mult = 1, m2_mult = 1) {
  date.range <- seq(as.Date(start), as.Date(end), "days")
  if (length(date.range) != length(temp)) stop("Error: climatic data and dates mismatch.")

  global.t <- seq(1, length(date.range))
  photoperiod <- daylength(lat = lat, doy = date.range)

  if (species == "haemonchus") {
    c_val <- ifelse(is.null(C_window), 4, round(C_window))
    Pt <- (4.95 * exp(0.062 * temp)) / 100
    PET <- 0.55 * ((photoperiod / 12)^2) * Pt * 25.4
    PET_win <- numeric(length(PET)); PET_m_win <- numeric(length(PET))
    Precip_win <- numeric(length(PET)); Precip_m_win <- numeric(length(PET))

    for (z in seq_along(PET)) {
      mn <- z; mx <- min(length(PET), z + c_val)
      PET_win[z] <- sum(PET[mn:mx]); Precip_win[z] <- sum(precip[mn:mx])
      mx <- z; mn <- max(1, z - c_val)
      PET_m_win[z] <- sum(PET[mn:mx]); Precip_m_win[z] <- sum(precip[mn:mx])
    }
    P.E.1 <- ifelse(PET_win > 0, Precip_win / PET_win, 1)
    P.E.2 <- ifelse(PET_m_win > 0, Precip_m_win / PET_m_win, 1)

    dev.1 <- pmin(1, pmax(0, -0.09746 + 0.01063 * temp))
    mu.1  <- pmin(1, exp(-1.47135 - 0.11444 * temp + 0.00327 * (temp^2)))
    mu.2  <- pmin(1, exp(-1.823   - 0.1418  * temp + 0.00405 * (temp^2)))
    mu.3  <- pmin(1, exp(-2.6308  - 0.14407 * temp + 0.00463 * (temp^2)))
    mu.4  <- pmin(1, exp(-3.68423 - 0.2535  * temp + 0.0074  * (temp^2)))
    mu.5  <- pmin(1, mu.3 * 6)

    h.mig <- ifelse(precip >= 2, 0.25, ifelse(P.E.2 >= 1, 0.051, 0))
    v.mig <- pmin(1, pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * (temp^2))) * m2_mult)
    egg.correction <- ifelse(P.E.1 < 1, 0.1, 1)
  }

  if (species == "teladorsagia") {
    c_val <- ifelse(is.null(C_window), 9, round(C_window))
    Pt <- (4.95 * exp(0.062 * temp)) / 100
    PET <- 0.55 * ((photoperiod / 12)^2) * Pt * 25.4
    PET_win <- numeric(length(PET)); PET_m_win <- numeric(length(PET))
    Precip_win <- numeric(length(PET)); Precip_m_win <- numeric(length(PET))

    for (z in seq_along(PET)) {
      mn <- z; mx <- min(length(PET), z + c_val)
      PET_win[z] <- sum(PET[mn:mx]); Precip_win[z] <- sum(precip[mn:mx])
      mx <- z; mn <- max(1, z - c_val)
      PET_m_win[z] <- sum(PET[mn:mx]); Precip_m_win[z] <- sum(precip[mn:mx])
    }
    P.E.1 <- ifelse(PET_win > 0, Precip_win / PET_win, 1)
    P.E.2 <- ifelse(PET_m_win > 0, Precip_m_win / PET_m_win, 1)

    dev.1 <- pmin(1, pmax(0, -0.02085 + 0.00467 * temp))
    mu.1  <- pmin(1, exp(-1.62026 - 0.17771 * temp + 0.00629 * (temp^2)))
    mu.2  <- mu.1
    mu.4  <- pmin(1, exp(-4.58817 - 0.13996 * temp + 0.00461 * (temp^2)))
    mu.3  <- pmin(1, 10 * mu.4)
    mu.5  <- pmin(1, mu.3 * u5_mult)

    h.mig <- ifelse(precip >= 2, 0.21, ifelse(P.E.2 >= 1, 0.025, 0))
    v.mig <- pmin(1, pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * (temp^2))) * m2_mult)
    egg.correction <- ifelse(P.E.1 < 1, 0.1, 1)
  }

  dev1rate <- safe_approxfun(dev.1)
  mu1rate  <- safe_approxfun(mu.1)
  mu2rate  <- safe_approxfun(mu.2)
  mu3rate  <- safe_approxfun(mu.3)
  mu4rate  <- safe_approxfun(mu.4)
  mu5rate  <- safe_approxfun(mu.5)
  hmigrate <- safe_approxfun(h.mig)
  vmigrate <- safe_approxfun(v.mig)

  # Only apply initial egg correction on a fresh baseline (carried-over eggs were already corrected)
  if (!is_continuation) {
    init[1] <- init[1] * egg.correction[1]
  }

  event <- data.frame(var = "E", time = global.t, value = input * egg.correction, method = "add")

  para.dyn <- function(t, para.init, para.par) {
    with(as.list(c(para.init, para.par)), {
      dev1 <- dev1rate(t); mu1 <- mu1rate(t); mu2 <- mu2rate(t)
      mu3  <- mu3rate(t);  mu4 <- mu4rate(t); mu5 <- mu5rate(t)
      m1   <- hmigrate(t); m2  <- vmigrate(t)
      dE <- -(dev1 * 2 + mu1) * E
      dL <- -(dev1 * 2 + mu2) * L + (dev1 * 2) * E
      dL3 <- -(mu3 + m1) * L3 + (dev1 * 2) * L
      dPasture <- -mu4 * (Pasture * (1 - m2)) - mu5 * (Pasture * m2) + m1 * L3
      return(list(c(dE = dE, dL = dL, dL3 = dL3, dPasture = dPasture)))
    })
  }

  para.init <- c(E = init[1], L = init[2], L3 = init[3], Pasture = (init[4] + init[5]))

  # On continuation runs, start ODE at t = 0 so a full 24h elapses between last week's end and Day 1
  sim_times <- if (is_continuation) c(0, global.t) else global.t
  para.sol <- lsoda(y = para.init, times = sim_times, func = para.dyn, parms = NULL, events = list(data = event))
  if (is_continuation) {
    para.sol <- para.sol[-1, , drop = FALSE]
  }

  Soil    <- para.sol[, "Pasture"] * (1 - vmigrate(global.t))
  Herbage <- para.sol[, "Pasture"] * vmigrate(global.t)
  return(cbind(para.sol, Soil, Herbage))
}

calculate_risk_score <- function(l3, species) {
  if (is.na(l3)) return(NA)
  if (species == "teladorsagia") {
    if (l3 < 50)       { score <- (l3 / 50) * 25 }
    else if (l3 < 200) { score <- 25 + ((l3 - 50) / 150) * 25 }
    else if (l3 < 400) { score <- 50 + ((l3 - 200) / 200) * 25 }
    else               { score <- 75 + min(25, ((l3 - 400) / 180) * 25) }
  } else if (species == "haemonchus") {
    if (l3 < 25)       { score <- (l3 / 25) * 25 }
    else if (l3 < 100) { score <- 25 + ((l3 - 25) / 75) * 25 }
    else if (l3 < 200) { score <- 50 + ((l3 - 100) / 100) * 25 }
    else               { score <- 75 + min(25, ((l3 - 200) / 90) * 25) }
  }
  return(round(score))
}

# ---- 5. Spatial Grid Execution ----
print("Fetching New Zealand geographic boundary and generating dashboard coordinates...")
nz_boundary   <- gadm(country = "NZ", level = 0, path = tempdir())
grid_template <- rast(nz_boundary, resolution = 0.05)
nz_raster     <- rasterize(nz_boundary, grid_template, field = 1)
land_coords   <- crds(nz_raster, df = TRUE, na.rm = TRUE)

saveRDS(land_coords, "processed_data/grid_coords.rds")

# ---- 6. Climate Extraction & Alignment ----
print("Loading ERA5 weather data...")
climate_raster <- rast(target_path)
temp_raster    <- climate_raster["t2m"]
precip_raster  <- climate_raster["tp"]

print("Aggregating hourly to daily...")
layer_dates <- as.Date(time(temp_raster))
keep_layers <- which(layer_dates >= start_date & layer_dates <= end_date)
temp_raster   <- temp_raster[[keep_layers]]
precip_raster <- precip_raster[[keep_layers]]
layer_dates   <- layer_dates[keep_layers]

daily_dates <- sort(unique(layer_dates))
date_index  <- match(layer_dates, daily_dates)

# Note: Preserved to match original 3-year threshold calibration; do not alter without recalibrating cutoffs
daily_temp   <- tapp(temp_raster, date_index, fun = mean, na.rm = TRUE) - 273.15
daily_precip <- tapp(precip_raster, date_index, fun = mean, na.rm = TRUE) * 24 * 1000

temp_ext   <- terra::extract(daily_temp, land_coords)
precip_ext <- terra::extract(daily_precip, land_coords)
temp_ext$ID <- NULL; precip_ext$ID <- NULL

monthly_shedding <- c(
  "01" = 80, "02" = 80, "03" = 60, "04" = 40,
  "05" = 10, "06" = 5,  "07" = 5,  "08" = 50,
  "09" = 100, "10" = 100, "11" = 80, "12" = 80
)
dynamic_input <- as.numeric(monthly_shedding[format(daily_dates, "%m")])

# ---- 7. Load or Initialize Biological State ----
if (is_reset_week) {
  print("Initializing fresh biological baseline...")
  baseline <- c(100, 0, 0, 0, 0)
  prev_state <- list(
    H = matrix(rep(baseline, nrow(land_coords)), ncol = 5, byrow = TRUE),
    T = matrix(rep(baseline, nrow(land_coords)), ncol = 5, byrow = TRUE)
  )
  history_db   <- list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
  reservoir_db <- list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
} else {
  print("Loading previous week's biological state...")
  prev_state   <- prev_state_peek
  history_db   <- if (file.exists(history_file)) readRDS(history_file) else list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
  reservoir_db <- if (file.exists(reservoir_file)) readRDS(reservoir_file) else list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
}

new_state <- list(
  H = matrix(NA_real_, nrow = nrow(land_coords), ncol = 5),
  T = matrix(NA_real_, nrow = nrow(land_coords), ncol = 5),
  last_date = max(daily_dates, na.rm = TRUE)
)

risk_scores_H      <- rep(NA_real_, nrow(land_coords))
risk_scores_T      <- rep(NA_real_, nrow(land_coords))
reservoir_scores_H <- rep(NA_real_, nrow(land_coords))
reservoir_scores_T <- rep(NA_real_, nrow(land_coords))

print("Running GLOWORM-FL simulation...")
sim_start <- min(daily_dates, na.rm = TRUE)
sim_end   <- max(daily_dates, na.rm = TRUE)
is_cont   <- !is_reset_week

for (i in 1:nrow(land_coords)) {
  if (i %% 1000 == 0) print(paste("Simulating square", i, "of", nrow(land_coords)))

  t_vec <- as.numeric(temp_ext[i, ])
  p_vec <- as.numeric(precip_ext[i, ])
  lat   <- land_coords[i, "y"]

  h_was_na <- any(is.na(prev_state$H[i, ]))
  t_was_na <- any(is.na(prev_state$T[i, ]))
  if (h_was_na) prev_state$H[i, ] <- c(100, 0, 0, 0, 0)
  if (t_was_na) prev_state$T[i, ] <- c(100, 0, 0, 0, 0)

  # Preserve biological state and record NA if weather is missing
  if (any(is.na(t_vec)) || any(is.na(p_vec))) {
    new_state$H[i, ] <- prev_state$H[i, ]
    new_state$T[i, ] <- prev_state$T[i, ]
    next
  }

  sim_H <- flsim(
    init = prev_state$H[i, ], start = sim_start, end = sim_end,
    input = dynamic_input, temp = t_vec, precip = p_vec,
    species = "haemonchus", lat = lat, is_continuation = (is_cont && !h_was_na)
  )
  new_state$H[i, ]      <- as.numeric(sim_H[nrow(sim_H), c("E", "L", "L3", "Soil", "Herbage")])
  risk_scores_H[i]      <- calculate_risk_score(sim_H[nrow(sim_H), "Herbage"], "haemonchus")
  reservoir_scores_H[i] <- sim_H[nrow(sim_H), "Soil"]

  sim_T <- flsim(
    init = prev_state$T[i, ], start = sim_start, end = sim_end,
    input = dynamic_input, temp = t_vec, precip = p_vec,
    species = "teladorsagia", lat = lat, is_continuation = (is_cont && !t_was_na)
  )
  new_state$T[i, ]      <- as.numeric(sim_T[nrow(sim_T), c("E", "L", "L3", "Soil", "Herbage")])
  risk_scores_T[i]      <- calculate_risk_score(sim_T[nrow(sim_T), "Herbage"], "teladorsagia")
  reservoir_scores_T[i] <- sim_T[nrow(sim_T), "Soil"]
}

# ---- 8. Save State, History, and Latest GeoTIFFs ----
print("Updating states and databases...")
saveRDS(new_state, state_file)

current_date <- as.character(sim_end)
history_db$H[[current_date]]   <- risk_scores_H
history_db$T[[current_date]]   <- risk_scores_T
reservoir_db$H[[current_date]] <- reservoir_scores_H
reservoir_db$T[[current_date]] <- reservoir_scores_T

saveRDS(history_db, history_file)
saveRDS(reservoir_db, reservoir_file)

print("Saving latest GeoTIFF maps...")
# Remove older dated GeoTIFFs so app.R always loads the single latest map
old_tifs <- list.files("processed_data", pattern = "\\.tif$", full.names = TRUE)
if (length(old_tifs) > 0) file.remove(old_tifs)

final_map_H <- rasterize(as.matrix(land_coords), grid_template, values = risk_scores_H)
final_map_T <- rasterize(as.matrix(land_coords), grid_template, values = risk_scores_T)

writeRaster(final_map_H, "processed_data/NZ_Risk_Map_H_contortus_latest.tif", overwrite = TRUE)
writeRaster(final_map_T, "processed_data/NZ_Risk_Map_T_circumcincta_latest.tif", overwrite = TRUE)

print(paste("Success! Update complete through", current_date))
