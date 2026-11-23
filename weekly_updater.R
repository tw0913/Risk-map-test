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


# ---- 1. API Credentials & Dynamic Dates ----
# [ACTION REQUIRED BY B+LNZ]: 
# These credentials authorize the script to download the weekly ERA5 weather data.
# To obtain them:
# 1. Register a free account at the Copernicus Climate Data Store (CDS): https://cds.climate.copernicus.eu/
# 2. Log in, navigate to your user profile, and scroll down to the "API key" section.
# 3. Replace the email below with your registered CDS email address.
# 4. Replace "your_api_key_here" with your Personal API Token (or use Sys.getenv("ECMWF_API_KEY") if using GitHub Secrets)

wf_set_key(user = "twang@kreavet.com", key = "98465c3a-9a6b-48e9-9118-acb9abcc232f")
Sys.setenv(WF_USER = "twang@kreavet.com")

# TEMPORARY TEST BASELINE (1 WEEK EARLIER)
end_date <- as.Date("2026-11-08")
start_date <- as.Date("2026-07-01")
is_reset_week <- TRUE

print(paste("Fetching ERA5 data from", start_date, "to", end_date))

# ---- 2. Download Weekly ERA5 Proxy Data ----
# Name the file based on the end date so each week gets a unique file
target_file <- paste0("NZ_ERA5_Weekly_", format(end_date, "%Y_%m_%d"), ".nc")

if (!file.exists(target_file)) {
  request_weekly <- list(
    dataset_short_name = "reanalysis-era5-land",
    variable = c("2m_temperature", "total_precipitation"),
    
    # Passing the exact dates handles month-crossovers automatically!
    date = paste(start_date, end_date, sep = "/"),
    time = sprintf("%02d:00", 0:23), 
    area = c(-34, 166, -48, 179),
    format = "netcdf",
    target = target_file
  )
  
  wf_request(user = "twang@kreavet.com", request = request_weekly, transfer = TRUE, path = getwd())
  
  zip_version <- sub("\\.nc$", ".zip", target_file)
  if(file.exists(zip_version)) {
    internal_file <- unzip(zip_version, list = TRUE)$Name[1]
    unzip(zip_version, overwrite = TRUE)
    file.rename(internal_file, target_file)
    file.remove(zip_version)
  }
}

# ---- 3. The Calibrated NZ GLOWORM-FL Engine & Risk Scaler ----
flsim = function (init, start, end, input, temp, precip, species, lat, C_window = NULL, u5_mult = 1, m2_mult = 1) 
{
  date.range = seq(as.Date(start), as.Date(end), "days")
  if (length(date.range) != length(temp)) stop("Error: climatic data and dates mismatch.")
  
  global.t = seq(1, length(date.range))
  photoperiod = daylength(lat = lat, doy = date.range)
  
  if (species == "haemonchus") {
    c_val <- ifelse(is.null(C_window), 4, round(C_window))
    Pt = (4.95 * exp(0.062 * temp))/100
    PET = 0.55 * ((photoperiod/12)^2) * Pt * 25.4
    PET_win = NA; PET_m_win = NA; Precip_win = NA; Precip_m_win = NA
    
    for (z in 1:length(PET)) {
      min = z; max = ifelse((z + c_val) > length(PET), length(PET), z + c_val)
      PET_win[z] = sum(x = PET[min:max]); Precip_win[z] = sum(x = precip[min:max])
      max = z; min = ifelse((z - c_val) < 1, 1, z - c_val)
      PET_m_win[z] = sum(x = PET[min:max]); Precip_m_win[z] = sum(x = precip[min:max])
      P.E.1 = Precip_win/PET_win; P.E.2 = Precip_m_win/PET_m_win
    }
    dev.1 = pmin(1, pmax(0, -0.09746 + 0.01063 * temp))
    mu.1 = pmin(1, exp(-1.47135 - 0.11444 * temp + 0.00327 * (temp^2)))
    mu.2 = pmin(1, exp(-1.823 - 0.1418 * temp + 0.00405 * (temp^2)))
    mu.3 = pmin(1, exp(-2.6308 - 0.14407 * temp + 0.00463 * (temp^2)))
    mu.4 = pmin(1, exp(-3.68423 - 0.2535 * temp + 0.0074 * (temp^2)))
    mu.5 = pmin(1, mu.3 * 6) 
    
    h.mig = ifelse(precip >= 2, 0.25, ifelse(P.E.2 >= 1, 0.051, 0))
    v.mig = pmin(1, (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * (temp^2)))) * m2_mult)
    
    dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
    mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
    mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
    mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
    mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
    mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
    hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
    vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
    
    init[5] = (init[4] * (1/v.mig[1])) - init[4]
    egg.correction = ifelse(P.E.1 < 1, 0.1, 1)
    init[1] = init[1] * egg.correction[1]
  }
  
  if (species == "teladorsagia") {
    c_val <- ifelse(is.null(C_window), 9, round(C_window))
    Pt = (4.95 * exp(0.062 * temp))/100
    PET = 0.55 * ((photoperiod/12)^2) * Pt * 25.4
    PET_win = NA; PET_m_win = NA; Precip_win = NA; Precip_m_win = NA
    
    for (z in 1:length(PET)) {
      min = z; max = ifelse((z + c_val) > length(PET), length(PET), z + c_val)
      PET_win[z] = sum(x = PET[min:max]); Precip_win[z] = sum(x = precip[min:max])
      max = z; min = ifelse((z - c_val) < 1, 1, z - c_val)
      PET_m_win[z] = sum(x = PET[min:max]); Precip_m_win[z] = sum(x = precip[min:max])
      P.E.1 = Precip_win/PET_win; P.E.2 = Precip_m_win/PET_m_win
    }
    dev.1 = pmin(1, pmax(0, -0.02085 + 0.00467 * temp))
    mu.1 = pmin(1, exp(-1.62026 - 0.17771 * temp + 0.00629 * (temp^2)))
    mu.2 = mu.1
    mu.4 = pmin(1, exp(-4.58817 - 0.13996 * temp + 0.00461 * (temp^2)))
    mu.3 = pmin(1, 10 * mu.4)
    mu.5 = pmin(1, mu.3 * u5_mult)
    
    h.mig = ifelse(precip >= 2, 0.21, ifelse(P.E.2 >= 1, 0.025, 0))
    v.mig = pmin(1, (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * (temp^2)))) * m2_mult)
    
    dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
    mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
    mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
    mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
    mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
    mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
    hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
    vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
    
    init[5] = (init[4] * (1/v.mig[1])) - init[4]
    egg.correction = ifelse(P.E.1 < 1, 0.1, 1)
    init[1] = init[1] * egg.correction[1]
  }
  
  event = data.frame(var = "E", time = global.t, value = input * egg.correction, method = "add")
  para.dyn = function(t, para.init, para.par) {
    with(as.list(c(para.init, para.par)), {
      dev1 = dev1rate(t); mu1 = mu1rate(t); mu2 = mu2rate(t)
      mu3 = mu3rate(t); mu4 = mu4rate(t); mu5 = mu5rate(t)
      m1 = hmigrate(t); m2 = vmigrate(t)
      dE = -(dev1 * 2 + mu1) * E
      dL = -(dev1 * 2 + mu2) * L + (dev1 * 2) * E
      dL3 = -(mu3 + m1) * L3 + (dev1 * 2) * L
      dPasture = -mu4 * (Pasture * (1 - m2)) - mu5 * (Pasture * m2) + m1 * L3
      return(list(c(dE = dE, dL = dL, dL3 = dL3, dPasture = dPasture)))
    })
  }
  para.init = c(E = init[1], L = init[2], L3 = init[3], Pasture = (init[4] + init[5]))
  para.sol = lsoda(y = para.init, times = global.t, func = para.dyn, parms = NULL, events = list(data = event))
  Soil = para.sol[, "Pasture"] * (1 - (vmigrate(global.t)))
  Herbage = para.sol[, "Pasture"] * vmigrate(global.t)
  return(cbind(para.sol, Soil, Herbage))
}

calculate_risk_score <- function(l3, species) {
  if (is.na(l3)) return(NA)
  
  if (species == "teladorsagia") {
    if (l3 < 50) { score <- (l3 / 50) * 25 }
    else if (l3 < 200) { score <- 25 + ((l3 - 50) / 150) * 25 }
    else if (l3 < 400) { score <- 50 + ((l3 - 200) / 200) * 25 }
    else { score <- 75 + min(25, ((l3 - 400) / 180) * 25) }
    
  } else if (species == "haemonchus") {
    if (l3 < 25) { score <- (l3 / 25) * 25 }
    else if (l3 < 100) { score <- 25 + ((l3 - 25) / 75) * 25 }
    else if (l3 < 200) { score <- 50 + ((l3 - 100) / 100) * 25 }
    else { score <- 75 + min(25, ((l3 - 200) / 90) * 25) }
  }
  return(round(score))
}

# ---- 4. Spatial Grid Execution ----
print("Fetching New Zealand geographic boundary and generating dashboard coordinates...")
nz_boundary <- gadm(country = "NZ", level = 0, path = tempdir())
grid_template <- rast(nz_boundary, resolution = 0.05)
nz_raster <- rasterize(nz_boundary, grid_template, field = 1)
land_coords <- crds(nz_raster, df = TRUE, na.rm = TRUE)

dir.create("processed_data", showWarnings = FALSE)
saveRDS(land_coords, "processed_data/grid_coords.rds")

# ---- 5. Climate Extraction & Simulation Loop ----
print("Loading ERA5 weather data...")
climate_raster <- rast(target_file)
temp_raster <- climate_raster["t2m"]
precip_raster <- climate_raster["tp"]

print("Aggregating hourly to daily...")
layer_dates <- as.Date(time(temp_raster))
daily_temp <- tapp(temp_raster, layer_dates, fun = mean, na.rm = TRUE) - 273.15
daily_precip <- tapp(precip_raster, layer_dates, fun = mean, na.rm = TRUE) * 24 * 1000

temp_ext <- terra::extract(daily_temp, land_coords)
precip_ext <- terra::extract(daily_precip, land_coords)
temp_ext$ID <- NULL; precip_ext$ID <- NULL

# Calculate the dynamic seasonal shedding curve for the exact daily window
print("Calculating seasonal pasture contamination inputs...")
monthly_shedding <- c(
  "01" = 80, "02" = 80, "03" = 60, "04" = 40, 
  "05" = 10, "06" = 5,  "07" = 5,  "08" = 50, 
  "09" = 100, "10" = 100, "11" = 80, "12" = 80
)
daily_dates <- unique(as.Date(time(temp_raster)))
sim_months <- format(daily_dates, "%m")
dynamic_input <- as.numeric(monthly_shedding[sim_months])

state_file <- "processed_data/model_state.rds"
history_file <- "processed_data/risk_history.rds"
reservoir_file <- "processed_data/reservoir_history.rds"

# Implement Annual Reset or Load State
if(is_reset_week) {
  print("Annual July 1st Reset: Initializing fresh biological baseline...")
  baseline <- c(100, 0, 0, 0, 0) 
  prev_state <- list(
    H = matrix(rep(baseline, nrow(land_coords)), ncol = 5, byrow = TRUE),
    T = matrix(rep(baseline, nrow(land_coords)), ncol = 5, byrow = TRUE)
  )
  history_db <- list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
  reservoir_db <- list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
} else if(file.exists(state_file)) {
  print("Loading previous week's biological state...")
  prev_state <- readRDS(state_file)
  history_db <- if(file.exists(history_file)) readRDS(history_file) else list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
  reservoir_db <- if(file.exists(reservoir_file)) readRDS(reservoir_file) else list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
} else {
  baseline <- c(100, 0, 0, 0, 0) 
  prev_state <- list(
    H = matrix(rep(baseline, nrow(land_coords)), ncol = 5, byrow = TRUE),
    T = matrix(rep(baseline, nrow(land_coords)), ncol = 5, byrow = TRUE)
  )
  history_db <- list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
  reservoir_db <- list(H = data.frame(row.names = 1:nrow(land_coords)), T = data.frame(row.names = 1:nrow(land_coords)))
}

new_state <- list(
  H = matrix(NA, nrow = nrow(land_coords), ncol = 5),
  T = matrix(NA, nrow = nrow(land_coords), ncol = 5)
)

risk_scores_H <- numeric(nrow(land_coords))
risk_scores_T <- numeric(nrow(land_coords))
reservoir_scores_H <- numeric(nrow(land_coords))
reservoir_scores_T <- numeric(nrow(land_coords))

print("Running GLOWORM-FL simulation...")

# Sync the simulation window to match the exact layers in the climate raster
start_date <- min(daily_dates, na.rm = TRUE)
end_date <- max(daily_dates, na.rm = TRUE)

for(i in 1:nrow(land_coords)) {
  if(i %% 1000 == 0) print(paste("Simulating square", i, "of", nrow(land_coords)))
  
  t_vec <- as.numeric(temp_ext[i, ])
  p_vec <- as.numeric(precip_ext[i, ])
  lat <- land_coords[i, "y"]
  
  # Ensure prev_state is properly populated before checking missing climate data
  if(any(is.na(prev_state$H[i, ]))) prev_state$H[i, ] <- c(100, 0, 0, 0, 0)
  if(any(is.na(prev_state$T[i, ]))) prev_state$T[i, ] <- c(100, 0, 0, 0, 0)
  
  # Prevent silent reset: If grid is missing weather, map NA but preserve biological state for next week
  if(any(is.na(t_vec)) || any(is.na(p_vec))) {
    risk_scores_H[i] <- NA; risk_scores_T[i] <- NA
    new_state$H[i, ] <- prev_state$H[i, ]
    new_state$T[i, ] <- prev_state$T[i, ]
    next
  }
  
  sim_H <- flsim(init = prev_state$H[i, ], start = start_date, end = end_date, 
                 input = dynamic_input, temp = t_vec, precip = p_vec, 
                 species = "haemonchus", lat = lat)
  
  new_state$H[i, ] <- as.numeric(sim_H[nrow(sim_H), c("E", "L", "L3", "Soil", "Herbage")])
  raw_L3_H <- sim_H[nrow(sim_H), "Herbage"]
  risk_scores_H[i] <- calculate_risk_score(raw_L3_H, "haemonchus")
  reservoir_scores_H[i] <- sim_H[nrow(sim_H), "Soil"] + raw_L3_H
  
  sim_T <- flsim(init = prev_state$T[i, ], start = start_date, end = end_date, 
                 input = dynamic_input, temp = t_vec, precip = p_vec, 
                 species = "teladorsagia", lat = lat)
  
  new_state$T[i, ] <- as.numeric(sim_T[nrow(sim_T), c("E", "L", "L3", "Soil", "Herbage")])
  raw_L3_T <- sim_T[nrow(sim_T), "Herbage"]
  risk_scores_T[i] <- calculate_risk_score(raw_L3_T, "teladorsagia")
  reservoir_scores_T[i] <- sim_T[nrow(sim_T), "Soil"] + raw_L3_T
}

print("Updating states and databases...")
saveRDS(new_state, state_file)

current_date <- as.character(end_date)
history_db$H[[current_date]] <- risk_scores_H
history_db$T[[current_date]] <- risk_scores_T
reservoir_db$H[[current_date]] <- reservoir_scores_H
reservoir_db$T[[current_date]] <- reservoir_scores_T

saveRDS(history_db, history_file)
saveRDS(reservoir_db, reservoir_file)

print("Saving GeoTIFF maps...")
final_map_H <- rasterize(as.matrix(land_coords), grid_template, values = risk_scores_H)
final_map_T <- rasterize(as.matrix(land_coords), grid_template, values = risk_scores_T)

writeRaster(final_map_H, paste0("processed_data/NZ_Risk_Map_H_contortus_", Sys.Date(), ".tif"), overwrite = TRUE)
writeRaster(final_map_T, paste0("processed_data/NZ_Risk_Map_T_circumcincta_", Sys.Date(), ".tif"), overwrite = TRUE)

print("Success! Update complete.")


