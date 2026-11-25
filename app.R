# ============================================================
# B+LNZ Parasite Risk Dashboard
# Bulletproof Base Shiny Version (Bypasses App Control Blocks)
# ============================================================

library(shiny)
library(leaflet)
library(terra)
library(plotly)
library(tidygeocoder)

# Keep server-side errors visible while the demo is being validated.
options(shiny.sanitize.errors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

DATA_DIR <- "processed_data"
RISK_HISTORY_FILE <- file.path(DATA_DIR, "risk_history.rds")
RESERVOIR_HISTORY_FILE <- file.path(DATA_DIR, "reservoir_history.rds")
GRID_COORDS_FILE <- file.path(DATA_DIR, "grid_coords.rds")

FILE_POLL_MS <- 60 * 1000
MAX_GRID_DISTANCE_DEG <- 0.08

RESERVOIR_HISTORY_IS_TOTAL <- FALSE

SPECIES_LABELS <- c(
  H = "Haemonchus contortus",
  T = "Teladorsagia circumcincta"
)

SPECIES_RASTER_PATTERNS <- c(
  H = "H_contortus",
  T = "T_circumcincta"
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

fail <- function(message) {
  structure(list(ok = FALSE, message = message), class = "dashboard_failure")
}

ok <- function(...) {
  c(list(ok = TRUE), list(...))
}

is_fail <- function(x) {
  inherits(x, "dashboard_failure") || (is.list(x) && identical(x$ok, FALSE))
}

safe_read_rds <- function(path) {
  if (!file.exists(path)) {
    return(fail(paste("File not found:", path)))
  }
  
  tryCatch(
    ok(value = readRDS(path)),
    error = function(e) fail(paste("Could not read", basename(path), "-", conditionMessage(e)))
  )
}

file_signature <- function(path) {
  if (!file.exists(path)) return("MISSING")
  i <- file.info(path)
  paste(as.numeric(i$mtime), i$size, sep = "_")
}

parse_history_dates <- function(x) {
  if (is.null(x)) return(as.Date(character(0)))
  
  x <- as.character(x)
  x <- sub("^X(?=\\d{4}[._-]\\d{2}[._-]\\d{2}$)", "", x, perl = TRUE)
  x <- gsub("[._]", "-", x)
  
  suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
}

date_column_map <- function(x) {
  nms <- colnames(x)
  
  if (is.null(nms) || length(nms) == 0) {
    return(data.frame(column = character(0), date = as.Date(character(0))))
  }
  
  d <- parse_history_dates(nms)
  keep <- !is.na(d)
  
  out <- data.frame(
    column = nms[keep],
    date = d[keep],
    stringsAsFactors = FALSE
  )
  
  if (nrow(out) > 0) {
    out <- out[order(out$date), , drop = FALSE]
    out <- out[!duplicated(out$date), , drop = FALSE]
    rownames(out) <- NULL
  }
  
  out
}

extract_numeric_row <- function(x, row_idx, columns) {
  if (length(columns) == 0) return(numeric(0))
  if (is.null(dim(x)) || row_idx < 1 || row_idx > nrow(x)) {
    return(rep(NA_real_, length(columns)))
  }
  
  if (is.data.frame(x)) {
    out <- vapply(
      x[row_idx, columns, drop = FALSE],
      function(z) suppressWarnings(as.numeric(as.character(z[[1]]))),
      numeric(1)
    )
  } else {
    out <- x[row_idx, columns, drop = TRUE]
    out <- suppressWarnings(as.numeric(as.character(out)))
  }
  
  out[!is.finite(out)] <- NA_real_
  unname(out)
}

get_species_table <- function(obj, species_key) {
  if (is.matrix(obj) || is.data.frame(obj)) return(ok(value = obj))
  
  if (!is.list(obj)) {
    return(fail("History RDS has an unsupported structure."))
  }
  
  candidates <- if (species_key == "H") {
    c("H", "H_contortus", "Haemonchus_contortus", "Haemonchus contortus")
  } else {
    c("T", "T_circumcincta", "Teladorsagia_circumcincta", "Teladorsagia circumcincta")
  }
  
  nms <- names(obj)
  
  if (!is.null(nms)) {
    hit <- candidates[candidates %in% nms]
    if (length(hit) > 0) return(ok(value = obj[[hit[1]]]))
  }
  
  return(fail(
    paste0(
      "Could not find species '", species_key,
      "' in history RDS. Available names: ",
      paste(nms %||% "<unnamed>", collapse = ", ")
    )
  ))
}

normalise_grid_coords <- function(coords) {
  if (inherits(coords, "SpatVector")) {
    xy <- as.data.frame(terra::crds(coords))
    names(xy)[1:2] <- c("x", "y")
    coords <- xy
  }
  
  if (!is.data.frame(coords) && !is.matrix(coords)) {
    return(fail("grid_coords.rds must contain a data.frame, matrix, or SpatVector."))
  }
  
  if (ncol(coords) < 2) {
    return(fail("grid_coords.rds contains fewer than two coordinate columns."))
  }
  
  nms <- tolower(colnames(coords) %||% rep("", ncol(coords)))
  lon_idx <- match(TRUE, nms %in% c("x", "lon", "lng", "longitude"), nomatch = 0)
  lat_idx <- match(TRUE, nms %in% c("y", "lat", "latitude"), nomatch = 0)
  
  if (lon_idx == 0 || lat_idx == 0) {
    lon_idx <- 1
    lat_idx <- 2
  }
  
  lon <- suppressWarnings(as.numeric(as.character(coords[, lon_idx])))
  lat <- suppressWarnings(as.numeric(as.character(coords[, lat_idx])))
  valid <- is.finite(lon) & is.finite(lat)
  
  if (!any(valid)) return(fail("No valid longitude/latitude values in grid_coords.rds."))
  
  ok(value = data.frame(
    source_row = seq_along(lon),
    lon = lon,
    lat = lat,
    valid = valid,
    stringsAsFactors = FALSE
  ))
}

status_box <- function(message) {
  div(
    message,
    style = paste(
      "background:#fff3cd;",
      "border:1px solid #ffe69c;",
      "color:#664d03;",
      "padding:8px 10px;",
      "border-radius:5px;",
      "font-size:0.82em;",
      "margin-bottom:8px;"
    )
  )
}

make_ticks <- function(labels, every = 14L) {
  if (length(labels) == 0) return(list(vals = character(0), text = character(0)))
  idx <- unique(c(seq(1, length(labels), by = every), length(labels)))
  list(vals = labels[idx], text = labels[idx])
}

# ------------------------------------------------------------
# UI (Standard Shiny Layout - No bslib)
# ------------------------------------------------------------

# ------------------------------------------------------------
# UI (Standard Shiny Layout - No bslib)
# ------------------------------------------------------------

ui <- fluidPage(
  tags$head(
    tags$style(HTML(paste0(
      "#map { height: calc(100vh - 160px) !important; min-height: 520px; margin-top: 8px; border-radius: 5px; border: 1px solid #ccc; }",
      ".side-panel { max-height: calc(100vh - 160px); min-height: 520px; overflow-y: auto; margin-top: 8px; }",
      ".plotly.html-widget { min-height:250px !important; }",
      ".js-plotly-plot { min-height:250px !important; }",
      ".custom-card { border: 1px solid #ddd; padding: 15px; border-radius: 5px; background: #fff; margin-bottom: 15px; }",
      ".card-title { margin-top: 0; font-weight: bold; font-size: 1.1em; border-bottom: 1px solid #eee; padding-bottom: 8px; margin-bottom: 15px; }",
      ".info-banner { background: #f8f9fa; border: 1px solid #dee2e6; border-left: 4px solid #1f4e79; border-radius: 5px; padding: 10px 14px; margin-bottom: 10px; font-size: 0.92em; }",
      ".info-top-row { display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between; gap: 8px; }",
      ".info-pills { display: flex; flex-wrap: wrap; gap: 6px; }",
      ".risk-pill { display: inline-block; padding: 2px 8px; border-radius: 12px; color: #fff; font-size: 0.82em; font-weight: 600; }",
      ".pill-low { background: #44ad49; }",
      ".pill-medhigh { background: linear-gradient(90deg, #f5a623, #e65c00); }",
      ".pill-vhigh { background: #8b0000; }",
      ".info-details { margin-top: 6px; }",
      ".info-details summary { display: list-item; cursor: pointer; color: #1f4e79; font-weight: 600; font-size: 0.9em; outline: none; }",
      ".info-details-body { margin-top: 10px; padding-top: 10px; border-top: 1px solid #e2e6ea; }",
      ".info-guide-list { margin: 0 0 10px 0; padding-left: 18px; }",
      ".info-guide-list li { margin-bottom: 4px; }",
      ".info-disclaimer { font-style: italic; color: #495057; background: #fff; border: 1px solid #e2e6ea; border-radius: 4px; padding: 8px 10px; margin-bottom: 8px; font-size: 0.9em; }",
      ".info-citation { font-size: 0.85em; color: #495057; margin: 0; }"
    )))
  ),
  
  titlePanel("Beef + Lamb New Zealand Worm Infection Risk Dashboard"),
  
  tags$div(
    class = "info-banner",
    tags$div(
      class = "info-top-row",
      tags$span(
        "This interactive map displays the weekly risk of pasture contamination by infective ",
        tags$em("Haemonchus contortus"), " and ", tags$em("Teladorsagia circumcincta"),
        " (L3) sheep larvae across New Zealand."
      ),
      tags$div(
        class = "info-pills",
        tags$span("Low (0–25)", class = "risk-pill pill-low"),
        tags$span("Medium–High (26–75)", class = "risk-pill pill-medhigh"),
        tags$span("Very High (76–100)", class = "risk-pill pill-vhigh")
      )
    ),
    tags$details(
      class = "info-details",
      tags$summary("Risk level guide, disclaimer & model citation"),
      tags$div(
        class = "info-details-body",
        tags$ul(
          class = "info-guide-list",
          tags$li(
            tags$strong("Low (0–25) — Baseline monitoring: "),
            "Continue routine flock observation and standard grazing management."
          ),
          tags$li(
            tags$strong("Medium–High (26–75) — Heightened alert: "),
            "Closely monitor weather trends, pasture conditions, and flock health."
          ),
          tags$li(
            tags$strong("Very High (76–100) — Action needed: "),
            "Indicates a peak contamination period; prompt diagnostic checking and strategic grazing management are recommended."
          )
        ),
        tags$div(
          class = "info-disclaimer",
          "Please note: Values represent a relative risk index driven by temperature and rainfall only. They are not absolute worm counts and do not account for farm management, treatments, or stocking density. Always interpret alongside local veterinary advice."
        ),
        tags$p(
          class = "info-citation",
          tags$strong("Model & Citation: "),
          "Risk maps are generated using the GLOWORM-FL mechanistic model, adapted and validated for New Zealand conditions: ",
          "Rose, H., Wang, T., van Dijk, J., & Morgan, E. R. (2015). GLOWORM-FL: a simulation model of the effects of climate and climate change on the free-living stages of gastrointestinal nematode parasites of ruminants. ",
          tags$em("Ecological Modelling"), ", 297, 232–245."
        )
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      class = "side-panel",
      
      selectInput(
        "species_select",
        "Select Species:",
        choices = c(
          "Haemonchus contortus" = "H",
          "Teladorsagia circumcincta" = "T"
        )
      ),
      
      hr(),
      h5("Location Search", style = "font-weight: bold;"),
      textInput("search_loc", NULL, placeholder = "e.g., Hamilton or 3204"),
      actionButton("search_btn", "Search Region/Postcode", class = "btn-primary", style = "width: 100%;"),
      
      hr(),
      h5("Local Grid Analysis", style = "font-weight: bold;"),
      p("Click the map or search above to extract the 0.05° grid time series.", style = "font-size:0.9em;"),
      
      uiOutput("local_stats"),
      br(),
      
      tags$div(
        class = "custom-card",
        tags$div("PASTURE WORM INGESTION RISK", class = "card-title"),
        uiOutput("risk_chart_status"),
        div(
          style = "height:250px; min-height:250px; width:100%; position:relative;",
          plotlyOutput("risk_chart", width = "100%", height = "250px")
        ),
        div(
          "The Risk Score (0-100) measures active ingestion risk based on the density of infective larvae present on the grass.",
          style = "background:#f8f9fa; padding:10px; border-radius:5px; font-size:0.85em; margin-top:10px;"
        )
      ),
      
      tags$div(
        class = "custom-card",
        tags$div("LARVAL RESERVOIR", class = "card-title"),
        uiOutput("reservoir_chart_status"),
        div(
          style = "height:250px; min-height:250px; width:100%; position:relative;",
          plotlyOutput("reservoir_chart", width = "100%", height = "250px")
        ),
        div(
          "Larval Location: 'Soil & Mat' tracks dormant reservoirs, while 'On Grass' tracks the immediate active risk.",
          style = "background:#f8f9fa; padding:10px; border-radius:5px; font-size:0.85em; margin-top:10px;"
        )
      )
    ),
    
    mainPanel(
      width = 8,
      leafletOutput("map")
    )
  )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------

server <- function(input, output, session) {
  
  active_location <- reactiveVal(list(lng = 175.28, lat = -37.78))
  
  # ----------------------------------------------------------
  # Data readers
  # ----------------------------------------------------------
  
  history_raw <- reactivePoll(
    FILE_POLL_MS,
    session,
    checkFunc = function() file_signature(RISK_HISTORY_FILE),
    valueFunc = function() safe_read_rds(RISK_HISTORY_FILE)
  )
  
  reservoir_raw <- reactivePoll(
    FILE_POLL_MS,
    session,
    checkFunc = function() file_signature(RESERVOIR_HISTORY_FILE),
    valueFunc = function() safe_read_rds(RESERVOIR_HISTORY_FILE)
  )
  
  coords_raw <- reactivePoll(
    FILE_POLL_MS,
    session,
    checkFunc = function() file_signature(GRID_COORDS_FILE),
    valueFunc = function() safe_read_rds(GRID_COORDS_FILE)
  )
  
  grid_coords <- reactive({
    x <- coords_raw()
    if (is_fail(x)) return(x)
    normalise_grid_coords(x$value)
  })
  
  selected_raster <- reactive({
    invalidateLater(FILE_POLL_MS, session)
    
    species_key <- input$species_select
    pattern_match <- SPECIES_RASTER_PATTERNS[[species_key]]
    
    validate(need(dir.exists(DATA_DIR), paste("Data directory not found:", DATA_DIR)))
    
    tif_files <- list.files(
      DATA_DIR,
      pattern = paste0(pattern_match, ".*\\.tif$"),
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    validate(need(length(tif_files) > 0, paste("No raster found for", SPECIES_LABELS[[species_key]])))
    
    info <- file.info(tif_files)
    keep <- !is.na(info$mtime)
    validate(need(any(keep), "Raster timestamps could not be read."))
    
    tif_files <- tif_files[keep]
    info <- info[keep, , drop = FALSE]
    latest_file <- tif_files[which.max(info$mtime)]
    
    r <- tryCatch(terra::rast(latest_file), error = function(e) NULL)
    validate(need(!is.null(r), paste("Could not read raster:", basename(latest_file))))
    
    r
  })
  
  # ----------------------------------------------------------
  # Map
  # ----------------------------------------------------------
  
  pal <- colorBin(
    palette = c("#44ad49", "#f5a623", "#e65c00", "#8b0000"),
    domain = c(0, 100),
    bins = c(0, 25, 50, 75, 100),
    na.color = "transparent"
  )
  
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE, minZoom = 5)) %>%
      addProviderTiles(providers$OpenStreetMap) %>%
      setView(lng = 172.5, lat = -41.0, zoom = 5) %>%
      addLegend(
        position = "bottomright",
        colors = c("#44ad49", "#f5a623", "#e65c00", "#8b0000"),
        labels = c("Low", "Medium", "High", "Very High"),
        title = "Risk Level",
        opacity = 1
      )
  })
  
  observe({
    r <- selected_raster()
    
    leafletProxy("map") %>%
      clearImages() %>%
      addRasterImage(
        r,
        colors = pal,
        opacity = 0.8,
        project = TRUE,
        layerId = "risk_layer"
      )
  })
  
  # ----------------------------------------------------------
  # Location selection
  # ----------------------------------------------------------
  
  observeEvent(input$map_click, {
    click <- input$map_click
    if (!is.null(click$lng) && !is.null(click$lat) && is.finite(click$lng) && is.finite(click$lat)) {
      active_location(list(lng = as.numeric(click$lng), lat = as.numeric(click$lat)))
    }
  })
  
  observeEvent(input$search_btn, {
    query <- trimws(input$search_loc %||% "")
    
    if (!nzchar(query)) {
      showNotification("Please enter a New Zealand town, region, or postcode.", type = "warning")
      return()
    }
    
    geo_res <- tryCatch(
      tidygeocoder::geo(
        address = paste(query, "New Zealand"),
        method = "osm",
        quiet = TRUE
      ),
      error = function(e) NULL
    )
    
    if (is.null(geo_res) || nrow(geo_res) == 0) {
      showNotification("Location search failed.", type = "error")
      return()
    }
    
    lat <- suppressWarnings(as.numeric(geo_res$lat))
    lon <- suppressWarnings(as.numeric(geo_res$long))
    valid <- which(is.finite(lat) & is.finite(lon))
    
    if (length(valid) == 0) {
      showNotification("Location not found.", type = "error")
      return()
    }
    
    i <- valid[1]
    active_location(list(lng = lon[i], lat = lat[i]))
  })
  
  selected_grid <- reactive({
    loc <- active_location()
    cobj <- grid_coords()
    
    if (is_fail(cobj)) return(cobj)
    coords <- cobj$value
    
    valid_idx <- which(coords$valid)
    if (length(valid_idx) == 0) return(fail("No valid grid coordinates are available."))
    
    d2 <- (coords$lon[valid_idx] - loc$lng)^2 + (coords$lat[valid_idx] - loc$lat)^2
    j <- which.min(d2)
    row <- valid_idx[j]
    
    ok(
      row_idx = coords$source_row[row],
      grid_lng = coords$lon[row],
      grid_lat = coords$lat[row],
      distance_deg = sqrt(d2[j])
    )
  })
  
  observeEvent(active_location(), {
    loc <- active_location()
    g <- selected_grid()
    
    proxy <- leafletProxy("map") %>% clearGroup("highlight")
    
    if (!is_fail(g) && is.finite(g$distance_deg) && g$distance_deg <= MAX_GRID_DISTANCE_DEG) {
      proxy <- proxy %>% addRectangles(
        lng1 = g$grid_lng - 0.025,
        lat1 = g$grid_lat - 0.025,
        lng2 = g$grid_lng + 0.025,
        lat2 = g$grid_lat + 0.025,
        color = "black",
        weight = 2,
        fill = FALSE,
        group = "highlight"
      )
    }
    
    proxy %>% flyTo(lng = loc$lng, lat = loc$lat, zoom = 10)
  }, ignoreInit = TRUE)
  
  # ----------------------------------------------------------
  # Current raster value
  # ----------------------------------------------------------
  
  output$local_stats <- renderUI({
    loc <- active_location()
    g <- selected_grid()
    
    if (is_fail(g)) {
      return(div(g$message, style = "color:#b02a37; font-weight:600;"))
    }
    
    if (!is.finite(g$distance_deg) || g$distance_deg > MAX_GRID_DISTANCE_DEG) {
      return(div("No model grid cell at this location.", style = "color:#b02a37; font-weight:600;"))
    }
    
    r <- selected_raster()
    
    pt <- terra::vect(
      data.frame(lon = loc$lng, lat = loc$lat),
      geom = c("lon", "lat"),
      crs = "EPSG:4326"
    )
    
    r_crs <- tryCatch(terra::crs(r), error = function(e) "")
    if (!is.null(r_crs) && nzchar(r_crs)) {
      pt <- tryCatch(terra::project(pt, r_crs), error = function(e) pt)
    }
    
    ex <- tryCatch(terra::extract(r, pt), error = function(e) NULL)
    
    if (is.null(ex) || ncol(ex) < 2) {
      return(div("Raster value could not be extracted.", style = "color:#b02a37; font-weight:600;"))
    }
    
    value <- suppressWarnings(as.numeric(ex[1, 2]))
    
    if (!is.finite(value)) {
      return(div("No data (Ocean / outside model mask)", style = "color:#b02a37; font-weight:600;"))
    }
    
    risk_label <- if (value <= 25) "Low" else if (value <= 50) "Medium" else if (value <= 75) "High" else "Very High"
    
    tagList(
      h5(
        paste0("Current Risk Score: ", round(value, 1), " (", risk_label, ")"),
        style = "color:#e65c00; font-weight:bold; margin-bottom:2px;"
      ),
      div(
        paste0("Selected grid: ", round(g$grid_lat, 3), ", ", round(g$grid_lng, 3)),
        style = "font-size:0.8em; color:#6c757d;"
      )
    )
  })
  
  # ----------------------------------------------------------
  # Build risk-series data ONCE, outside Plotly
  # ----------------------------------------------------------
  
  risk_series <- reactive({
    g <- selected_grid()
    if (is_fail(g)) return(fail(g$message))
    
    if (!is.finite(g$distance_deg) || g$distance_deg > MAX_GRID_DISTANCE_DEG) {
      return(fail("No historical grid cell at this location."))
    }
    
    raw <- history_raw()
    if (is_fail(raw)) return(fail(raw$message))
    
    sp <- get_species_table(raw$value, input$species_select)
    if (is_fail(sp)) return(sp)
    
    tab <- sp$value
    if (is.null(dim(tab))) return(fail("Risk history is not a 2-D table."))
    
    if (g$row_idx > nrow(tab)) {
      return(fail(paste0("Grid row ", g$row_idx, " is outside risk_history (", nrow(tab), " rows).")))
    }
    
    dm <- date_column_map(tab)
    if (nrow(dm) == 0) {
      examples <- head(colnames(tab) %||% character(0), 8)
      return(fail(paste0("No date columns detected. Example columns: ", paste(examples, collapse = ", "))))
    }
    
    vals <- extract_numeric_row(tab, g$row_idx, dm$column)
    
    df <- data.frame(
      Date = dm$date,
      DateLabel = format(dm$date, "%Y-%m-%d"),
      Risk = vals,
      stringsAsFactors = FALSE
    )
    
    df <- df[order(df$Date), , drop = FALSE]
    
    df <- tail(df, 15)
    
    if (!any(is.finite(df$Risk))) {
      return(fail("Date columns were found, but this grid row contains no numeric historical risk values."))
    }
    
    ok(data = df)
  })
  
  output$risk_chart_status <- renderUI({
    x <- risk_series()
    if (is_fail(x)) status_box(x$message) else NULL
  })
  
  output$risk_chart <- renderPlotly({
    x <- risk_series()
    req(!is_fail(x))
    
    df <- x$data
    ticks <- make_ticks(df$DateLabel, 1L)
    short_dates <- format(as.Date(ticks$vals), "%d %b")
    
    bar_colors <- sapply(df$Risk, function(r) {
      if (is.na(r)) return("#cccccc")
      if (r <= 25) return("#44ad49")
      if (r <= 50) return("#f5a623")
      if (r <= 75) return("#e65c00")
      return("#8b0000")
    })
    
    risk_labels <- sapply(df$Risk, function(r) {
      if (is.na(r)) return("Unknown")
      if (r <= 25) return("Low")
      if (r <= 50) return("Medium")
      if (r <= 75) return("High")
      return("Very High")
    })
    
    plot_ly(
      x = df$DateLabel,
      y = df$Risk,
      type = "bar",
      marker = list(color = bar_colors),
      text = paste0(format(df$Date, "%d %b %Y"), "<br>Risk: ", round(df$Risk, 1), " (", risk_labels, ")"),
      hoverinfo = "text",
      textposition = "none"
    ) %>%
      layout(
        autosize = TRUE,
        xaxis = list(
          title = "",
          type = "category",
          categoryorder = "array",
          categoryarray = df$DateLabel,
          tickmode = "array",
          tickvals = ticks$vals,
          ticktext = short_dates,
          tickangle = -45,
          tickfont = list(size = 10),
          automargin = TRUE
        ),
        yaxis = list(
          title = "Risk Score",
          range = c(0, 100),
          automargin = TRUE
        ),
        margin = list(l = 50, r = 15, t = 10, b = 60),
        plot_bgcolor = "rgba(0,0,0,0)",
        paper_bgcolor = "rgba(0,0,0,0)"
      ) %>%
      config(displayModeBar = FALSE, responsive = TRUE)
  })
  
  # ----------------------------------------------------------
  # Build reservoir-series data ONCE, outside Plotly
  # ----------------------------------------------------------
  
  reservoir_series <- reactive({
    g <- selected_grid()
    if (is_fail(g)) return(fail(g$message))
    
    risk_x <- history_raw()
    res_x <- reservoir_raw()
    
    if (is_fail(risk_x)) return(fail(risk_x$message))
    if (is_fail(res_x)) return(fail(res_x$message))
    
    risk_sp <- get_species_table(risk_x$value, input$species_select)
    res_sp <- get_species_table(res_x$value, input$species_select)
    
    if (is_fail(risk_sp)) return(risk_sp)
    if (is_fail(res_sp)) return(res_sp)
    
    risk_tab <- risk_sp$value
    res_tab <- res_sp$value
    
    if (is.null(dim(risk_tab)) || is.null(dim(res_tab))) {
      return(fail("One of the history objects is not a 2-D table."))
    }
    
    if (g$row_idx > nrow(risk_tab) || g$row_idx > nrow(res_tab)) {
      return(fail("Selected grid row is outside one of the history databases."))
    }
    
    rm <- date_column_map(risk_tab)
    sm <- date_column_map(res_tab)
    
    if (nrow(rm) == 0 || nrow(sm) == 0) {
      return(fail("No valid date columns were found in one of the history databases."))
    }
    
    rm$key <- as.character(rm$date)
    sm$key <- as.character(sm$date)
    
    common <- merge(
      rm[, c("key", "column")],
      sm[, c("key", "column")],
      by = "key",
      suffixes = c("_risk", "_reservoir"),
      all = FALSE,
      sort = TRUE
    )
    
    if (nrow(common) == 0) return(fail("Risk and reservoir histories contain no matching dates."))
    
    common$Date <- as.Date(common$key)
    common <- common[order(common$Date), , drop = FALSE]
    
    grass <- extract_numeric_row(risk_tab, g$row_idx, common$column_risk)
    reservoir <- extract_numeric_row(res_tab, g$row_idx, common$column_reservoir)
    
    grass[is.finite(grass) & grass < 0] <- NA_real_
    reservoir[is.finite(reservoir) & reservoir < 0] <- NA_real_
    
    if (RESERVOIR_HISTORY_IS_TOTAL) {
      soil <- pmax(0, reservoir - grass)
    } else {
      soil <- reservoir
    }
    
    total <- grass + soil
    valid <- is.finite(total) & total > 0
    
    gp <- rep(NA_real_, length(total))
    sp <- rep(NA_real_, length(total))
    
    gp[valid] <- grass[valid] / total[valid] * 100
    sp[valid] <- soil[valid] / total[valid] * 100
    
    df <- data.frame(
      Date = common$Date,
      DateLabel = format(common$Date, "%Y-%m-%d"),
      Grass = gp,
      Soil = sp,
      stringsAsFactors = FALSE
    )
    df <- tail(df, 15)
    if (!any(is.finite(df$Grass) | is.finite(df$Soil))) {
      return(fail("Matching dates exist, but this grid row contains no usable reservoir values."))
    }
    
    ok(data = df)
  })
  
  output$reservoir_chart_status <- renderUI({
    x <- reservoir_series()
    if (is_fail(x)) status_box(x$message) else NULL
  })
  
  output$reservoir_chart <- renderPlotly({
    x <- reservoir_series()
    req(!is_fail(x))
    
    df <- x$data
    ticks <- make_ticks(df$DateLabel, 1L)
    short_dates <- format(as.Date(ticks$vals), "%d %b")
    
    plot_ly(x = df$DateLabel) %>%
      add_bars(
        y = df$Soil,
        name = "Soil & Mat (Dormant)",
        marker = list(color = "#8d6e63"),
        text = paste0(format(df$Date, "%d %b %Y"), "<br>Soil & Mat: ", round(df$Soil, 1), "%"),
        hoverinfo = "text",
        textposition = "none"
      ) %>%
      add_bars(
        y = df$Grass,
        name = "On Grass (Infective)",
        marker = list(color = "#4caf50"),
        text = paste0(format(df$Date, "%d %b %Y"), "<br>On Grass: ", round(df$Grass, 1), "%"),
        hoverinfo = "text",
        textposition = "none"
      ) %>%
      layout(
        autosize = TRUE,
        barmode = "stack",
        xaxis = list(
          title = "",
          type = "category",
          categoryorder = "array",
          categoryarray = df$DateLabel,
          tickmode = "array",
          tickvals = ticks$vals,
          ticktext = short_dates,
          tickangle = -45,
          tickfont = list(size = 10),
          automargin = TRUE
        ),
        yaxis = list(
          title = "Proportion (%)",
          range = c(0, 100),
          automargin = TRUE
        ),
        legend = list(
          orientation = "h",
          x = 0.5,
          xanchor = "center",
          y = -0.5,
          yanchor = "top"
        ),
        margin = list(l = 50, r = 15, t = 10, b = 120),
        plot_bgcolor = "rgba(0,0,0,0)",
        paper_bgcolor = "rgba(0,0,0,0)"
      ) %>%
      config(displayModeBar = FALSE, responsive = TRUE)
  })
}

shinyApp(ui, server)
