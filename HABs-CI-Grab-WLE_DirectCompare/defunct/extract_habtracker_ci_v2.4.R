#!/usr/bin/env Rscript
renv::activate()

'Extract Cyanobacteria Index (CI) Values Predicted by HABTracker

Usage:
  extract_habtracker_ci.R --lat latitude --lon longitude -d date [options]

Options:
  -l, --lat=LATITUDE       Latitude in decimal degrees.
  -L, --lon=LONGITUDE      Longitude in decimal degrees.
  -d, --date=DATE          Date of sample collection (format: YYYY-MM-DD).
  -t, --time=TIME          Time of sample collection in EDT (format: hh:mm). DEFAULT = 12:00.
  -r, --rad=DISTANCE       Radial distance in meters to calculate average CI values. DEFAULT = nearest CI value.
  -i, --csv=CSV            CSV containing lat/lon/date/time/rad for multiple samples, by row.
  -o, --out=OUTPUT         Output file.
  -h, --help               Show this help screen.
' -> doc

###########################
### Dependencies ####
###########################
library(docopt)
library(PBSmapping)
library(geosphere)
library(data.table)
library(dplyr)
library(stringr)
library(lubridate)
library(raster)
library(ncdf4)
library(fields)
library(sf)
library(png)
library(readr)

source("code/hab_tracker_2017_code_share/hab_tracker/habtracker_functions.R")

# Parse arguments
arguments <- docopt(doc)
print(arguments)

# Core function: does everything for one sample
process_sample <- function(lat, lon, collection_date, collection_time = "12:00", radial = NULL) {
  ## Find appropriate tif file
  tiff_paths <- system("ls LErie_HABTracker_archive_FromYizhen/*/gtiff/*.tif",intern = TRUE) %>% 
    tibble(tif_path = .) %>%
    mutate(
      model_name = str_remove(tif_path, ".*HABTracker_archive_FromYizhen/") %>% str_remove("/gtiff.*"),
      model_year = str_sub(model_name, 1, 4),
      model_month = str_sub(model_name, 5, 6),
      model_day = str_sub(model_name, 7, 8),
      model_hour = str_sub(model_name, 9, 10),
      model_time_edt = ymd_h(paste(model_year, model_month, model_day, model_hour), tz = "America/New_York"),
      forecast_name = str_remove(str_remove(tif_path, ".*gtiff/"), ".tif"),
      forecast_year = str_sub(forecast_name, 1, 4),
      forecast_month = str_sub(forecast_name, 5, 6),
      forecast_day = str_sub(forecast_name, 7, 8),
      forecast_hour = str_sub(forecast_name, 9, 10),
      forecast_time_utc = ymd_h(paste(forecast_year, forecast_month, forecast_day, forecast_hour)),
      model_forecast_diff = abs(forecast_time_utc - model_time_edt)) %>%
    group_by(forecast_time_utc) %>%
    slice_min(model_forecast_diff) %>%
    ungroup()
  
  # Find nearest forecast time
  sample_time <- paste0(collection_date, "T", collection_time) %>% ymd_hm(tz = "America/New_York")
  sample_time_utc <- with_tz(sample_time, tz = "UTC")
  
  nearest_time_row <- tiff_paths %>%
    mutate(sample_forecast_diff = abs(forecast_time_utc - sample_time_utc)) %>%
    slice_min(sample_forecast_diff)
  
  tif_image_path <- nearest_time_row$tif_path[1]
  
  # Open the geotiff
  r1 <- raster(tif_image_path)
  
  # CRS, projection, etc.
  crs <- crs(r1)
  prj1 <- wkt(crs)
  
  dfrc <- as.data.frame(r1, xy=TRUE)
  names(dfrc) <- c("X","Y","ci")
  dfrc <- projLatLon(dfrc/1000, prj1)
  dfrc$ci <- dfrc$ci*1000
  dfrc <- dfrc[, -c(1:2)] # remove X and Y
  dfrc <- dfrc[dfrc$ci <= 250,]
  dfrc$chl <- 0
  dfrc$chl[dfrc$ci == 0] <- 0
  xx <- which(dfrc$ci %in% 1:250)
  CI <- 10^(dfrc$ci[xx]/100 - 4)
  dfrc$chl[xx] <- 12570*CI + 10
  
  setDT(dfrc)
  
  all_coord <- matrix(c(dfrc$lon, dfrc$lat), ncol = 2)
  dfci_results <- distHaversine(all_coord, c(lon,lat))
  dfci_results_df <- setNames(data.frame(dfci_results), "dist_m")
  dfci_final <- cbind(dfrc, dfci_results_df)
  
  # -------- output processing
  if (!is.null(radial)) {
    out <- as.data.frame(dfci_final %>%
                           filter(dist_m <= as.numeric(radial)) %>%
                           mutate(ci_mean = mean(ci),
                                  chl_mean = mean(chl)) %>%
                           dplyr::select(ci, chl, dist_m, ci_mean, chl_mean) %>%
                           dplyr::arrange(dist_m))
  } else {
    out <- as.data.frame(dfci_final %>%
                           slice_min(dist_m) %>%
                           dplyr::select(ci, chl, dist_m))
  }
  return(out)
}

# -------------- BATCH OR SINGLE SAMPLE ENTRY POINT --------------

if (!is.null(arguments$csv)) {
  input <- read_csv(arguments$csv, show_col_types=FALSE)
  if(!"time" %in% names(input)) input$time <- "12:00"
  if(!"rad" %in% names(input)) input$rad <- NA
  results_list <- vector("list", nrow(input))
  for (i in seq_len(nrow(input))) {
    lat <- input$lat[i]
    lon <- input$lon[i]
    date <- input$date[i]
    time <- input$time[i]
    rad <- input$rad[i]
    cat(sprintf("Processing sample %d: %f, %f, %s, %s, %s\n", i, lat, lon, date, time, rad))
    result <- tryCatch(
      process_sample(lat, lon, date, time, if(is.na(rad)) NULL else rad),
      error = function(e) data.frame(error=as.character(e))
    )
    result$lat <- lat
    result$lon <- lon
    result$date <- date
    result$time <- time
    result$rad <- rad
    results_list[[i]] <- result
  }
  final_results <- dplyr::bind_rows(results_list)
  if (!is.null(arguments$o)) {
    write.csv(final_results, arguments$o, row.names=FALSE)
  } else {
    print(final_results)
  }
  quit(save = "no")
}

# SINGLE SAMPLE MODE
lat <- as.numeric(arguments$lat)
lon <- as.numeric(arguments$lon)
collection_date <- arguments$d
collection_time <- ifelse(!is.null(arguments$t), arguments$t, "12:00")
output_file <- arguments$o
radial <- arguments$`--rad`

res <- process_sample(lat, lon, collection_date, collection_time, radial)
if (!is.null(output_file)) {
  write.csv(res, output_file, row.names=FALSE)
} else {
  print(res)
}
       