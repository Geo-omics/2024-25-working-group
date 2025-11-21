#!/usr/bin/env Rscript
renv::activate()

'Extract Cyanobacteria Index (CI) Values Predicted by HABTracker

Usage:
  extract_habtracker_ci.R -l LATITUDE -L LONGITUDE -d DATE [options]
  extract_habtracker_ci.R -c CSV [options]
  extract_habtracker_ci.R --csv CSV [options]

Options:
  -l LATITUDE --lat=LATITUDE       Latitude in decimal degrees.
  -L LONGITUDE --lon=LONGITUDE     Longitude in decimal degrees.
  -d DATE --date=DATE              Date of sample collection (YYYY-MM-DD).
  -t TIME --time=TIME              Time of sample collection [default: 12:00].
  -r DISTANCE --rad=DISTANCE       Radial distance (m) [default: nearest CI value].
  -c CSV --csv=CSV                 CSV input file path (for batch mode).
  -o OUTPUT --out=OUTPUT           Output file path.
  -h --help' -> doc  

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
library(ncdf4) #new
library(fields) #new
library(sf)
library(png)
library(readr)
source("code/hab_tracker_2017_code_share/hab_tracker/habtracker_functions.R")

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

## for testing interactively
#arguments <- docopt(doc, args = c(" -l 41.8077 -L -83.3435 -d 2021-08-19 -t 14:00 -r 1000 -o /geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2024/MODIS_acquisition/worker/test.csv"))
#arguments <- docopt(doc, args = c(" -l 41.8077 -L -83.3435 -d 2021-08-19 -t 14:00 -r 1000 -csv /geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/MODIS_acquisition/test.csv"))
arguments <- docopt(doc, args = c("--csv /geomicro/data2/pdenuyl2/CIGLR_bioinformatics/2025/MODIS_acquisition/test_nat.csv"))

#print(arguments)

###########################
### Script starts here ####
###########################

#Get list of tif files and their directories
tiff_paths <- system("ls LErie_HABTracker_archive_FromYizhen/*/gtiff/*.tif",intern = TRUE) %>% 
  tibble(tif_path = .) %>% 
  mutate(model_name = tif_path %>% str_remove(".*HABTracker_archive_FromYizhen/") %>% str_remove("/gtiff.*"),
         model_year = model_name %>% str_sub(1, 4),
         model_month = model_name %>% str_sub(5, 6),
         model_day = model_name %>% str_sub(7, 8),
         model_hour = model_name %>% str_sub(9, 10),
         model_time_edt = ymd_h(paste(model_year, model_month, model_day, model_hour), tz = "America/New_York"),
         forecast_name = tif_path %>% str_remove(".*gtiff/") %>% str_remove(".tif"),
         forecast_year = forecast_name %>% str_sub(1, 4),
         forecast_month = forecast_name %>% str_sub(5, 6),
         forecast_day = forecast_name %>% str_sub(7, 8),
         forecast_hour = forecast_name %>% str_sub(9, 10),
         forecast_time_utc = ymd_h(paste(forecast_year, forecast_month, forecast_day, forecast_hour)),
         model_forecast_diff = abs(forecast_time_utc - model_time_edt)) %>%
  group_by(forecast_time_utc) %>%
  slice_min(model_forecast_diff) %>% #select the smallest time difference between model run and forecast time for each forecast time
  ungroup()

######################################
#SAMPLE STREAM FOR BULK SAMPLES! 
######################################

results_list <- list()
all_final_results <- tibble()

if (!is.null(arguments$csv)) {
#Define features
  input <- read_csv(arguments$csv, show_col_types=FALSE)

for (i in seq_len(nrow(input))) {

#clean up time (if NA, adjust) - 12:00, rad = NA
collection_lat <- as.numeric(input$lat[i])   #latitude
collection_lon <- as.numeric(input$lon[i]) #longitude
collection_date <- ymd(input$date[i])
collection_time_prep <- hms(ifelse(!is.na(input$time[i]), as.character(input$time[i]), "12:00:00"))
collection_time <- sprintf("%02d:%02d", hour(collection_time_prep), minute(collection_time_prep))
rad <- ifelse(!is.na(input$rad[i]), input$rad[i], NA)
output_file <- arguments$o  #output path    

# Find nearest forecast_time using tidyverse
sample_time <- paste0(collection_date, "T", collection_time) %>% ymd_hm(tz = "America/New_York")   #Combine sample information
# Convert EDT time to UTC
sample_time_utc <- with_tz(sample_time, tz = "UTC")

nearest_time_row <- tiff_paths %>%
  mutate(sample_forecast_diff = abs(forecast_time_utc - sample_time_utc)) %>%
  dplyr::slice_min(sample_forecast_diff)

tif_image_path <- nearest_time_row$tif_path

#Next section of script can be divided into two parts:
#
#part 1: 
#Mark Rowe 7-2-2015
#this R script will read a geotiff CI file from Stumpf/Wynne
#and write a CI and chlorophyll *.txt file
#
#part 2:
#take the output from part 1 and compare coordinates from samples of interest
#the output will provide distance between points of ci/chlorophyll readings and sample

########
#Part 1#
########

#open the geotiff
r1 <- raster(tif_image_path)

# Extract the CRS from 'r1'
crs <- crs(r1)
# Convert the CRS to PROJ string format
prj1 <- wkt(crs)

##THIS IS WHERE THRE"S A PROBLEM
dfrc_raw <- as.data.frame(r1, xy=TRUE)
names(dfrc_raw) <- c("X","Y","ci")
dfrc_process <- projLatLon(dfrc_raw/1000, prj1)
dfrc_process$ci <- dfrc_process$ci*1000
dfrc_process2 <- dfrc_process[, -c(1:2)] #delete columns 1 through 2 (X and Y)

dfrc <- dfrc_process2[dfrc_process2$ci <= 250,] # remove land and cloud pixels

dfrc$chl <- 0
dfrc$chl[dfrc$ci == 0] <- 0
xx <- which(dfrc$ci %in% 1:250)
CI <- 10^( dfrc$ci[xx]/100 -4)
dfrc$chl[xx] <- 12570*CI + 10

#View(dfrc)

########
#Part 2#
########

setDT(dfrc)

all_coord <- matrix(c(dfrc$lon, dfrc$lat), ncol = 2)

dfci_results <- distHaversine(all_coord, c(collection_lon,collection_lat))

dfci_results_df <- setNames(data.frame(dfci_results), "dist_m")
dfci_final <- cbind(dfrc, dfci_results_df)

if(!is.na(rad)) {
       result <- as.data.frame(dfci_final %>%
         filter(dist_m <= as.numeric(rad)) %>%
         mutate(ci_mean = mean(ci),
                chl_mean = mean(chl),
                model_image = tif_image_path,
                measurement_lat = lat,
                measurement_lon = lon,
                collection_lat = collection_lat,
                collection_lon = collection_lon,
                radius_limit_m = as.numeric(rad)) %>%
         dplyr::select(collection_lat, collection_lon, measurement_lat, measurement_lon, ci, chl, dist_m, radius_limit_m, ci_mean, chl_mean, model_image) %>%
         dplyr::arrange(dist_m))
 } else {
  #get nearest calculated values
  result <- as.data.frame(dfci_final %>% 
    slice_min(dist_m) %>%
    dplyr::select(ci, chl, dist_m))}
  results_list[[i]] <- result
}}
all_final_results <- dplyr::bind_rows(results_list) 

# Write to output if requested, else print
if (!is.null(arguments$o)) {
  write.csv(all_final_results, arguments$o, row.names = FALSE)
} else {
  print(all_final_results)} 

if (is.null(arguments$csv)) {
  ######################################
  #SAMPLE STREAM FOR INDIVIDUAL SAMPLES! 
  ######################################
  #Define features
  #Coordinates - Format: decimal
  collection_lat <- as.numeric(arguments$lat)   #latitude
  collection_lon <- as.numeric(arguments$lon) #longitude
  collection_date <- arguments$d   #date of sample collection
  collection_time <- ifelse(!is.null(arguments$t), arguments$t, "12:00")    #time of sample collection (tz = EDT)
  output_file <- arguments$o  #output path
  
  # Find nearest forecast_time using tidyverse
  sample_time <- paste0(collection_date, "T", collection_time) %>% ymd_hm(tz = "America/New_York")   #Combine sample information
  # Convert EDT time to UTC
  sample_time_utc <- with_tz(sample_time, tz = "UTC")
  
  nearest_time_row <- tiff_paths %>%
    mutate(sample_forecast_diff = abs(forecast_time_utc - sample_time_utc)) %>%
    slice_min(sample_forecast_diff)
  
  tif_image_path <- nearest_time_row$tif_path
  
  #Next section of script can be divided into two parts:
  #
  #part 1: 
  #Mark Rowe 7-2-2015
  #this R script will read a geotiff CI file from Stumpf/Wynne
  #and write a CI and chlorophyll *.txt file
  #
  #part 2:
  #take the output from part 1 and compare coordinates from samples of interest
  #the output will provide distance between points of ci/chlorophyll readings and sample
  
  ########
  #Part 1#
  ########
  
  #open the geotiff
  r1 <- raster(tif_image_path)
  
  # Extract the CRS from 'r1'
  crs <- crs(r1)
  # Convert the CRS to PROJ string format
  prj1 <- wkt(crs)
  
  dfrc <- as.data.frame(r1, xy=TRUE)
  names(dfrc) <- c("X","Y","ci")
  dfrc <- projLatLon(dfrc/1000, prj1)
  dfrc$ci <- dfrc$ci*1000
  dfrc <- dfrc[, -c(1:2)] #delete columns 1 through 2 (X and Y)
  
  dfrc <- dfrc[dfrc$ci <= 250,] # remove land and cloud pixels
  
  dfrc$chl <- 0
  dfrc$chl[dfrc$ci == 0] <- 0
  xx <- which(dfrc$ci %in% 1:250)
  CI <- 10^( dfrc$ci[xx]/100 -4)
  dfrc$chl[xx] <- 12570*CI + 10
  
  #View(dfrc)
  
  ########
  #Part 2#
  ########
  
  setDT(dfrc)
  
  all_coord <- matrix(c(dfrc$lon, dfrc$lat), ncol = 2)
  
  dfci_results <- distHaversine(all_coord, c(collection_lon,collection_lat))
  
  dfci_results_df <- setNames(data.frame(dfci_results), "dist_m")
  dfci_final <- cbind(dfrc, dfci_results_df)
  #View(dfci_final)
  
  #export explained:
  #ci is cyanobacterial index
  #lon is longitude of measurement point of MODIS image
  #lat is latitude of measurement point of MODIS image
  #chl is calculated chlorophyll
  
  #remove any values further than 1000m (1km)
  #dfci_final_rm_1km <- dfci_final[!dfci_final[["dist_m"]] >= 1000, ]
  #ci_mean_1km <- mean(dfci_final_rm_1km$ci)
  #chl_mean_1km <- mean(dfci_final_rm_1km$chl)
  
  #remove any values further than 5000m (5km)
  #dfci_final_rm_5km <- dfci_final[!dfci_final[["dist_m"]] >= 5000, ]
  #ci_mean_5km <- mean(dfci_final_rm_5km$ci)
  #chl_mean_5km <- mean(dfci_final_rm_5km$chl)
  
  if (!is.null(rad) & !is.na(rad)) {
    result <- as.data.frame(dfci_final %>%
                              filter(dist_m <= as.numeric(rad)) %>%
                              mutate(
                                ci_mean = mean(ci),
                                chl_mean = mean(chl),
                                model_image = tif_image_path,
                                # define these if you have them:
                                # measurement_lat = lat,
                                # measurement_lon = lon,
                                collection_lat = collection_lat,
                                collection_lon = collection_lon,
                                radius_limit_m = rad
                              ) %>%
                              dplyr::select(
                                collection_lat, collection_lon, 
                                # measurement_lat, measurement_lon, # comment/uncomment as needed
                                ci, chl, dist_m, radius_limit_m, ci_mean, chl_mean, model_image
                              ) %>%
                              dplyr::arrange(dist_m))
  } else {
    result <- as.data.frame(dfci_final %>%
                              slice_min(dist_m) %>%
                              dplyr::select(ci, chl, dist_m))
  }
  
  if(!is.null(arguments$o)) {
    write.csv(result, output_file, row.names = FALSE)
  } else {
    print(result)}
}




  
  
  
  
  