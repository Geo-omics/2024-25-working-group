# Mark Rowe 7-2-2015
# This R script will read a geotiff CI file from Stumpf/Wynne
# and write a chlorophyll *.txt file for input to habtracker
# This script is meant to be called by habtracker.R

setwd("Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/")

# run the next two lines to load up variables to run this script outside of habtracker.R
 #dryrun <- TRUE
 #source("habtracker.R")

print(" load common data and functions:  habtracker_functions.R")
source("habtracker_functions.R")

print("make standard format date strings from the Wynne Stumpf file names")

# Aug. 31 Aqua
string1 <- strsplit("CI Processed/terra/Sept2.terra.2018245.0902.1602_1605_1740C.L3.GL1.v930044d4967_1_1.CI_merge_Kdcorr-E.tiff", "/")[[1]]
string1 <- string1[length(string1)] 
string2 <- strsplit(string1, "\\.")[[1]]
hr1 <- substr(string2[4],1,4)
if(is.na(as.numeric(hr1))){ # if hour is not provided in the filename
  hr1 <- "1700" 
}
datestring <- paste0(string2[2],hr1,"00")
year <- substr(string2[2], 1, 4)

imageDate <- as.POSIXct(datestring, format="%Y%j%H%M", tz="UTC")

datestring1 <- substr(datestring, 1, 7)

  fntxt <- paste0("/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/CI Processed/", format(imageDate, "M%Y%j%H%M00.txt") )

  # open the geotiff
  r1 <- raster("CI Processed/terra/Sept2.terra.2018245.0902.1602_1605_1740C.L3.GL1.v930044d4967_1_1.CI_merge_Kdcorr-E.tiff")
  prj1 <- CRS(r1@crs@projargs)
  if(substr(prj1@projargs, 1, 12) == " +proj=sinu "){
  # I adjusted +x and +y to get it to align with shorelineLL, it was shifted north by 18.5 km or so
#  prj1 <- CRS("+proj=sinu +lon_0=-82 +x_0=1000 +y_0=18500 +ellps=WGS84 +units=m +no_defs")
  }
  
  # crop it to Lake Erie erie extent
  nodes <- data.frame(
    lon=c(-83.50791, -78.83501), 
    lat=c(41.35162, 42.93132)
    )
  nodes <- projxy(nodes, prj1)
  xlim <- nodes$X
  ylim <- nodes$Y
  
 #  plot(r1, xlim=xlim, ylim=ylim)
  
  e <- extent(min(xlim), max(xlim), min(ylim), max(ylim))
  rc <- crop(r1, e)
  #plot(rc)
  
  dfrc <- as.data.frame(rc, xy=TRUE)
  names(dfrc) <- c("X","Y","ci")
  dfrc <- projLatLon(dfrc/1000, prj1)
  dfrc$ci <- dfrc$ci*1000
  
  dfrc <- dfrc[dfrc$ci <= 250,] # remove land and cloud pixels
  
View(dfrc)
write.csv(dfrc, file = "dfci.csv") 

  ###  convert CI to chlorophyll 
#   Rick sent me the following explanation for interpretation of the pixel values on the MODIS images. 
#   Do the same conventions apply to the images that you sent?
#   
#   0 is no data, >250 is other stuff,  clouds, land. 
#   1-249 is CI
#   
#   CI = 10^[DN/100 - 4]
#   
#   where DN is the pixel value (1-249) ?
#   
#   chl nominally = 12,570*CI + 10
#   
#   so in my example the correct answer is for
#   
#   a DN of 100 = Ci of 0.001 the chl of 23 ug/L
#     
#   Timothy Wynne - NOAA Federal <timothy.wynne@noaa.gov>
#   yes.
#   it should be the same.
#   252 is land
#   253 is clouds  
#  250 appear to be scum areas (pink)
  
#   cival = 250
#   col1 = rc@legend@colortable[cival+1]
#   xx <- which(dfrc$ci == cival)
#   plot(dfrc$lon[xx], dfrc$lat[xx], pch=15, xlim=xlimLL, ylim=ylimLL, col=col1)
  
#   cols <- fields::color.scale( as.numeric(dfrc$ci), col=rc@legend@colortable)
#   plot(dfrc$lon, dfrc$lat, pch=15, xlim=xlimLL, ylim=ylimLL, col=rc@legend@colortable[dfrc$ci])
  
  #dfrc <- dfrc[dfrc$ci != 252,] # remove land pixels
  dfrc <- dfrc[dfrc$ci <= 250,] # remove land and cloud pixels
  dfrc$chl <- 0
  dfrc$chl[dfrc$ci == 0] <- 0
  xx <- which(dfrc$ci %in% 1:250)
  CI <- 10^( dfrc$ci[xx]/100 -4)
  dfrc$chl[xx] <- 12570*CI + 10
  dfrc$temp <- NA
  
#   plot(rc)
#   par(def.par)
#    plotPolys(shorelineLL, xlim=xlimLL, ylim=ylimLL)
#    xx <- which(dfrc$ci == 0)
#    points(dfrc$lon, dfrc$lat, pch=15, cex=0.5)
  
  # write the txt file in the format that habtracker reads
   ptsmat <- as.matrix(dfrc[,c("lon","lat","ci","chl","temp")])
#  ptsmat <- sprintf(fmt="%12.6f %12.6f %10.6f %10.6f", ptsmat[,1], ptsmat[,2] , ptsmat[,3], ptsmat[,4] )
#  write.table(ptsmat, file = fntxt, append = FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)

  write.csv(ptsmat, file = "Sept2_terra_dfci_chl.csv")  
# # create a uniform concentration image file
# dfrc <- data.frame(
#   lon = coordsn$lon
#   ,lat = coordsn$lat
#   ,chl = 19
#   ,temp = NA
#   )

