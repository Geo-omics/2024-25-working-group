# Mark Rowe 6-16-2015
# load data and functions that are used in habtracker pre and post processing
#install.packages("PBSmapping")
require(PBSmapping)
require(ncdf4)
require(fields)
require(raster)

# assign values here only if not calling from habracker.R script 
if(exists("xlimLL")==FALSE){
  xlimLL <- c(-83.47336, -82.01775)
  ylimLL <- c(41.36373, 42.08565)
}

# define functions to project coordinates
source("code/hab_tracker_2017_code_share/hab_tracker/proj_functions.R")



# load FVCOM grid and polygons
load(file="code/hab_tracker_2017_code_share/hab_tracker/erie05lcc_grid_coords.Rdata")
load(file="code/hab_tracker_2017_code_share/hab_tracker/erie05lcc_grid_polys.Rdata")
shoreline$X <- shoreline$X
shoreline$Y <- shoreline$Y
psTces$X <- psTces$X
psTces$Y <- psTces$Y
coordsn$X <- coordsn$X/1000
coordsn$Y <- coordsn$Y/1000
coordse$X <- coordse$X/1000
coordse$Y <- coordse$Y/1000

# make a land mask polygon
mask1 <- shorelineLL[1:4,]
mask1$X <- c(-100,-100,393.328,393.328)
mask1$Y <- c(-100,200.000,200.000,-100)
mask2 <- joinPolys(mask1, shoreline, operation="DIFF")
# plotPolys(mask1, col="tan")
# addPolys(shoreline, col="darkblue")
# addPolys(mask2, col="tan")

zlab <- expression(paste("Chlorophyll-a ",mu,g," ",L^{-1}))

lims <- projLcc(data.frame(lon=xlimLL, lat=ylimLL) , prj_new)
xlimLcc <- lims$X
ylimLcc <- lims$Y


## Function to draw a logo on a plot from a png image
require(png)
require(grid)
logo.draw <- function(z, x, y, size, ...) {
  asp <- dim(z)[1]/dim(z)[2]
  logo <- rasterGrob(image = z,
                     x = unit(x, "npc") , y = unit(y, "npc"),
                     width = unit(size, "cm"), height = unit(size*asp, "cm"),
                     just = c("left", "bottom")
                     , gp = gpar(...))
  grid.draw(logo)
}

# plotPolys(shoreline, xlim=xlimLcc, ylim=ylimLcc)
# ## select the PNG file
# z <- readPNG("noaa_glerl_ciler_logo.png")
# logo.draw( z, 0.1, 0.1, 10.0, fontsize = 8)

# function to add cities to a plot

# a few cities to add to the plot base map
cities <- read.csv("code/hab_tracker_2017_code_share/hab_tracker/cities.csv")
cities <- projLcc(cities, prj_new)
plotCities <- function(cex.cities=0.6){
  # add some cities to the base map
  #pos: 1, 2, 3 and 4,  below, to the left of, above and to the right of the specified coordinates.
  #cex.cities <- 0.6 
  points(cities$X, cities$Y, pch=16)
  xx <- which(cities$name %in% c("Port Clinton","Sandusky","Lorain","Cleveland","Ashtabula","Mentor"))
  text(cities$X[xx], cities$Y[xx], cities$name[xx], pos=1, cex=cex.cities)
  xx <- which(cities$name %in% c("Rockwood"))
  text(cities$X[xx]-8.000, cities$Y[xx], cities$name[xx], pos=1, cex=cex.cities)
  xx <- which(cities$name %in% c("Monroe"))
  text(cities$X[xx]-8.000, cities$Y[xx], cities$name[xx], pos=1, cex=cex.cities)
  xx <- which(cities$name %in% c("Toledo","Windsor","Amherstburg"))
  text(cities$X[xx]+0.000, cities$Y[xx], cities$name[xx], pos=4, cex=cex.cities)
  xx <- which(cities$name %in% c("Leamington","Detroit"))
  text(cities$X[xx], cities$Y[xx], cities$name[xx], pos=2, cex=cex.cities)
  xx <- which(cities$name %in% c("Port Stanley"))
  text(cities$X[xx]-2, cities$Y[xx]+0.5, cities$name[xx], pos=2, cex=cex.cities)
}
# plot(cities$X, cities$Y)
# text(cities$X, cities$Y, cities$name)

# function to draw a scalebar
x <- -86.5
y <- 42
dist <- 10
scalebar <- function(x, y, dist=10, lwd=1, vertical=FALSE, cex=1.0){
  y2 <- y
  require(sp)
  pts <- matrix(c(x, y),nrow=1,ncol=2)
  pt <- matrix(c(x+1, y),nrow=1,ncol=2)
  kmperlon <- 1.0# spDistsN1(pts, pt, longlat=TRUE)
  x2 <- x+dist/kmperlon
  lines(c(x,x2),c(y,y2),lwd=lwd, lend=2)
  if(vertical){
    pt <- matrix(c(x, y+1),nrow=1,ncol=2)
    kmperlat <- 1.0# spDistsN1(pts, pt, longlat=TRUE)
    y2 <- y+dist/kmperlat
    lines(c(x,x),c(y,y2),lwd=lwd, lend=2)
  }
  voffset <- par("cxy")[2]*cex #/2.5
  text((x+x2)/2, y-voffset, paste(dist, "km"), cex=cex)
}


processChlData <- function(chlData){
  require(PBSmapping)
  
  # project modis pixel locations to prj_new
  chlData <- projLcc(chlData, prj_new)
  
  # to speed up model simulations, only initialize particles west of specified X value
  # easternLimit <- xlimLL[2]
  xx <- which(chlData$lon < easternLimit)
  chlData <- chlData[xx,]
  
  # deal with no data values
  chlData$temp[chlData$temp > 40] <- NA
  chlData$chl[chlData$chl == -999] <- NA
  chlData$temp[chlData$temp == -999] <- NA
  # chlData[is.na(chlData$chl+chlData$temp) == TRUE,c("chl","temp")] <- NA # if either chl or temp is missing, set both to missing
  chlData$chl[chlData$chl < 0.0] <- 0.1 # set negative chl to low value
  
  # select satellite values within model domain 
  chlData$EID <- 1:nrow(chlData)
  fp <- findPolys(events=chlData, polys=shoreline, maxRows = 1e+07)
  chlData <- chlData[chlData$EID %in% fp$EID,]
  
  fp <- findPolys(events=chlData, polys=psTces, maxRows = 1e+07)
  chlData <- merge(chlData, fp[,c("EID","PID")])
  # PID is the FVCOM node associated with the satellite obs
  chlData$depth <- coordsn$depth[chlData$PID]
  
  return(chlData)
}

if(exists("buffer_rad")==FALSE){
  buffer_rad <- 1.5
}
# define a function to apply a buffer, so that the same buffer can be used when the skill
# stats are calculated in habplot_ptraj_conc2.R
bufferFun <- function(chlData, buffer_rad, buffer_rad2, fill=FALSE){

    # To avoid artifacts caused by bottom reflectance or other edge effects
    # don't use pixels that are shallower than a specified depth, but only
    # apply this outside of the southern, western basin
    # lat1 <- 41.91
    lat1 <- 41.93906 # Stony point
    # lon1 <- -82.60    # just east of islands
    lon1 <- -82.85499 # just west of the islands
    #  dCritNW <- 5 # Exclude pixels shallower than dCritNW North of lat1 in Western Basin
    #  dCritE <- 6 # Exclude pixels shallower than dCritE East of lon1 (Central and Eastern basins)
    #  xx <- which(chlData$lat > lat1 & chlData$depth < dCritNW)
    #  xx2 <- which(chlData$lon > lon1 & chlData$depth < dCritE)
    #  xx <- unique(c(xx, xx2))
    #  #chlData <- chlData[!(1:nrow(chlData) %in% xx),]
    #  chlData$chl[xx] <- NA
    
    # don't use pixels within a buffer distance from land
    
    # buffer_rad outside of the southern western basin
    if(buffer_rad > 0){
      xx <- which(chlData$lat > lat1 )
      xx2 <- which(chlData$lon > lon1 )
      xx <- unique(c(xx, xx2))
    rad <- buffer_rad # 1.5
    for(ro in xx){
      dist <- sqrt((chlData$X[ro]-shoreline$X)^2 + (chlData$Y[ro]-shoreline$Y)^2)
      if(min(dist, na.rm=TRUE) <= rad){
          chlData$chl[ro] <- NA          
      }
    }
    } # end if buffer_rad
    
    # buffer_rad2 inside of the southern western basin
    if(buffer_rad2 > 0){
      xx <- which(chlData$lat <= lat1 & chlData$lon <= lon1 )
      rad <- buffer_rad2 # 1.5
      for(ro in xx){
        dist <- sqrt((chlData$X[ro]-shoreline$X)^2 + (chlData$Y[ro]-shoreline$Y)^2)
        if(min(dist, na.rm=TRUE) <= rad){
          chlData$chl[ro] <- NA          
        }
      }
    } # end if buffer_rad
    
    if(fill){
      # fill in buffer data points with nearest neighbor
      xx <- which(is.na(chlData$chl ))
      chlData2 <- chlData[!(is.na(chlData$chl )),]
      for(ro in xx){
        dist <- sqrt((chlData$X[ro]-chlData2$X)^2 + (chlData$Y[ro]-chlData2$Y)^2)
        mindist <- min(dist, na.rm=TRUE)
          chlData$chl[ro] <- chlData2[which(dist == mindist)[1],"chl"]     
      }
    }
    # remove NA values  
    chlData <- chlData[is.na(chlData$chl)==FALSE,]

} # end buffer function

# Define a function to make a classified colorbar
getColors <- function(values){
  require(fields)
  # values is a numeric vector 
  chl1 <- values # chlData$chl[is.na(chlData$chl)==FALSE]
  if(length(chl1 > 0)){
    brkm <- ceiling(max(chl1, na.rm=TRUE))+1
    if(brkm <= 100){brkm <- 101}
  }else{
    brkm <- 101
  }
  
#  brk <- c(0,18,seq(20,100,20), brkm)
  brk <- c(0,23,30,seq(40,100,20), brkm)
  label <- zlab
  
  brkLabs <- round(brk, 1)
  values <- cut(chl1, breaks=brk)
  nlevel <- length(brk)-1
  zlim <- c(0,nlevel)
  colorbar <- tim.colors(nlevel) # tim.colors(256)
  cols <- fields::color.scale( as.numeric(values), col=colorbar
                               , zlim =zlim 
                               ,eps= 0.0 , transparent.color="gray"
  )
  return(list(cols=cols, zlim=zlim, label=label, brk=brk, colorbar=colorbar))
}

# define a function to get prediction accuracy stats
# Stumpf email 8-4-2015:  Our CI to cells was adjusted to 10^8 cells/mL per CI (so CI = 0.001 = 10^5)
# , recognizing we did not have more precision in the original tuning. This was validated for multiple
# lakes by Lunetta et al., 2014.   So, 10^5 cells is consistent with WHO, based on cells, and that is
# more appropriate.  The WHO chl was based on an assumption of chlorophyll per cell. Also, between all
# of our papers using 0.001 as the threshold, and our observations that that is a good estimate of 
# where the bloom is a problem, that's the best number. 
# By the way, WHO contact for children is 10 ug/L, so if we use their rule of thumb, that puts us to 25 ug/L chl. 
# The other consideration is the nuisance factor at > 10^5 cells the cells clearly discolor the water, and begin to develop scums, so even when they are relatively non-toxic, they are a nuisance. 
chlCrit <-  23 # 18 # define critical chl concentration for a hab. CI = 0.001 = 22.6 ug/L chl
getPstats <- function(chlData, chlCrit=chlCrit){
  # chlData is a dataframe with columns "chl" and "pred"
  pred <- rep("nodata",nrow(chlData))
  pred[chlData$pred >= chlCrit] <- "hab"
  pred[chlData$pred < chlCrit] <- "nohab"  
  obs <- rep("nodata",nrow(chlData))
  obs[chlData$chl >= chlCrit] <- "hab"
  obs[chlData$chl < chlCrit] <- "nohab"
  nobs <- length(which(obs %in% c("hab","nohab"))) # number of observations
  ngnh <- length(which(pred %in% c("nohab","nodata") & obs == "nohab")) # number of good no hab predictions
  ngh <- length(which(pred == "hab" & obs == "hab")) # number of good hab predictions
  nObsHab <-   length(obs[obs=="hab"])
  nObsNoHab <-  length(obs[obs=="nohab"])
  
  # calculate skill stats for chlorophyll concentration
  #remove NAs
  dat <- data.frame(mod=chlData$pred, meas=chlData$chl)
  dat$mod[is.na(dat$mod)==TRUE] <- 0 # assume cells with zero particles in the habtracker model means zero concentration
  dat$sum <- dat$mod + dat$meas
  dat <- dat[is.na(dat$sum)==FALSE,]
  
  if(nrow(dat) > 0){
    mod <- dat$mod
    meas <- dat$meas
    
    bias1 <- sum(mod-meas)/sum(meas)*100
    rmse1 <- (mean((mod-meas)^2))^0.5
    cor1 <- cor(meas,mod, method = "pearson")
  }else{
    bias1 <- NA
    rmse1 <- NA
    cor1 <- NA
  }
  
  stats <- data.frame(
    nObsHab =    nObsHab
    ,nObsNoHab =  nObsNoHab
    ,nPredHab =   length(pred[pred=="hab"])
    ,nPredNoHab = length(pred[pred %in% c("nohab","nodata")])
    ,pctAccHab =  ngh/nObsHab*100
    ,pctAccNoHab= ngnh/nObsNoHab*100
    ,pctAccAll =  (ngh + ngnh)/nobs*100
    ,chlPctBias = bias1
    ,chlRmse = rmse1
    ,chlCor = cor1
  )
  return(stats)
} # end getPstats function


getIDW <- function(pts, IDWdist ){
  print("start IDW")
  # loop through NA values, and assign chl value based on 
  # inverse dist weighting or neareast neighbor
  nn <- 1 # number of nearest points to use: set to 1 to use nearest neighbor, > 1 to use IDW
  idp <- 1 # power in IDW formula
  chl1 <- pts$chl
  
  mindist <- rep(1, nr)
  goodData <- pts[which(is.na(pts$chl)==FALSE),]
  X <- goodData$X
  Y <- goodData$Y
  chlData1 <- goodData$chl
  
  # determine whether there are goodData within Maumee Bay
  mbData <- FALSE
  mbX <- 13.7 # X and Y defining NE corner of Maumee Bay region
  mbY <- 46.56
  xx <- which(X < mbX & Y < mbY)
  if(length(xx) > 0){mbData <- TRUE}
  
  napts <- which(is.na(pts$chl) == TRUE)
  ro <- napts[1]
  for(ro in napts){
    x1 <- pts[ro,1]
    y1 <- pts[ro,2]
    dist <- ((x1-X)^2 + (y1-Y)^2)^0.5
    #  ind1 <- which(dist==min(dist))
    o <- rank(dist)
    ind1 <- which(o <= nn)
    dist1 <- dist[ind1]
    chl11 <- chlData1[ind1]
    weights <- dist1^(-idp) # IDW weights
    mindist[ro] <- mean(dist[ind1])
    # only fill in areas with good data closer than IDWdist
    if(min(dist[ind1]) <= IDWdist ){
      chl1[ro] <- sum(chl11*weights/sum(weights))
    }
    # if there are any goodData within Maumee Bay, fill in Maumee Bay nodes by IDW
    if(x1 < mbX & y1 < mbY & mbData){
      chl1[ro] <- sum(chl11*weights/sum(weights))
    }
  }
  
  pts$chl <- chl1
  
  print("finish IDW")
  print(paste(Sys.time()))
  return(pts)

} # end getIDW <- function(pts ){


setupPlot <- function(fn){
  graphics.off()
  width=4.5
  height=3.5
 # png(filename = fn, width = width, height = height, units = "in", pointsize = 10, res=300)
  # Joe wants 568x448 images for the website
  png(filename = fn, width = 568, height = 448, units = "px", pointsize = 14)
  
  par(mar=c(2, 0, 2, 0) + 0.1)#c(5, 4, 4, 2) + 0.1 c(bottom, left, top, right)
  mat <- matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE)
 # widths <- c(6,1.5)
  widths <- c(6,1.1)
  widths <- widths/sum(widths)
  heights <- c(1)
  heights <- heights/sum(heights)
  n <- layout(mat,widths=widths,heights=heights)
  #  layout.show(n)
}

logo1 <- readPNG("code/hab_tracker_2017_code_share/hab_tracker/noaa_glerl_ciler_logo.png")


# define a function to calculate percent HAB by zone polygons
getPctHab <- function(pDataHab){
  # pDataHab is a data frame with columns 
  # "chl", chlorophyll concentration
  # "PID", the FVCOM node

  # calculate the HAB area in each polygon
  pDataHab$hab <- 0
  pDataHab$hab[which(pDataHab$chl >= chlCrit)] <- 1
  df3 <- data.frame(
    art1 = art1
    ,poly = polys3_nodes
    ,hab = rep(NA,length(art1))
    ,chl = rep(NA,length(art1))
    ,stringsAsFactors = FALSE
  )
  df3$hab[pDataHab$PID] <- pDataHab$hab 
  df3$chl[pDataHab$PID] <- pDataHab$chl
  
  #   cols <- getColors(df3$chl)
  #   cols$cols[which(df3$hab==0)] <- "pink"
  #   cols$cols[pDataHab$PID] <- "pink"
  #   plotPolys(psTces, col=cols$cols, xlim=xlimLcc, ylim=ylimLcc)
  #   xx <- which(polys3_nodes=="81")
  #   points(coordsn$X[xx], coordsn$Y[xx])
  
  df3$hab[is.na(df3$hab)] <- 0
  df3$habArea <- df3$hab*df3$art1
  habArea <- aggregate(habArea ~ poly, data=df3, sum)
  habArea <- merge(habArea, polys3_areas)
  habArea$pctHab <- habArea$habArea/habArea$area*100
  return(habArea)
}

print("end habtracker_functions.R")
