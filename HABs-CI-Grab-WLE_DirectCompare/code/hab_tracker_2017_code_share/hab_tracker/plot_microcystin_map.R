# Mark Rowe 7-7-2017
# Make a map plot for HAB Tracker with Microcystin values at stations

source("habtracker_functions.R")

# Plot a map showing the 1D station locations

# plot surface temperature at the image time on the map plot
#values <- get.var.ncdf(nc2, "temp", start=c(1,1,xxi-1), count=c(-1,1,1))
# 
# stas1d2 <- projLcc(stas1d[stas2plot,], prj_new)
# xlim = c(xlimLcc[1], max(stas1d2$X)+20)
# ylim = c(0, 81.477247)
# #temp <- values

# read microcystin values from weekly data share file
datastoreDir <- "/mnt/datastore/OHH/WLE/"
#yyyy <- "2017"
fnwd <- Sys.glob(paste0(datastoreDir,yyyy,"/",yyyy,"_WLE_Weekly_Datashare*.csv"))
mtime <- file.mtime(fnwd)
fnwd <- fnwd[which(mtime==max(mtime, na.rm=TRUE))] # if there are multiple files, use the newest
print(fnwd)
if(file.exists(fnwd)){
  
  wdat <- read.csv(fnwd, header=TRUE, colClasses = "character" )
  #wdat <- wdat[,c(1,2,5,8,9,23,24)]
  names1 <- c("Date","Site","Sample_Depth_category","Lat_deg","Long_deg","Particulate_Microcystin_ugL.1","Dissolved_Microcystin_ugL.1")
  wdat <- wdat[,names1]
  names(wdat) <- c("datestring","station","depthclass","lat","lon","mcp","mcd")
  wdat <- wdat[is.na(wdat$datestring)==FALSE,]
  
  # find out if year has century
  nchar1 <- nchar(strsplit(wdat$datestring[nrow(wdat)],"/")[[1]][3])
  if(nchar1==4){
    wdat$date <- as.POSIXct(wdat$datestring, format="%m/%d/%Y")
  }else{
    wdat$date <- as.POSIXct(wdat$datestring, format="%m/%d/%y")
  }

  wdat$mcp[wdat$mcp %in% c("ns","nd")] <- NA
  wdat$mcp[wdat$mcp %in% c("bdl")] <- "0"
  wdat$mcd[wdat$mcd %in% c("ns","nd")] <- NA
  wdat$mcd[wdat$mcd %in% c("bdl")] <- "0"
  wdat$mcp <- as.numeric(wdat$mcp)
  wdat$mcd <- as.numeric(wdat$mcd)
  wdat$mct <- wdat$mcp + wdat$mcd
  brks <- c(-10,1e-6,0.3,1.6,20,1e32)
  wdat$mccut <- base::cut(wdat$mct, breaks = brks)
  wdat$lon <- as.numeric(wdat$lon)
  wdat$lat <- as.numeric(wdat$lat)
  
  date1 <- wdat$date[which(wdat$date==max(wdat$date,na.rm=TRUE))[1]] # as.POSIXct("2016-08-22")
  wdat1 <- wdat[wdat$date == date1,]
  wdat1 <- wdat1[is.na(wdat1$station)==FALSE,]
  #wdat1 <- merge(wdat1, stas1d[,c("station","lon","lat")])
  
  # use the actual lat lon in the weekly data share file
  # it's only filled in for the first depth class for each station so copy it down
  # for each depth class
  ustas <- unique(wdat1$station)
  for(stai in 1:length(ustas)){
    wdat2 <- wdat1[wdat1$station==ustas[stai],]
    lat <- wdat2$lat[is.na(wdat2$lat)==FALSE][1]
    if(is.na(lat)==TRUE){lat <- stas1d$lat[which(stas1d$station==ustas[stai])]}
    wdat1$lat[wdat1$station==ustas[stai]] <- lat
    lon <- wdat2$lon[is.na(wdat2$lon)==FALSE][1]
    if(is.na(lon)==TRUE){lon <- stas1d$lon[which(stas1d$station==ustas[stai])]}
    wdat1$lon[wdat1$station==ustas[stai]] <- lon
  }
  
  #wdat1 <- wdat1[wdat1$depthclass %in% c("Surface","Bottom"),]
  wdat1 <- wdat1[is.na(wdat1$date)==FALSE,]
  wdat1 <- projxy(wdat1, prj_new=prj_new)
  wdat1$X <- wdat1$X/1000
  wdat1$Y <- wdat1$Y/1000
  
  fn <- paste(plotPath, "microcystin_stationMap.png",sep="")
  
  graphics.off()
  width=4.5
  height=3.5
  # png(filename = fn, width = width, height = height, units = "in", pointsize = 10, res=300)
  # Joe wants 568x448 images for the website
  png(filename = fn, width = 568, height = 448, units = "px", pointsize = 14)
  
  par(mar=c(2, 0, 1, 0) + 0.1)#c(5, 4, 4, 2) + 0.1 c(bottom, left, top, right)
  mat <- matrix(c(1,2), nrow=2, ncol=1, byrow = TRUE)
  # widths <- c(6,1.5)
  widths <- c(1)
  widths <- widths/sum(widths)
  heights <- c(2.5,1)
  heights <- heights/sum(heights)
  n <- layout(mat,widths=widths,heights=heights)
  #  layout.show(n)
  
  xlim <- c(-7, 90)
  ylim <- c(13, 65)
  
  plotPolys(shoreline, projection=TRUE
            # ,plt = c(0.11, 0.98, 0.12, 0.88)
            ,plt = c(0.01, 0.99, 0.0, 0.9)
            # ,col=tim.colors(7)[1] # fill the background with the low-HAB color in case there are gaps
            ,xaxt="n",yaxt="n",bty="n",ylab="",xlab=""
            ,xlim=xlim, ylim=ylim
  )
  #box("figure")
  datePlot <- format(wdat1$date[1], "%Y-%m-%d")
  mtext(paste(datePlot), side=3, line=0.5, cex=2.0)
  
  cex1 <- 1.2
  dy <- 2
  cols <- c("blue","green","yellow","orange","red")
  wdat2 <- wdat1[wdat1$depthclass=="Surface",]
  cols1 <- cols[as.numeric(wdat2$mccut)]
  points(wdat2$X, wdat2$Y+dy, col="black", pch=16, cex=cex1+0.5)
  points(wdat2$X, wdat2$Y+dy, col=cols1, pch=16, cex=cex1)
  wdat2 <- wdat1[wdat1$depthclass=="Bottom",]
  cols1 <- cols[as.numeric(wdat2$mccut)]
  points(wdat2$X, wdat2$Y-dy, col="black", pch=17, cex=cex1+0.5)
  points(wdat2$X, wdat2$Y-dy, col=cols1, pch=17, cex=cex1)
  wdat2 <- wdat1[wdat1$depthclass=="Scum",]
  cols1 <- cols[as.numeric(wdat2$mccut)]
  points(wdat2$X, wdat2$Y, col="black", pch=15, cex=cex1+0.5)
  points(wdat2$X, wdat2$Y, col=cols1, pch=15, cex=cex1)
  
  addPolys(mask2, col="tan")
  plotCities(cex.cities=0.8)
  scalebar(2.9, 21.2, lwd=4)
  
  # plot the 1D station locations and labels
  cex1 <- 1.0
  # #xxb <- which(stas1d2$X == max(stas1d2$X))
  # xxb <- which(stas1d2$station %in% c("WE6","45176"))
  # xxnb <- which(!(1:nrow(stas1d2) %in% xxb))
  # x <- stas1d2$X
  # y <- stas1d2$Y
  # x[which(stas1d2$station == "WE25")] <- x[which(stas1d2$station == "WE25")] +1 # shift plotting position of WE25 slightly
  # points(x, y, pch=3)
  # text(stas1d2$X[xxnb], stas1d2$Y[xxnb], stas1d2$station[xxnb], font=2, pos=4, cex=cex1)
  # text(stas1d2$X[xxb], stas1d2$Y[xxb], stas1d2$station[xxb], font=2, pos=2, cex=cex1)
  
  ustas <- unique(wdat1$station)
  for(stai in 1:length(ustas)){
    wdat2 <- wdat1[wdat1$station==ustas[stai][1],]
    x <- wdat2$X
    y <- wdat2$Y
    points(x, y, pch=3)
    pos=4
    if(wdat2$station[1] == c("WE9")){pos=2}
    # only label the standard stations to avoid overlapping labels
    if(wdat2$station[1] %in% c("WE2","WE6","WE8","WE9","WE4","WE12","WE13","WE15")){
      text(x,y, wdat2$station, font=2, pos=pos, cex=cex1)
    }
  }
  
  logo.draw( logo1, 0.45, 0.02, 11.0)
  
  par(mar=c(0, 0, 0, 0) + 0.1)#c(5, 4, 4, 2) + 0.1 c(bottom, left, top, right)
  plot(c(0,10),c(0,10),type="n", yaxt="n", xaxt="n", xlab=""
       , ylab="", bty="n")
  
  mlab <- expression(paste("Microcystins ", mu, g, " ", L^{-1}))
  text(-0.1,9.0,mlab,pos=4, cex=1.4)
  #plot(NA,NA)
  legend("bottomleft", bty="n", pt.cex = 1.4, ncol=2,
         legend=c("Scum","Surface","Bottom","","","not-detected","detect-0.29","0.3-1.59","1.6-19.9",">= 20")
         ,pch=c(22,1,2,NA,NA,rep(16,length(cols)))
         ,col=c("black","black","black",NA,NA,cols)
  )
  
  par(def.par)
  graphics.off()
  
}else{
  print(paste("file does not exist",fnwd))
}
