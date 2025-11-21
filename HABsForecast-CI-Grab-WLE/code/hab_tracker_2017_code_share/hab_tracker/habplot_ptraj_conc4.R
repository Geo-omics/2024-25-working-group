# Mark Rowe 12-16-2014
# R script to plot output from habtracker model
# using ptraj

#.libPaths("~/zebra2/users/rowe/rlib")
require(ncdf)
require(raster)
require(sp)
require(PBSmapping)
require(fields)
require(plotrix)
def.par <- par(no.readonly = TRUE) # save default, for resetting...

par(def.par)
graphics.off()

# run the next two lines to load up variables to run this script outside of habtracker.R
#  dryrun <- TRUE
#  source("habtracker.R")

workdir <- habtrackerDir # paste0(getwd(),"/") # '~/zebra2/users/rowe/hab/hab_tracker/parallel1/'
# fntrknc <- 'erienc-ptraj.nc'
# fntrkfc <- 'eriefc-ptraj.nc'

#workdir <- '.'
#plotPath <- paste(workdir,"p3d_plot/tmp/", sep="")
#plotPath <- paste("./p3d_plot/tmp/", sep="")

print(" load common data and functions:  habtracker_functions.R")
source("habtracker_functions.R")

#fvcomdirnc="~/zebra2/users/rowe/hab/fvcom_netcdf/"

# # these numbers should be the same as the ones that were used to initialize 
# # the number of particles in image2ini.pro script
# concdep=.1  		#;assume depth (meters) for image resolution (????)
# massperpart=1e+6	#;scale # particles to ug (mass) of algae

# # read parameter values set in hindcast_main.R
# p <- read.table("parameters.txt")
# massperpart <- p[which(p[,1]=="massperpart"),2]
# concdep <- p[which(p[,1]=="concdep"),2]
# concdep <- 0.1 # assume depth in m of grid cell (only used for doIDL)
# massperpart <- 1e+9 # 1e+8   # microgram of chl mass per particle

if(doSkill){
  imagePath <- paste0(zebraDir,"/users/rowe/hab/hab_text/")
  wqPath <- paste0(zebraDir,"/users/rowe/hab/wfiles/")
  pPath <- paste0(zebraDir,"/users/rowe/hab/persistence_composite/") 
  # path containing plots of the satellite data to place alongside the skill stat images
  #imagePlotPath <- "../Rimages/" 
  imagePlotPath <- paste0(zebraDir,"/users/rowe/hab/Rimages/")
  
  # get image file names
  # the zone polys, polys3, are loaded here
  load(paste0(imagePath,"ImageInventory.Rdata"))
  imageFiles <- paste(df1$fnplot)
  #imageFiles <- list.files(imagePlotPath, "*.png")
  #imageDates <- substr(imageFiles, start=13, stop=21)
  #imageFiles2 <- list.files("~/zebra2/users/rowe/hab/hab_text/", "*.txt")
  imageFiles2 <- paste(df1$fntxt)
  imageDates <- substr(imageFiles2, nchar(imageFiles2)-16, nchar(imageFiles2)-8)
  
  # get water quality file names
  wfiles <- Sys.glob(paste(wqPath, "W*.txt",sep=""))
  datestringsw <- substr(wfiles, nchar(wfiles)-16, nchar(wfiles)-10)
  
} # end if doSkill
# 
# # limits for map plots
# xlimLL <- c(-83.47336, -82.01775)
# ylimLL <- c(41.36373, 42.08565)

#Get file start time
# fn <- paste("timestamp.txt",sep="")
# print(fn)
# date1 <- read.table(fn )
year1 <- as.numeric(yyyy) # date1[,1]
day1 <- as.numeric(ddd) #date1[,2]
hour1 <- as.numeric(hh) #date1[,3]
#datestring0 <- paste(sprintf("%04.0f", year1),sprintf("%03.0f", day1),sprintf("%02.0f", hour1),"0000",sep="")
datestring0 <- paste(yyyy,ddd,hh,"0000",sep="")
print(paste("datestring0 = ",datestring0))
year <- year1
day <- day1
hour <- hour1-1

# this is done in habtracker.R
# # # delete old image files
#  string1=paste('rm ',plotPath,'*.png',sep="") # ;mdr
#  system(string1)

framecount=0
skillFiles <- c()
statscount=0

############################################
# loop though nowcast forecast to load up nc files from parallel dirs
ic = 1
for(ic in 1:2){
  #for(ic in 1:1){ #plot nowcast only
  
  if(ic == 1){
    #;=================================================================
    #  ;    				Nowcast particles
    #;=================================================================
    print('------Loading Nowcast------')
    
    # fntrk <- paste(workdir,'erienc-ptraj.nc',sep="")
    # fntrk <- paste(hindcastArchive3,'erienc-ptraj.nc',sep="")
    fntrk <- fntrknc
    
  }else{
    #;=================================================================
    #  ;      			Forecast particles
    #;=================================================================
    print('------Loading Forecast------')
    
    #  fntrk <- paste(workdir,'eriefc-ptraj.nc',sep="")
    fntrk <- fntrkfc
    
  }
  
  trkList <- list() # empty list to hold parallel nc files
  for(ri in 1:npdirs){
    pdir <- paste0(habtrackerDir, pdirs[ri],"/")
    fntrk1 <- paste0(pdir,fntrk)
    if(file.exists(fntrk1)){
      print(paste("reading",fntrk1))
      trkList[[ri]] <- open.ncdf(fntrk1, readunlim=FALSE)
    }else{
      print(paste("nc file does not exist, fntrk1=",fntrk1))
      break
    }
  } # end ri
  
  # save the initial particle positions to plot displacement arrows
  if( ic==1){
    trkListnc <- trkList
  }else{
    trkListfc <- trkList
  }
} # end ic

############################################
ic=1
# loop though nowcast forecast to make plots
for(ic in 1:2){
  # select nowcast or forecast output to plot
  if( ic==1){
    trkList <- trkListnc 
    print('------Plotting Nowcast------')
  }else{
    trkList <- trkListfc 
    print('------Plotting Forecast------')
  }
  
  it=1
  tind <- 0 # time index in output file
  
  elev = get.var.ncdf(trkList[[1]], "elev", start=c(1,1), count=c(1,-1))
  ncrecords <- length(elev)
  for(it in 1:ncrecords){
    
    #   # apply growth rate (assumes hourly records)
    #   mu <- 0.1 # growth rate 1/day
    #   dt <- 1.0/24.0 # time increment, days
    #   massperpart <- massperpart + massperpart*mu*dt
    
    hour <- hour + 1
    tind <- tind + 1
    if(hour >= 24){
      hour <- hour - 24
      day <- day +1
    }
    date <- as.POSIXct(paste(year,day,hour), format="%Y %j %H", tz="GMT")
    tzoffset <- -4*3600 # time zone offset from GMT for plot time stamp, seconds
    tzstring <- "EDT"   # name of the time zone for plot time stamp
    datePlot <- paste(format(date+tzoffset, "%Y-%m-%d %H:%M"), tzstring)
    datestring <- paste(sprintf("%04.0f", year),sprintf("%03.0f", day),sprintf("%02.0f", hour),"0000",sep="")
    datestring1 <- substr(datestring, start=1, stop=9)
    
    if(doSkill){
      # check if there is an image for this date from which to calculate skill stats
      # if not, skip to the next hour in "it" loop
      #  if(!(datestring1 %in% imageDates) & it+ic != 2){next} # comment this line out to plot each hour
    } # end if doSkill
    
    print(paste("make plots for",datestring1))
    
    # include sequential framecount in the filename for use in animation
    framecount <- framecount + 1
    fileChar <- as.character(framecount)
    char1 <- paste(substr("000",start=1,stop=3-nchar(fileChar)),fileChar,sep="")
    datestringa <- paste(substr(datestring,start=1,stop=nchar(datestring)-nchar(char1)),char1,sep="")
    fn <- paste(plotPath, datestringa, ".png",sep="")
    
    # read 3D concentration from CONCOUT ptraj output
    # conc is output in particles/ m3
    # loop over output files from parallel runs and average the conc output arrays
    nc <- trkList[[1]] 
    conc <- get.var.ncdf(nc, "conc", start=c(1,1,tind), count=c(-1,-1,1))*0.0 #get array size
    for(ri in 1:npdirs){
      nc <- trkList[[ri]]
      conc <- conc + get.var.ncdf(nc, "conc", start=c(1,1,tind), count=c(-1,-1,1))*massperpart/1000
    }
    conc <- conc/ri
    
    if(plot2D){ 
      # particles were initiated on surface only, so plot as if all particles are in the surface layer
      conc <- apply(conc, MARGIN=1, FUN=sum)
    }else{
      # plot the surface layer concentration 
      conc <- conc[,1]
    }
    
    if(doSkill){
      # check if there is an image for this date, and if so, calculate skill stats
      #   if(datestring1 %in% imageDates){
      #     fi <- which(imageDates == datestring1)[1]
      fis <- which(imageDates == datestring1)
      fi <- fis[1]
      if(length(fis) > 0){
        #  break
        if(exists("stats2")){rm(stats2)}
        
        for(fi in fis){
          
          imageFile <- imageFiles[fi]
          if(substr(imageFile, 12,12) == "W"){
            imageFile2 <- wfiles[which(datestringsw == substr(datestring1, 1,7)) ]
            wqonly <- TRUE
            idata1 <- 2
          }else{
            wqonly <- FALSE
            idata1 <- 1
            imageFile2 <- imageFiles2[fi]
          }
          if(file.exists(imageFile2)==FALSE){imageFile2 <- gsub("~/zebra2/", zebraDir, imageFile2) }
          
          path1 <- paste(imagePlotPath,imageFile,".png",sep="")
          
          # get the initial image for the persistence forecast
          fip <- which(imageDates == substr(datestring0,1,9))[1]
          imageFilep <- imageFiles2[fip]
          nchar1 <- nchar(imageFilep)
          imageFilep <- paste0(pPath, substr(imageFilep, nchar1-17, nchar1))
          
          print("Matched model date with image date")
          print(paste("Model date",datestring1))
          print(paste("Image file",imageFile))
          print(paste("txt file",imageFile2))
          print(paste("persistence file",imageFilep))
          
          if(exists("pData")==FALSE){
            print(paste("reading persistence data, imageFile = ",imageFilep))
            pathp <- imageFilep
            pData <- read.table(pathp, header=FALSE)
            names(pData) <- c("lon","lat","chl","temp")
            pData$PID <- seq(1,nrow(pData))
            pData <- pData[is.na(pData$chl)==FALSE,]
            names(pData)[which(names(pData)=="chl")] <- "pred"
          }
          
          skillFiles <- c(skillFiles, fn, path1)
          
          # check if there are water quality data for this date
          fiw <- which(datestringsw == substr(datestring1, 1,7))[1]
          wfile <- wfiles[fiw]
          if(is.na(wfile)==FALSE & wqonly == FALSE){
            ndata <- 2 
          }else{
            ndata <- 1  
          }
          
          # calculate skill stats first for satellite data, then for water quality data
          idata=1
          for(idata in idata1:ndata){
            
            if(idata == 1){ # satellite data: current image vs forecast
              print(paste("reading satellite data, imageFile = ",imageFile2))
              # path2 <- paste(imagePath,substr(imageFile, 12, nchar(imageFile)-4),".txt", sep="")
              path2 <- imageFile2
              chlData <- read.table(path2, header=FALSE)
              names(chlData) <- c("lon","lat","chl","temp")

              # apply quality criteria and subset values within model domain 
              chlData <- processChlData(chlData)
              
              # screen data using the buffer function defined in image2ini_ptraj_ug.R
              chlData <- bufferFun(chlData, buffer_rad, buffer_rad2, fill=TRUE)
         
              if(nrow(chlData) < 1){break}
              
              #                             plotPolys(shorelineLL, projection=TRUE, bg="wheat")
              #                             cols <- getColors(chlData$chl)
              #                             points(chlData$lon, chlData$lat, col=cols$cols, pch=16, cex=0.5)
              
              
              #                          # graphics.off()
              #                           plotPolys(shoreline, projection=TRUE, bg="wheat", xlim=xlimLcc, ylim=ylimLcc)
              #                           cols <- getColors(chlData$chl)
              #                           points(chlData$X, chlData$Y, col=cols$cols, pch=15, cex=0.8)
              #                           plotPolys(shoreline, projection=TRUE, bg="wheat", xlim=xlimLcc, ylim=ylimLcc)
              #                           cols <- getColors(pData$pred)
              #                           points(pData$X, pData$Y, col=cols$cols, pch=15, cex=0.8)
              #                          plotPolys(shoreline, projection=TRUE, bg="wheat", xlim=xlimLcc, ylim=ylimLcc)
              #                          cols <- getColors(pData$chl)
              #                          points(pData$X, pData$Y, col=cols$cols, pch=15, cex=0.8)
              
              # aggregate  chl data so there is one mean value per FVCOM node
              chlData <- aggregate(chl ~ PID, data=chlData, mean)
              
              # associate persistence chl data with current image chl data

              pDatas <- merge(pData, chlData[,c("PID","chl")])
              
              statsp <- getPstats(pDatas, chlCrit=chlCrit) # persistence forecast stats
              statsp$type <- "Persistence"
              
            }else{ # water quality data, idata=2
             # break
              print(paste("reading water quality data, wfile = ",wfile))
              chlData <- read.table(wfile, header=FALSE)
              names(chlData) <- c("lon","lat","chl","temp")
              chlData <- projLcc(chlData, prj_new)
              chlData$EID <- 1:nrow(chlData)
              fp <- findPolys(events=chlData, polys=psTces)
              chlData <- merge(chlData, fp[,c("EID","PID")])
              chlData <- aggregate(chl ~ PID, data=chlData, mean)
              # PID is the FVCOM node associated with the satellite obs
              pDatao <- merge(pData, chlData[,c("PID","chl")])
              statsp <- getPstats(pDatao, chlCrit=chlCrit) # persistence forecast stats
              statsp$type <- "PersistenceObs"
            }
            
            if(nrow(chlData) > 0){
              print(paste("calculate skill stats"))
              
              #  chlPred <- extract(rpd, chlData[,c("X","Y")])
              chlPred <- conc[chlData$PID]
              chlData$pred <- chlPred
              
              stats <- getPstats(chlData, chlCrit=chlCrit)
              if(idata == 1){
                stats$type <- "Satellite"
              }
              if(idata == 2){
                stats$type <- "Observed"
              }
              
              if(exists("stats2")){
                stats2 <- rbind(stats2, statsp, stats)
              }else{
                stats2 <- rbind(statsp,stats)          
              }
              
              if(statscount==0){
               # statscount <- 1
                statsOut0 <- data.frame(
                  year = year
                  ,day = day
                  ,hour = hour
                  ,datestring0 = datestring0
                  ,datestring = datestring
                  ,fcday = day - day1
                  #  ,idata = idata
                )
                statsOut <- cbind(statsOut0, stats)
                statsOutp <- cbind(statsOut0, statsp)
                statsOut <- rbind(statsOut, statsOutp)
                chlDataOut1 <- chlData
                chlDataOut1$fcday <- day - day1
                chlDataOut <- chlDataOut1
                
                if(idata==1){
                  pDataOut1 <- pDatas                  
                }else{
                  pDataOut1 <- pDatao                  
                }
                pDataOut1$fcday <- day - day1
                pDataOut <- pDataOut1
                
              }else{
                statsOut0 <- data.frame(
                  year = year
                  ,day = day
                  ,hour = hour
                  ,datestring0 = datestring0
                  ,datestring = datestring
                  ,fcday = day - day1
                  #  ,idata = idata
                )
                statsOut1 <- cbind(statsOut0, stats)
                statsOutp <- cbind(statsOut0, statsp)
                statsOut <- rbind(statsOut, statsOut1)
                statsOut <- rbind(statsOut, statsOutp)
                
                chlDataOut1 <- chlData
                chlDataOut1$fcday <- day - day1
                chlDataOut <- rbind(chlDataOut, chlDataOut1)
                
                if(idata==1){
                  pDataOut1 <- pDatas                  
                }else{
                  pDataOut1 <- pDatao                  
                }
                pDataOut1$fcday <- day - day1
                pDataOut <- rbind(pDataOut, pDataOut1)
                
              } # end if(statscount==0){
              
              # for pct HAB by polygon
              if(idata ==1){
              if(statscount==0){
                statscount <- 1
                pDataHab <- pData
                pDataHab$chl <- pDataHab$pred
                habAreap <- getPctHab(pDataHab)
                concData <- data.frame(PID=1:length(conc), chl=conc)
                habAream <- getPctHab(concData)
                
                # only use pct HAB data when NA area is less than 20%
                xx <- which(pctNaOut > 20, arr.ind=TRUE)
                pctHabOut[xx] <- NA
                
                df4 <- data.frame(t(habAreap$pctHab))
                names(df4) <- habAreap$poly
                df4$type <- "Persistence"
                df4 <- cbind(statsOut0, df4)
                
                df5 <- data.frame(t(habAream$pctHab))
                names(df5) <- habAreap$poly
                df5$type <- "Model"
                df5 <- cbind(statsOut0, df5)
                df4 <- rbind(df4, df5)
                
                # find the row in pctHabOut corresponding to the date
                # then get the observed pctHab and pctNa
                ro1 <- which(substr(datestrings,1,7)==substr(datestring,1,7))[1]
                df5 <- data.frame(t(pctHabOut[ro1,]))
                names(df5) <- habAreap$poly
                df5$type <- "Satellite"
                df5 <- cbind(statsOut0, df5)
                df4 <- rbind(df4, df5)
                pctHabSkill <- df4

              }else{
                
                concData <- data.frame(PID=1:length(conc), chl=conc)
                habAream <- getPctHab(concData)
                
                df5 <- data.frame(t(habAream$pctHab))
                names(df5) <- habAreap$poly
                df5$type <- "Model"
                df5 <- cbind(statsOut0, df5)
                df4 <- df5
                
                # find the row in pctHabOut corresponding to the date
                # then get the observed pctHab and pctNa
                ro1 <- which(substr(datestrings,1,7)==substr(datestring,1,7))[1]
                df5 <- data.frame(t(pctHabOut[ro1,]))
                names(df5) <- habAreap$poly
                df5$type <- "Satellite"
                df5 <- cbind(statsOut0, df5)
                df4 <- rbind(df4, df5)
                pctHabSkill <- rbind(pctHabSkill, df4)
                
              }} # end pct HAB polygon
              
              
            } # end if nrow chlData > 0
          } # end idata
        } # end loop over fis
      } # end  if(datestring1 %in% imageDates){
    } # end if doSkill
    
    # plot the image
    print(paste("plot",datestringa))
    
    
    values <- conc
    values[values < 0.1] = 0.1 # don't plot zero values as no-data 
    
    zoom=TRUE
    if(zoom ){ #zoom into an area by clicking on the plot
      # you have to make the plot first, and don't run par(def.par) after making the plot
      # xx <- locator(2) # click upper left, and lower right of zoom region
      xlim = xlimLcc # c(min(xx$x),max(xx$x))
      ylim = ylimLcc # c(min(xx$y),max(xx$y))
      psClip <- clipPolys(psTces, xlim=xlim, ylim=ylim)
      temp <- values[unique(psClip$PID)]   
    }else{
      xlim=NULL
      ylim=NULL
      temp <- values
      psClip <- psTces
    }
    
    cols <- getColors(temp)
    
    brk <- cols$brk
    brkm <- max(brk)
    label <- zlab
    brkLabs <- round(brk, 0)
    
    nlevel <- length(brk)-1
    zlim <- cols$zlim
    colorbar <- cols$colorbar
    
    require(akima)
    plotres <- 0.5 # 0.2 resolution of interpolation grid, km
    xlength <- (xlim[2]-xlim[1])/plotres
    ylength <- (ylim[2]-ylim[1])/plotres
    xx <- which(coordsn$X <= xlim[2] & coordsn$X >= xlim[1] & coordsn$Y <= ylim[2] & coordsn$Y >= ylim[1])
    x=coordsn$X[xx]
    y=coordsn$Y[xx]
    z=values[xx]
    data.interp <- interp(x,y,z, duplicate="mean"
                          ,xo=seq(min(x), max(x), length = xlength)
                          ,yo=seq(min(y), max(y), length = ylength)
                          , extrap=FALSE
                          , linear=TRUE
    )
    
    setupPlot(fn)
    
    plotPolys(shoreline, projection=TRUE
              # ,plt = c(0.11, 0.98, 0.12, 0.88)
              ,plt = c(0.01, 1, 0.07, 1)
              ,col=colorbar[1] # fill the background with the low-HAB color in case there are gaps
              ,xaxt="n",yaxt="n",bty="n",ylab="",xlab=""
              ,xlim=xlim, ylim=ylim
    )
    mtext(datePlot, side=3, line=1, cex=2.0)
    
    #box("plot")
    #box("figure")
    #graphics.off()
    
    brk1 <- brk
    brk1[length(brk1)] <- 1e6
    image(data.interp, add=TRUE,  col=cols$colorbar, breaks=brk1)
    
    #box()
    #    addPolys(psTces, col=cols$cols, border=cols$cols, bg="gray")
    #    addPolys(psClip, col=cols$cols, border=cols$cols, bg="gray",lwd=0.8) #plot TCE polygons
    addPolys(mask2, col="tan")
    #   addPolys(mask2, col="transparent")
    #points(chlData$lon, chlData$lat, pch=16, cex=0.3, col=cols)
    #addPolys(shoreline, density=0)
    plotCities(cex.cities=0.6)
    scalebar(188, 15, lwd=4)
    
    if(!doSkill){ # plot the logo if not doing hindcast skill assessment
      logo1 <- readPNG("noaa_glerl_ciler_logo.png")
      #  logo.draw( logo1, 0.01, 0.01, 5.0)
      logo.draw( logo1, 0.02, 0.02, 11.0)
    }
    
    if(doSkill ){
      if( datestring1 %in% imageDates){
        # add stats table to bottom of plot
        if(exists("stats2")){
          statsf <- stats2
          #       # if there is more than one row, only plot the row with most obs
          #       xx2 <- which(statsf$nObsNoHab == max(statsf$nObsNoHab, na.rm=TRUE))[1]
          #       if(length(xx2) < 1){xx2 <- 1}
          #       statsf <- statsf[xx2,]
          
          statsf$chlCor <- sprintf("%.2f", statsf$chlCor)
          statsf[,1:9] <- apply(statsf[,1:9], MARGIN=2, FUN=function(x){sprintf("%.0f", x)})
          
          #   addtable2plot(x=5000,y=-20000,statsf[1,],cex=0.55)
          # addtable2plot(x=xlim[1]+(xlim[2]-xlim[1])/100,y=ylim[1]-(ylim[2]-ylim[1])/14,statsf,cex=0.45)
          cex.table=0.55
          addtable2plot(x=xlim[1]+(xlim[2]-xlim[1])/100,y=ylim[1]-(ylim[2]-ylim[1])/10,statsf,cex=cex.table)
          # rm(stats2)
        } # end if exists stats2
      } # end if( datestring1 %in% imageDates)
    } # end if(doSkill )
    
    # par(mar = c(10, 1, 10, 7) + 0.0)  #c(bottom, left, top, right)
    #  plot(1,1)
    plot(c(0,10),c(0,10),type="n", yaxt="n", xaxt="n", xlab=""
         , ylab="", bty="n")
    #box()
    # text(1,5,label,srt=90, cex=2)
    text(1,5,label,srt=90, cex=1.5)
    
    brk[1] <- 0
    lab.breaks=sprintf("%3s" ,brkLabs) #format(brkLabs)
    if(!doSkill){lab.breaks[length(lab.breaks)] <- "   "}
    
    image.plot(zlim=zlim
               ,col=colorbar 
               ,breaks= seq(zlim[1],zlim[2])
               ,lab.breaks=lab.breaks
               , legend.width = 6
               # ,legend.args=list(cex=2.0)
               #  ,legend.mar=18 #15.1
               ,legend.mar=20 #15.1
               , legend.only=T, add=F)
    
    par(def.par)
    graphics.off()
    
  } # end it loop over hours
  
} # end ic loop over nowcast, forecast

################################################################
if(doAnimation){
  # convert still images to an animation
  require(animation)
  ## Use animation package to gernerate gif file
  # http://www.r-bloggers.com/animated-images-of-arctic-sea-ice-extent-decline/
  ani.options(convert="/opt/local/ImageMagick-6.7.9/bin/convert")
  ani.options(interval= 0.20)
  
  # convert still images to animation  
  stillFiles <- Sys.glob(paste0(plotPath, "*.png"))
  # files to exclude from the animation
  excludeFiles <- c()
  excludeFiles <- c(excludeFiles, Sys.glob(paste0(plotPath, "depth-time*.png")))
  excludeFiles <- c(excludeFiles, Sys.glob(paste0(plotPath, "*fskill*.png")))
  excludeFiles <- c(excludeFiles, Sys.glob(paste0(plotPath, "*displacement*.png")))
  excludeFiles <- c(excludeFiles, Sys.glob(paste0(plotPath, "*conc*.png")))
  excludeFiles <- c(excludeFiles, Sys.glob(paste0(plotPath, "*initial*.png")))
  excludeFiles <- c(excludeFiles, Sys.glob(paste0(plotPath, "1D*.png")))
  stillFiles <- stillFiles[!(stillFiles %in% excludeFiles)]
  #path = "/Users/rowe/Documents/mi_gem_archive/plots/050920143/xcShorts/"
  ani.options(outdir = substr(plotPath,1,nchar(plotPath)-1))    # direct gif to directory
  fileID <- substr(datestring0, 1, 9)
  anifn <- paste( fileID, "_animation.gif", sep="")
  system("rm *animation.gif")
  im.convert(stillFiles, anifn)
  string1 <- paste("mv",anifn,plotPath)
  system(string1)
} # end if doAnimation
################################################################
print(" make a plot of particle displacement over the run")

# only a sample of particle displacements are plotted, so only
# open the first parallel run file

fntrk <- fntrknc
pdir <- paste0(habtrackerDir, pdirs[1],"/")
fntrk1 <- paste0(pdir,fntrk)
if(file.exists(fntrk1)){
  print(paste("reading",fntrk1))
  nc0 <- open.ncdf(fntrk1, readunlim=FALSE)
}else{
  print(paste("nc file does not exist, fntrk1=",fntrk1))
  break
}

fntrk <- fntrkfc
fntrk1 <- paste0(pdir,fntrk)
if(file.exists(fntrk1)){
  print(paste("reading",fntrk1))
  nc1 <- open.ncdf(fntrk1, readunlim=FALSE)
}else{
  print(paste("nc file does not exist, fntrk1=",fntrk1))
  break
}
nrecfc <- length(get.var.ncdf(nc1, "time", start=c(1), count=c(-1)))

trksi <- data.frame(
  lon = get.var.ncdf(nc0, "x", start=c(1,1), count=c(-1,1)) - 360.0
  ,lat = get.var.ncdf(nc0, "y", start=c(1,1), count=c(-1,1))
  ,z1 = get.var.ncdf(nc0, "z", start=c(1,1), count=c(-1,1))
  ,elev = get.var.ncdf(nc0, "elev", start=c(1,1), count=c(-1,1))
)
# add projected coordinates to data frame
trksi <- projLcc(trksi, prj_new)

skip <- floor(nrow(trksi)/50) # number of particles to skip when plotting arrows
skip <- max(skip,1)
xxArrows <- seq(1,nrow(trksi),skip)

# loop over output files from parallel runs and average the conc output arrays
# nc <- trkListfc[[1]] # use the final conc value for the plot
# conc <- get.var.ncdf(nc, "conc", start=c(1,1,tind), count=c(-1,-1,1))*0.0 #get array size
# for(ri in 1:npdirs){
#   nc <- trkListfc[[ri]]
#   conc <- conc + get.var.ncdf(nc, "conc", start=c(1,1,tind), count=c(-1,-1,1))*massperpart/1000
# }
# conc <- conc/ri

nc <- trkListnc[[1]] # use the INITIAL conc value for the plot
conc <- get.var.ncdf(nc, "conc", start=c(1,1,1), count=c(-1,-1,1))*0.0 #get array size
for(ri in 1:npdirs){
  nc <- trkListnc[[ri]]
  conc <- conc + get.var.ncdf(nc, "conc", start=c(1,1,1), count=c(-1,-1,1))*massperpart/1000
}
conc <- conc/ri

if(plot2D){ 
  # particles were initiated on surface only, so plot as if all particles are in the surface layer
  conc <- apply(conc, MARGIN=1, FUN=sum)
}else{
  # plot the surface layer concentration 
  conc <- conc[,1]
}

values <- conc
values[values < 0.1] = 0.1 # don't plot zero values as no-data

zoom=TRUE
if(zoom ){ #zoom into an area by clicking on the plot
  # you have to make the plot first, and don't run par(def.par) after making the plot
  # xx <- locator(2) # click upper left, and lower right of zoom region
  xlim = xlimLcc # c(min(xx$x),max(xx$x))
  ylim = ylimLcc # c(min(xx$y),max(xx$y))
  psClip <- clipPolys(psTces, xlim=xlim, ylim=ylim)
  temp <- values[unique(psClip$PID)]   
}else{
  xlim=NULL
  ylim=NULL
  temp <- values
  psClip <- psTces
}

cols <- getColors(temp)

brk <- cols$brk
brkm <- max(brk)
label <- zlab
brkLabs <- round(brk, 0)

nlevel <- length(brk)-1
zlim <- cols$zlim
colorbar <- cols$colorbar

require(akima)
#plotres <- 200 # resolution of interpolation grid, m
xlength <- (xlim[2]-xlim[1])/plotres
ylength <- (ylim[2]-ylim[1])/plotres
xx <- which(coordsn$X <= xlim[2] & coordsn$X >= xlim[1] & coordsn$Y <= ylim[2] & coordsn$Y >= ylim[1])
x=coordsn$X[xx]
y=coordsn$Y[xx]
z=values[xx]
data.interp <- interp(x,y,z, duplicate="mean"
                      ,xo=seq(min(x), max(x), length = xlength)
                      ,yo=seq(min(y), max(y), length = ylength)
                      , extrap=FALSE
                      , linear=TRUE
)

date <- as.POSIXct(datestring0, format="%Y%j%H%M%S", tz="GMT")
fcday <- 0
for(fcday in seq(0,min(5,nrecfc/24),1)){
  

  
  fn <- paste(plotPath, datestring0, "displacement",fcday,".png",sep="")
  
  setupPlot(fn)
  
  plotPolys(shoreline, projection=TRUE
            # ,plt = c(0.11, 0.98, 0.12, 0.88)
            ,plt = c(0.01, 1, 0.07, 1)
            ,col=colorbar[1] # fill the background with the low-HAB color in case there are gaps
            ,xaxt="n",yaxt="n",bty="n",ylab="",xlab=""
            ,xlim=xlim, ylim=ylim
  )

  datePlot <- paste(format(date+tzoffset+fcday*24*3600+(24-hour1)*3600, "%Y-%m-%d %H:%M"), tzstring)
  mtext(datePlot, side=3, line=1, cex=2.0)
  
  
  brk1 <- brk
  brk1[length(brk1)] <- 1e6
  image(data.interp, add=TRUE,  col=cols$colorbar, breaks=brk1)
  
  addPolys(mask2, col="tan")
  plotCities(cex.cities=0.6)
  scalebar(188, 15, lwd=4)
  
  logo1 <- readPNG("noaa_glerl_ciler_logo.png")
  #  logo.draw( logo1, 0.01, 0.01, 5.0)
  logo.draw( logo1, 0.02, 0.02, 11.0)
  
  # get the final particle positions
  tind <- fcday*24 # time record to read
  if(tind==0){tind=1}
  
  # nc1 <- open.ncdf("/mnt/projects/hpc/hab/hab_tracker_2017/hab_tracker/parallel5/eriefc-ptraj.nc", readunlim=FALSE)
  # plot(get.var.ncdf(nc1, "y", start=c(1,1), count=c(1,-1)))
  # concx <- get.var.ncdf(nc1, "conc", start=c(1,1,65), count=c(-1,-1,1))
  # range(concx)
  
  trksf <- data.frame(
    lon = get.var.ncdf(nc1, "x", start=c(1,tind), count=c(-1,1)) - 360.0
    ,lat = get.var.ncdf(nc1, "y", start=c(1,tind), count=c(-1,1))
    ,z1 = get.var.ncdf(nc1, "z", start=c(1,tind), count=c(-1,1))
    ,elev = get.var.ncdf(nc1, "elev", start=c(1,tind), count=c(-1,1))
  )
  # add projected coordinates to data frame
  trksf <- projLcc(trksf, prj_new)
  # summary(trksf)
  # 
  # # add displacement arrows to the plot
  
  arrows(x0=trksi$X[xxArrows],y0=trksi$Y[xxArrows],x1=trksf$X[xxArrows],y1=trksf$Y[xxArrows]
         ,col="white", length=0.06, lwd=2)
  arrows(x0=trksi$X[xxArrows],y0=trksi$Y[xxArrows],x1=trksf$X[xxArrows],y1=trksf$Y[xxArrows]
         ,col="black", length=0.06, lwd=1)
         
         
  # 
  # # plot selected particle trajectories
  # parti <- xx[1]
  # for(parti in xx){
  # lon0 = get.var.ncdf(nc0, "x", start=c(parti,1), count=c(1,-1)) - 360.0
  # lat0 = get.var.ncdf(nc0, "y", start=c(parti,1), count=c(1,-1))
  # lon1 = get.var.ncdf(nc1, "x", start=c(parti,1), count=c(1,-1)) - 360.0
  # lat1 = get.var.ncdf(nc1, "y", start=c(parti,1), count=c(1,-1))
  # track1 <- data.frame(
  #   lon = c(lon0,lon1)
  #   ,lat= c(lat0,lat1)
  #   )
  # track1 <- projLcc(track1, prj_new)
  # nrec <- nrow(track1)
  # points(track1$X, track1$Y, type="l"
  #        ,col="white", lwd=1)
  # points(track1$X, track1$Y, type="l"
  #        ,col="black",  lwd=0.5)
  # arrows(x0=track1$X[nrec-1],y0=track1$Y[nrec-1],x1=track1$X[nrec],y1=track1$Y[nrec]
  #        ,col="white", length=0.04, lwd=1)
  # arrows(x0=track1$X[nrec-1],y0=track1$Y[nrec-1],x1=track1$X[nrec],y1=track1$Y[nrec]
  #        ,col="black", length=0.04, lwd=0.5)
  # }
  
  # par(mar = c(10, 1, 10, 7) + 0.0)  #c(bottom, left, top, right)
  #  plot(1,1)
  plot(c(0,10),c(0,10),type="n", yaxt="n", xaxt="n", xlab=""
       , ylab="", bty="n")
  #box()
  # text(1,5,label,srt=90, cex=2)
  text(1,5,label,srt=90, cex=1.5)
  
  brk[1] <- 0
  lab.breaks=sprintf("%3s" ,brkLabs) #format(brkLabs)
  if(!doSkill){lab.breaks[length(lab.breaks)] <- "   "}
  image.plot(zlim=zlim
             ,col=colorbar 
             ,breaks= seq(zlim[1],zlim[2])
             ,lab.breaks=lab.breaks
             , legend.width = 6
             # ,legend.args=list(cex=2.0)
             #  ,legend.mar=18 #15.1
             ,legend.mar=20 #15.1
             , legend.only=T, add=F)
  
  par(def.par)
  graphics.off()
  
} # end fcday loop


##############################################################

#rm(trks, trks1, trksi, trksf) # delete the large arrays

if(doSkill){
  # write statsout file
  fnstatsout <- paste(plotPath, "statsout.txt",sep="")
  write.table(statsOut, file = fnstatsout, append = FALSE, row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  # write the chlDataOut, paired sim and obs
  fnstatsout <- paste(plotPath, "chlDataOut.Rdata",sep="")
  save(chlDataOut, file=fnstatsout)
  fnstatsout <- paste(plotPath, "pDataOut.Rdata",sep="")
  save(pDataOut, file=fnstatsout)
  fnstatsout <- paste(plotPath, "pctHabSkill.Rdata",sep="")
  save(pctHabSkill, file=fnstatsout)
  
  # combine pngs 
  
  print( "loop  through skillFiles and combine pngs")
  #skillFiles <- gsub("txt","png",skillFiles)
  skillFiles <- matrix(skillFiles, ncol=2, byrow=TRUE)
  #print(skillFiles)
  
  combinePngs <- function(file1, file2, sidebyside=TRUE){
    require(png)
    require(abind)
    if(sidebyside){
      along=2
    }else{
      along=1
    }
    
    #   img1 <- readPNG(file1)
    #   img2 <- readPNG(file2)
    # on mac png(bg="white") writes a 4-layer png file, 
    # but on wings X11 device it will write 3-layer png without transparency layer
    # so only read three layers so this will work on both devices; won't support transparent background  
    img1 <- readPNG(file1)[,,1:3]
    img2 <- readPNG(file2)[,,1:3]
    img3 <- abind(img1, img2, along=along)
    
    return(img3)
  }
  
  ro <- 1
  for(ro in 1:nrow(skillFiles)){
    file1 <- skillFiles[ro,1]
    file2 <- skillFiles[ro,2]
    print(paste(file1,file2))
    img3 <- combinePngs(file1, file2, sidebyside=TRUE)
    
    if(ro == 1){
      img4 <- img3
    }else{
      img4 <- abind(img4, img3, along=1)
    }
    
  } # end ro loop
  
  
  fnimg1 <- paste(plotPath, "fskill_",datestring0, ".png",sep="")
  writePNG(img4, fnimg1)
} # end if doSkill

#############################################################################
print(" Make depth-time colorbar plots at stations")

if(doSkill){
  # load data frame of stations that were visited more than once 2008-2014
  #save(stasdf, file="/Users/rowe~/zebra2/users/rowe/hab/hab_tracker/ErieWQstas.Rdata")
  load("ErieWQstas.Rdata")
  load("ErieSummary2008-2014.Rdata")
  wq <- ErieSummary[[2]]
  wq <- wq[wq$year == year1,c("date","year","jday","station","lon","lat","temp","chl")]
  wqtime <- as.POSIXct(format(wq$date, "%Y-%m-%d %H:00")) - 4*3600
}else{
  # make plots at the 1d stations
  stasdf <- stas1d[is.na(stas1d$station)==FALSE,]
} # end if doSkill

# fntrk <- paste(workdir,'erienc-ptraj.nc',sep="")
# nc0 <- open.ncdf(fntrk, readunlim=FALSE)
# fntrk <- paste(workdir,'eriefc-ptraj.nc',sep="")
# nc <- open.ncdf(fntrk, readunlim=FALSE)

nc0 <- trkListnc[[1]]
nc <- trkList[[1]]

zvals <- get.var.ncdf(nc0, "zvals", start=1, count=-1)
zres <- zvals[1]-zvals[2]

date0 <- as.POSIXct(paste(year1,day1,hour1), format="%Y %j %H", tz="GMT")

# open hydrodynamic files
# if this script is run by habtracker.R fvcomnc is defined there, if run alone define it here
if(exists("fvcomnc")==FALSE){
  fvcomdirnc <- paste0(zebraDir,"/users/rowe/hab/fvcom_netcdf/")
  fncur0 <- paste(fvcomdirnc,substr(datestring0,start=1,stop=7),"00.nc",sep="")
  fvcomnc <- fncur0
}
fncur0 <- fvcomnc # paste(fvcomdirnc,fncur,sep="")
nccur0 <- open.ncdf(fncur0, readunlim=FALSE)

if(exists("fvcomfc")==FALSE){
  fncur1 <- paste(fvcomdirnc,"forecast.nc",sep="")
  fvcomfc <- fncur1
}
fncur1 <- fvcomfc  
nccur1 <- open.ncdf(fncur1, readunlim=FALSE)
siglay <- get.var.ncdf(nccur0, "siglay", start=c(1,1), count=c(1,-1))

stai <- 5
for(stai in 1:nrow(stasdf)){
  sta <- stasdf$station[stai]
  node <- stasdf$node[stai]
  if(is.na(node)==TRUE){
    diff <- (coordsn$lon - stas1d$lon[stai])^2 + (coordsn$lat - stas1d$lat[stai])^2
    node <- coordsn$node[which(diff==min(diff))[1]]
  }
  print(paste("station=",sta," node=",node))
  
  #     conc0 <- get.var.ncdf(nc0, "conc", start=c(node,1,1), count=c(1,-1,-1))*massperpart/1000
  #     conc1 <- get.var.ncdf(nc, "conc", start=c(node,1,1), count=c(1,-1,-1))*massperpart/1000
  
  # loop over output files from parallel runs and average the conc output arrays
  nc <- trkListnc[[1]] 
  conc <- get.var.ncdf(nc, "conc", start=c(node,1,1), count=c(1,-1,-1))*0.0 #get array size
  for(ri in 1:npdirs){
    nc <- trkListnc[[ri]]
    conc <- conc + get.var.ncdf(nc, "conc", start=c(node,1,1), count=c(1,-1,-1))*massperpart/1000
  }
  conc0 <- conc/ri
  
  nc <- trkListfc[[1]] 
  conc <- get.var.ncdf(nc, "conc", start=c(node,1,1), count=c(1,-1,-1))*0.0 #get array size
  for(ri in 1:npdirs){
    nc <- trkListfc[[ri]]
    conc <- conc + get.var.ncdf(nc, "conc", start=c(node,1,1), count=c(1,-1,-1))*massperpart/1000
  }
  conc1 <- conc/ri
  conc1 <- cbind(conc0, conc1)
  
  timi <- ncol(conc1)
  time1 <- seq(date0, date0 + timi*3600, by="hour")
  
  nrecnc <- ncol(conc0)
  xxnc <- seq((24-nrecnc+1),24)
  #    xxnc <- (hour1-1):(ncrecords+hour1-1)
  
  temp0 <- get.var.ncdf(nccur0, "temp", start=c(node,1,1), count=c(1,-1,-1))
  temp1 <- get.var.ncdf(nccur1, "temp", start=c(node,1,1), count=c(1,-1,-1))
  temp1 <- cbind(temp0[,xxnc], temp1)
  #    temp1 <- temp0[,xxnc]
  
  el0 <- get.var.ncdf(nccur0, "zeta", start=c(node,1), count=c(1,-1))
  el1 <- get.var.ncdf(nccur1, "zeta", start=c(node,1), count=c(1,-1))
  el1 <- c(el0[xxnc], el1)    
  #    el1 <- el0[xxnc]
  
  # plot the water depth as constant set at mean of the surface elevation
  el1 <- el1*0 + mean(el1)
  
  ds <- zvals - zres/2
  depth <- -coordsn$depth[node] - max(el1, na.rm=TRUE)
  xx <- which(ds > depth)
  ds <- ds[seq(1,(max(xx)+1))]
  nd <- length(ds)
  dmax <- min(ds)-zres/2 # depth at bottom edge of image
  
#  detm <- t(conc1)
#  mat2 <- detm[1:timi,nd:1]
mat2 <- conc1[1:nd, 1:timi]
  
#   zlim <- c(0, max(mat2, na.rm=TRUE))
#   if(zlim[2]==0){zlim[2] <- 1.0}
#   #    zlim <- c(min(mat2, na.rm=TRUE), max(mat2, na.rm=TRUE))
#   # zlim <- c(0,135)
#   cols <- tim.colors(256)
  
  # plot with the standard HAB Tracker color bar chlorophyll scale
  mat2 <- mat2 + 1.0e-6 # don't give NA for zero values
  cols <- getColors(mat2)
  brk <- cols$brk
  brkm <- max(brk)
  label <- zlab
  brkLabs <- round(brk, 0)
  nlevel <- length(brk)-1
  zlim <- cols$zlim
  colorbar <- cols$colorbar
  mat2 <- as.numeric(cut(mat2, breaks=brk))
  mat2 <- matrix(mat2, nrow=nd, ncol=timi, byrow=FALSE)
  
  rt <- raster(nrows=nrow(mat2), ncol=ncol(mat2), xmn=0, xmx=1.0, ymn=dmax, ymx=0)
  rt <- setValues(rt, mat2)
  
  fn <- paste(plotPath,"depth-time_",sta,".png",sep="")
#  png(file = fn, width=7, height=4, units="in", res=300, bg="white", pointsize=10)
#  par(mar = c(6, 4, 3, 5) + 0.1) #c(bottom, left, top, right)
  png(filename = fn, width = 568, height = 448, units = "px", pointsize = 14
  #, type="Xlib", antialias="subpixel"
  ) # type="Xlib" will prevent occasional stray vertical white line on image plot, bear uses type="cairo" by default 
  par(mar = c(8.5, 4, 2.5, 8.0) + 0.1) #c(bottom, left, top, right)
  
#    image(y=rev(ds), z=mat2, zlim=zlim, col=colorbar , xaxt="n", yaxt="n"
    image(rt, zlim=zlim, col=colorbar , xaxt="n", yaxt="n"
        ,xlab=" "
        ,ylab="Depth, m"
        ,ylim=c(depth[1], 0)
        ,main=paste("Station",sta)
        #  ,cex.main=0.8
  )
    
     hrs <- format(time1-3600*4, "%H")
     xx <- which(hrs %in% c("00","06","12","18")) # seq(1, nhr, 6)
     xx1 <- which(hrs %in% c("00")) # seq(1, nhr, 6)
     if(length(xx) <= 20){
       labels <- format(time1[xx]-3600*4, "%H")
     }else{
       labels <- FALSE
     }
     axis(1, at=(xx-1)/(length(hrs)-1), labels=labels, las=2)
     axis(1, at=(xx1-1)/(length(hrs)-1), labels=format(time1[xx1]-3600*4, "%Y-%m-%d %H"), las=2)
    
#   # skip <- timi/12 # number of times to skip for axis labels
#   # xx <- seq(1,timi,skip)
#   xx <- seq(as.POSIXct(format(time1[1], "%Y-%m-%d"),tz="GMT"),time1[length(time1)],by="day")
#   xx <- which(time1 %in% xx)
#   xx <- c(xx,timi)
#   ts <- time1[xx]
#   at=(xx-1)/(max(xx)-1)
#   labels=format(ts, "%m-%d")
#   xx <- 1:(length(labels)-1)
#   axis(1, at=at[xx], labels=labels[xx], las=2)
  
  dlabs <- seq(0,depth-2*zres,-zres)
  if(length(dlabs) <= 5){   dlabs <- seq(0,depth-2*zres,-zres*0.5)}
  if(length(dlabs) > 10){   dlabs <- seq(0,depth-2*zres,-zres*2)}
  at <- dlabs
  #  axis(2, at=0:(nd-1)/(nd-1), labels=format(round(rev(ds),0), nsmall=0), las=2)
  axis(2, at=at, labels=format(round(dlabs,1), nsmall=1), las=2)
#  mtext(year, side = 1, line = 4, cex=1.0)
#  mtext(sta, side = 3, line = 2, cex=1.0)
  # mtext(nutLab, side = 3, line = 0.5)
#  mtext(zlab, side = 3, line = 0.5)
  
#   # mask out plot area below the bottom
#   depth <- -coordsn$depth[node] - el1
#   depth <- c(depth, depth[length(depth)])
#   depth <- depth[1:timi]
#   timi1 <- timi+1
#   #  xx <- c(1:timi,timi:1)
#   xx <- c(0:timi1,timi1:0)
#   at=(xx-1)/(timi-1)
#   # y <- c(depth,(depth*0+dmax))
#   y <- c(depth[1],depth,depth[length(depth)],dmax,(depth*0+dmax),dmax)
#   polygon(at,y, col="black")
  
#   # plot temperature contours
#   depth <- (coordsn$depth[node] + mean(el1))*siglay
#   xx <- c(1:(timi-1))
#   at=(xx-1)/(timi-1)
#   detm <- t(temp1[,xx])
#   mat2 <- detm[1:nrow(detm),ncol(detm):1]
#   contour(x=at, y=rev(depth), z=mat2, add=TRUE, col="gray", lwd=1, vfont=c("sans serif", "bold"))
  
  if(doSkill){
    #     # plot obs 
    tind <- which(wqtime %in% time1)
    wq1 <- wq[tind,]
    wqtime1 <- wqtime[tind]
    pts <- matrix(c(-wq1$lon, wq1$lat), nrow=nrow(wq1), ncol=2)
    pt <- c(coordsn[node, "lon"], coordsn[node, "lat"])
    dist <- spDistsN1(pts, pt, longlat = TRUE)
    tind <- which(dist < 3)
    wqtime1 <- wqtime1[tind]
    wq1 <- wq1[tind,] # only plot water quality data within 3km of the node
    
    if(length(tind) > 0){
      cols2 <- fields::color.scale(wq1$chl, col=cols, zlim=zlim)
      int <- difftime(time1[timi],time1[1],unit="secs") 
      xvals <- as.numeric(difftime(wq1$date, time1[1], unit="secs"))/as.numeric(int)
      yvals <-  -0.5
      cex <- 0.8
      points(xvals, yvals, col=cols2, pch=19, cex=cex)
      points(xvals, yvals, col="white", cex=cex)
      graphics::text(xvals, -0.75, round(wq1$temp,1), col="gray", cex=cex )
    }
  } # end if doSkill
  
  #  logo.draw( logo1, 0.01, 0.01, 5.0)
  logo.draw( logo1, 0.67, 0.005, 6.0)

  mtext("Time, EDT", side=1, line=7)
  
#  image.plot(zlim=zlim, legend.only=TRUE)
#  plot(1,1)
  
  mtext(zlab, side=4, line=1.4)
  
  brk[1] <- 0
  lab.breaks=sprintf("%3s" ,brkLabs) #format(brkLabs)
  if(!doSkill){lab.breaks[length(lab.breaks)] <- "   "}
  image.plot(zlim=zlim
             ,col=colorbar 
             ,breaks= seq(zlim[1],zlim[2])
             ,lab.breaks=lab.breaks
             , legend.width = 1.5
             # ,legend.args=list(cex=2.0)
             #  ,legend.mar=18 #15.1
             ,legend.mar=5 #15.1
             , legend.only=T
             #, add=F
             )
  
  par(def.par)
  graphics.off()
} # end stai loop



if(doSkill){
# make a plot of pctHabSkill

fn <- paste(plotPath,"pctHab_mod_meas.png",sep="")
width=3.5
height=4.0
png(filename = fn, width = width, height = height, units = "in", pointsize = 10, res=300)


plot(NA, NA, xlim=c(-100,100), ylim=c(-100,100), asp=1
     ,xlab="Model"
     ,ylab="Observed"
     ,main="Change in Pct. HAB by zone"
     )
abline(0,1, col="gray", lty=2)
abline(h=0, col="gray")
abline(v=0, col="gray")

dfm <- pctHabSkill[pctHabSkill$type=="Model",]
dfo <- pctHabSkill[pctHabSkill$type=="Satellite",]

fcdays <- unique(dfm$fcday)
fcdays <- fcdays[2:length(fcdays)]
polys <- polys3_areas$poly[polys3_areas$poly != "99"]
xx <- which(names(dfm) %in% polys) 
cols <- tim.colors(5)
cols <- two.colors(n=5, start="orange",end="red")
cols <- c("#fef0d9"
  ,"#fdd49e"
  ,"#fdbb84"
  ,"#fc8d59"
  ,"#e34a33"
  ,"#b30000")
for(fcday in fcdays){
  sim <- dfm[dfm$fcday==fcday,xx] - dfm[dfm$fcday==0,xx]
  obs <- dfo[dfo$fcday==fcday,xx] - dfo[dfo$fcday==0,xx]
  points(sim, obs, pch=paste(fcday), col=cols[fcday])
}

graphics.off()
} # if doSkill

print("habbplot complete")


