# Mark Rowe 6-30-2015
# this R script will be called by habtracker.R
# it will do a series of 1D vertical mixing simulations using ptraj
# at specified station locations and make plots

def.par <- par(no.readonly = TRUE) # save default, for resetting...

# run the next two lines to load up variables to run this script outside of habtracker.R
#    dryrun <- TRUE
#    if(dryrun){print("Doing a dry run")}
#    source("habtracker.R")

print(" load common data and functions:  habtracker_functions.R")
setwd(habtrackerDir)
source("habtracker_functions.R")

#===================================
#  Define the stations to be used for the 1D simulations
#===================================

stas1d <- read.csv(paste0(habtrackerDir,"stations_1d.csv"))
stas1d <- stas1d[stas1d$station %in%
  c("WE2","WE6","WE8","WE4","WE12","WE13","WE25","NDBC45005","DWE","SouthPass","WE14","WE15","45176")
                 ,] # select the stations to use

# # add additional locations to estimate the mixing depth, but these won't be plotted
stas1d2 <- data.frame(
  station=NA
  ,station_depth=NA
  ,lon=c(-83.3876520246101, -83.3018028330233, -83.2717556159679, -83.2159536414365, -83.1472742881671, -83.0657175561597, -82.9111890113035, -82.9197739304622, -83.1301044498498, -83.082887394477, -82.9541136070969, -82.9884532837316, -82.7309057089712, -82.5162827300043, -82.563499785377, -82.1428387466018, -82.0484046358563, -82.1514236657605, -83.070010015739, -83.1944913435398)
  ,lat=c(41.772499642295, 41.9012734296752, 41.8497639147231, 41.7639147231364, 41.6737730719702, 41.6608956932322, 41.5922163399628, 41.708112748605, 41.9356131063099, 41.8454714551438, 41.8454714551438, 41.9313206467306, 41.720990127343, 41.5106596079554, 41.6780655315496, 42.1201888682215, 41.8583488338818, 41.5879238803835, 41.785377021033, 41.6694806123909)
  ,node=NA
)
stas1d <- rbind(stas1d, stas1d2)
# 
# plotPolys(shorelineLL, xlim=xlimLL, ylim=ylimLL, projection=TRUE)
# points(stas1d$lon, stas1d$lat, pch=3, col="red")
# points(stas1d2$lon, stas1d2$lat, pch=3, col="blue")

#===================================
#  Acquire GLCFS currents (FVCOM) for 1D runs
#  Concatenate the nowcast-24hr, nowcast, and forecast files
#  into one netcdf file
#===================================

print( 'Acquire GLCFS (FVCOM) currents for 1D runs...')

ddd1 <- sprintf("%03.0f", as.numeric(ddd)-1)
fvcomnc1 <- paste0(fvcomdirnc,yyyy,ddd1,"00.nc")
fvcom1dfiles <- paste(c(fvcomnc1, fvcomnc, fvcomfc), collapse=" ")
print(paste("fvcom1dfiles = ",fvcom1dfiles))

# fvcom1d <- paste0(fvcomdirfc, "fvcom1d.nc")
fvcom1d <- paste0(tmpcursdir, "fvcom1d.nc")

if(!dryrun & !skipCur){
  file.remove(fvcom1d)
#  string1 <- paste("ncrcat -F", fvcom1dfiles, fvcom1d)
# MDR, 6-27-2016, list the variables to contatenate in case nowcast and forecast files have different variables
  string1 <- paste("ncrcat -F -vu,v,omega,zeta,km,temp,Times,nv,h,siglev,siglay,art1,a1u,a2u,aw0,awx,awy", fvcom1dfiles, fvcom1d)
  system(string1)
}else{
  print("Skipping creation of the 1D currents file")
}

#===================================
#  Do the 1D runs
#===================================
print( 'Running 1D ...')
ptm <- proc.time() # track run time

npdirs1d <- nrow(stas1d)
system(paste("echo ",npdirs1d," > npdirs1d.txt"))
# set up directories to run a series of parallel runs
pdirs1d <- c()
for(ro in 1:npdirs1d){
  pdirs1d <- c(pdirs1d, paste0("dir1d",ro))
}

if(!dryrun){
  ri=1
  for(ri in 1:length(pdirs1d)){
    
    pdir <- paste0(habtrackerDir, pdirs1d[ri],"/")
    if(file.exists(pdir) == FALSE){dir.create(pdir)}
    
    setwd(pdir)
    # clean up old files
    rmFiles <- Sys.glob(paste("job-1.*",sep=""))
    file.remove(rmFiles)
    rmFiles <- Sys.glob(paste("*.out",sep=""))
    file.remove(rmFiles)
    rmFiles <- Sys.glob(paste("*.nc",sep=""))
    file.remove(rmFiles)
    rmFiles <- Sys.glob(paste("*.end",sep=""))
    file.remove(rmFiles)
    rmFiles <- Sys.glob(paste("*.ini",sep=""))
    file.remove(rmFiles)
    rmFiles <- Sys.glob(paste("*.dat",sep=""))
    file.remove(rmFiles)
    
    # write the run parameter file
    fncur <- fvcom1d
    fnini <- "erie1d-ptraj.ini"
    source(paste0(habtrackerDir, "write_ptraj_runfile.R"))
    
    # write the ini file
    totparts <- 1000 # 1000
    maxdepth <- -1.0
    ttrelease <- 0.0
    pdepth <- seq(0,maxdepth, by=maxdepth/totparts)
    
     lat1d <- stas1d$lat[ri]
     lon1d <- stas1d$lon[ri]

    # # assign vertical buoyant velocity based on diameters measured by flowcam and empirical relation from 
    # # Nakamura 1993
    # df <- read.csv("/Users/rowe/Documents/HABproject/ToxicToledoFlowcam04082014/WE12at20%dilution4Aug2015.csv",skip=3)
    df <- read.csv(paste0(habtrackerDir, "WE12at20%dilution4Aug2015.csv"),skip=3)
    dia <- sample(df$Diameter..ESD., size=totparts, replace = TRUE)
    # #### buoyant velocity and colony diameter from Nakamura et al. 1993 Wat. Res. Vol. 27, No. 6, pp. 979-983, 1993 Fig. 3
    # # from DigitizeNakamuraFig.R :
    # wb1 <- 10^(log10(dia)*1.471-1.384) # 88*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 3-Aug-1990 sample
    # wb2 <- 10^(log10(dia)*1.1694-0.4058) # 36*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 18-Sep-1990 sample
    # revised regression 9-30-2015
    wb1 <- 10^(log10(dia)*1.4439 -0.9194) # 88*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 3-Aug-1990 sample
    #wb1 <- 10^(log10(dia)*1.471  -1.380 ) # 36*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 18-Sep-1990 sample
    
    ptsmat <- matrix(0.0, nrow=totparts, ncol=6)
    ptsmat[,1] <- seq(1,totparts)
    ptsmat[,2] <- lon1d
    ptsmat[,3] <- lat1d
    ptsmat[,4] <- pdepth[1:totparts]
    ptsmat[,5] <- ttrelease
    ptsmat[,6] <-  wb1*1e-6
    
    ptsmat1 <- sprintf(fmt="%12.0f %12.6f %12.6f %10.6f %10.6f %1.6E", ptsmat[,1], ptsmat[,2] , ptsmat[,3], ptsmat[,4], ptsmat[,5], ptsmat[,6] )
    #    ptsmat1 <- sprintf(fmt="%12.0f %12.6f %12.6f %10.2f %10.2f", ptsmat[,1], ptsmat[,2] , ptsmat[,3], ptsmat[,4], ptsmat[,5] )
    
    write.table(totparts, file = fnini, append = FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(ptsmat1, file = fnini, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    print(paste("totparts = ",totparts))
    print(paste(Sys.time()))

    setwd(habtrackerDir)
    
  } # end ri loop to set up parallel directories
  
  # write the job file
  # 5-16-2016: mdr, modified version to start multiple runs from one job script, on one node, using the normal slurm queue
  # time trials with 1000 particles, 1D runs, IRW   = 2, DTRW  = 10, WB1 =  -1.0, RWMETHOD = 3, SMOOTH = 2, DO1D = T
  # npdirs   run time, s
  # 1        6.2
  # 2        6.78
  # 4        6.9-7.0
  # 8        7.15-7.8
  # 16       7.7-9.4
  # 32       8.4-10 
  # 56       8.3-12.3
  

    write.table("#!/bin/bash", file = jobfile, append = FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(paste("#SBATCH -D",habtrackerDir_hpc), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    if(machineName == "bear"){write.table("#SBATCH --partition=ops", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)}
    write.table("#SBATCH -N 1", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table("#SBATCH --error=job-1.err", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table("#SBATCH --output=job-1.out", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)

    ri=1
    for(ri in 1:length(pdirs1d)){
#    for(ri in 1:1){
      
      pdir <- paste0(habtrackerDir, pdirs1d[ri],"/")
      if(file.exists(pdir) == FALSE){dir.create(pdir)}
      
      write.table(paste0("cd ",paste0(habtrackerDir_hpc, pdirs1d[ri],"/")), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
      write.table(paste("/cm/shared/apps/slurm/current/bin/srun -n 1 -N 1", exefile, casename," >",fnrunout,"&"), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)

    } # end ri loop
    
    write.table("wait", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    # start ptraj runs bear
    system("/cm/shared/apps/slurm/current/bin/sbatch job-1", wait=FALSE) 
  
  
  
  ###########################################################################
  # check to see if runs are finished, and don't proceed until all have finished.
  finished <- rep(0, npdirs1d)
  nfinished <- sum(finished)
  while(nfinished < npdirs1d){
    Sys.sleep(10)
    for(ri in 1:npdirs1d){
      pdir <- paste0(habtrackerDir, pdirs1d[ri],"/")
      print(pdir)
      setwd(pdir)
      string1 <- paste("tail -1",fnrunout)
      xx <- system(string1,intern=TRUE)
      if(length(xx) < 1){xx <- "nofile"}
      print(paste("ri = ",ri,xx))
      if(xx =="  TADA!"){
        finished[ri] <- 1
      }
      
    } # end ri
    nfinished <- sum(finished)
    print(paste("finished = ",paste(finished,collapse=" ")))
  } # end checking for run completion
  setwd(habtrackerDir)
  
  etime <- (proc.time()[3]-ptm[3])/3600
  print(paste("elapsed 1D run time =",round(etime,2),"hrs"))
  
} # end if(!dryrun)

#===================================
#  Plot the 1D runs
#===================================
if(!dryrun){
print( 'Plotting 1D ...')
ptm <- proc.time() # track run time

ri=1
stas2plot <- which(is.na(stas1d$station)==FALSE)
for(ri in stas2plot){
  station <- stas1d$station[ri]
  scenario <- station
  print(paste("Post-process 1D station ",station))
  
  pdir <- paste0(habtrackerDir, pdirs1d[ri],"/")
  
  ncfile <- fnout1d
  skip = 0 # hours to skip for spinup
  
  nc <- open.ncdf(paste(pdir,ncfile,sep=""), readunlim=FALSE)
  
  att1 <- att.get.ncdf(nc, 0, "dtout")
  dtout1 <- as.numeric(att1$value)/3600 #10/3600 # output interval, hours
  
  # read particle xy arrays
  start=c(1,skip+1) 
  OFFSET   = skip +1
  host <- get.var.ncdf(nc, "host", start=start, count=c(-1,-1))
  inactive <- which(host<=0, arr.ind=TRUE)
  host1 <- host[1,1]
  
  #  x <- get.var.ncdf(nc, "x", start=start, count=c(-1,-1)) - 360
  #  y <- get.var.ncdf(nc, "y", start=start, count=c(-1,-1))
  z <- get.var.ncdf(nc, "z", start=start, count=c(-1,-1))
  elev <- get.var.ncdf(nc, "elev", start=start, count=c(-1,-1))
  depth <- get.var.ncdf(nc, "depth", start=start, count=c(-1,-1))
  npts <- nrow(z)
  #z <- z[,1:500]
  z[inactive] <- NA
  z <- z-elev # adjust z so that the surface is at zero
  zmin <- min(z, na.rm=T)
  nt <- dim(z)[2]
  nt1 <- nt-1
  ztime <- seq(0,nt1)*dtout1
  nhr <- floor(nt*dtout1)
  nhr1 <- nhr-1
  
  # get temperature and kh profile from forcing file
  nc2 <- open.ncdf(fvcom1d, readunlim=FALSE)
  nv <- get.var.ncdf(nc2, "nv", start=c(1,1), count=c(-1,-1))
  h <- get.var.ncdf(nc2, "h", start=c(1), count=c(-1))
  Times <- get.var.ncdf(nc2, "Times", start=c(1,1), count=c(-1,-1))
  time2 <- as.POSIXct(Times, format="%Y-%m-%dT%H:%M:%S", tz="UTC")
  
#  node1 <- nv[host1,1]
  node1 <- stas1d$node[ri]
  if(is.na(node1)==TRUE){
    diff <- (coordsn$lon - stas1d$lon[ri])^2 + (coordsn$lat - stas1d$lat[ri])^2
    node1 <- coordsn$node[which(diff==min(diff))[1]]
  }
  temp1 <- get.var.ncdf(nc2, "temp", start=c(node1,1,OFFSET), count=c(1,-1,nhr1))
  h <- abs(min(z)) # get.var.ncdf(nc2, "h", start=c(node1), count=c(1)) 
  siglev <- get.var.ncdf(nc2, "siglev", start=c(node1,1), count=c(1,-1))
  siglay <- get.var.ncdf(nc2, "siglay", start=c(node1,1), count=c(1,-1))
  depths <- as.numeric(h)*siglay
  levels <- as.numeric(h)*siglay
  kh <- get.var.ncdf(nc2, "km", start=c(node1,1,OFFSET), count=c(1,length(depths),nhr1))
  
  
  # plot z vs time with raster
  fn <- paste(plotPath, "1D_",scenario, "_concProfile.png",sep="")
  width=5
  height=3
  # png(filename = fn, width = width, height = height, units = "in", pointsize = 10, res=300)
  png(filename = fn, width = 568, height = 448, units = "px", pointsize = 14)
  
  par(mar = c(8.5, 4, 2.5, 8) + 0.1) #c(bottom, left, top, right)
  times <- ztime # seq(1,ncol(x))
  nt <- length(times)
  nd <-  length(depths) #*2 Use 20 rows for ratio calculation 5-19-2015
  dmin <- min(levels) #zmin # ylim[1] # zmin
  dmax <- 0 #ylim[2]
  r1 <- raster(nrows=nd, ncol=nt, xmn=0, xmx=max(times), ymn=dmin, ymx=dmax)
  xy <- cbind(times, as.vector(t(z)))
  r2 <- rasterize(xy, r1, fun='count')
  H <- abs(dmin)
  r2 <- r2*nd/H # convert particles per vertical cell to particles per meter of depth
  r2[is.na(r2)] <- 0.0 # set concentration =0 where there are no particles
  r2 <- r2*H/npts # /mean(r2@data@values) # normalize concentration to mean
  #zlim=c(min(r2@data@values),max(r2@data@values))
  zlim <- c(0,5)
  r2[r2 > zlim[2]] <- zlim[2] # set max values to max of zlim to prevent overplotting
  cols <- tim.colors(256)
  
  imageTime <- as.POSIXct(paste(yyyy,ddd,hh), format="%Y %j %H", tz="UTC") 
  xxi <- which(time2 == imageTime)
  
  image(r2, col=cols
        ,ylab="Depth, m"
        ,xlab="" #Time, hr"
        ,xaxt="n"
        ,main=paste("Station",scenario)
        #, cex.main=0.8
        , las=2
        , xlim=c(xxi, ncol(r2)-1) # start the plot at the image time
        # ,xlim=c(200,245)
        # ,zlim=c(quantile(r2, probs=0.001),quantile(r2, probs=0.999))
        , zlim=zlim
  )
  
  hrs <- format(time2-3600*4, "%H")
  xx <- which(hrs %in% c("00","06","12","18")) # seq(1, nhr, 6)
  xx1 <- which(hrs %in% c("00")) # seq(1, nhr, 6)
  if(ncol(r2)/6 <= 20){
    labels <- format(time2[xx]-3600*4, "%H")
  }else{
    labels <- FALSE
  }
  axis(1, at=xx-1, labels=labels, las=2)
  axis(1, at=xx1-1, labels=format(time2[xx1]-3600*4, "%Y-%m-%d %H"), las=2)
  mtext("Normalized concentration", side=4, line=1)
  mtext("Time, EDT", side=1, line=7)
  

#  abline(v=xxi-1, lwd=2, lty=2)
#  mtext("Image time", side=3, at=xxi-1 , line=0.24, cex=0.8)
  
  #  logo.draw( logo1, 0.01, 0.01, 5.0)
  logo.draw( logo1, 0.67, 0.005, 6.0)
  
  image.plot(r2
             ,legend.mar=6
             # ,legend.lab="Normalized concentration"
             ,legend.only=TRUE
             ,add=TRUE
             # ,main=""
             # ,xlim=c(200,245)
             # ,zlim=c(quantile(r2, probs=0.001),quantile(r2, probs=0.999))
             , zlim=zlim
  )
  
  graphics.off()
  
  # plot temperature profile
  
  temp2 <- temp1 # get.var.ncdf(nc2, "temp", start=c(node1,1,OFFSET), count=c(1,-1,nt))
  
  rt <- raster(nrows=nrow(temp2), ncol=ncol(temp2), xmn=0, xmx=max(times), ymn=dmin, ymx=dmax)
  rt <- setValues(rt, temp2)
  
  fn <- paste(plotPath,"1D_",scenario, "_TempProfile.png",sep="")
  png(filename = fn, width = 568, height = 448, units = "px", pointsize = 14)
  par(mar = c(8.5, 4, 2.5, 8) + 0.1) #c(bottom, left, top, right)
  
  zlim=c(min(rt@data@values),max(rt@data@values))
  image(rt, col=cols
        ,ylab="Depth, m"
        ,xlab="" #Time, hr"
        ,xaxt="n"
        ,main=paste("Station",scenario)
        , xlim=c(xxi, ncol(rt)) # start the plot at the image time
        #, cex.main=0.8
        , las=2
        # ,xlim=c(200,245)
        # ,zlim=c(quantile(r2, probs=0.001),quantile(r2, probs=0.999))
        , zlim=zlim
  )
  
  hrs <- format(time2-3600*4, "%H")
  xx <- which(hrs %in% c("00","06","12","18")) # seq(1, nhr, 6)
  xx1 <- which(hrs %in% c("00")) # seq(1, nhr, 6)
  if(ncol(rt)/6 <= 20){
    labels <- format(time2[xx]-3600*4, "%H")
  }else{
    labels <- FALSE
  }
  axis(1, at=xx-1, labels=labels, las=2)
  axis(1, at=xx1-1, labels=format(time2[xx1]-3600*4, "%Y-%m-%d %H"), las=2)
  mtext(expression(paste("Temperature, ",degree,"C")), side=4, line=1)
  mtext("Time, EDT", side=1, line=7)
  
  imageTime <- as.POSIXct(paste(yyyy,ddd,hh), format="%Y %j %H", tz="UTC") 
  xxi <- which(time2 == imageTime)
#  abline(v=xxi-1, lwd=2, lty=2)
#  mtext("Image time", side=3, at=xxi-1 , line=0.25, cex=0.8)
  
  #  logo.draw( logo1, 0.01, 0.01, 5.0)
  logo.draw( logo1, 0.67, 0.005, 6.0)
  
  image.plot(rt
             ,legend.mar=6
             # ,legend.lab="Normalized concentration"
             ,legend.only=TRUE
             ,add=TRUE
             # ,main=""
             # ,xlim=c(200,245)
             # ,zlim=c(quantile(r2, probs=0.001),quantile(r2, probs=0.999))
             , zlim=zlim
  )
  
  graphics.off()
  
  # plot km profile
  
  temp2 <- kh # get.var.ncdf(nc2, "km", start=c(node1,1,OFFSET), count=c(1,20,nt))
  temp2 <- log10(temp2+1e-6)
  
  rt <- raster(nrows=nrow(temp2), ncol=ncol(temp2), xmn=0, xmx=max(times), ymn=dmin, ymx=dmax)
  rt <- setValues(rt, temp2)
  
  fn <- paste(plotPath,"1D_",scenario, "_KmProfile.png",sep="")
  png(filename = fn, width = 568, height = 448, units = "px", pointsize = 14)
  par(mar = c(8.5, 4, 2.5, 8) + 0.1) #c(bottom, left, top, right)
  
  zlim=c(min(rt@data@values),max(rt@data@values))
  image(rt, col=cols
        ,ylab="Depth, m"
        ,xlab="" #Time, hr"
        ,xaxt="n"
        , xlim=c(xxi, ncol(rt)) # start the plot at the image time
        ,main=paste("Station",scenario)
        #, cex.main=0.8
        , las=2
        # ,xlim=c(200,245)
        # ,zlim=c(quantile(r2, probs=0.001),quantile(r2, probs=0.999))
        , zlim=zlim
  )
  
  hrs <- format(time2-3600*4, "%H")
  xx <- which(hrs %in% c("00","06","12","18")) # seq(1, nhr, 6)
  xx1 <- which(hrs %in% c("00")) # seq(1, nhr, 6)
  if(ncol(rt)/6 <= 20){
    labels <- format(time2[xx]-3600*4, "%H")
  }else{
    labels <- FALSE
  }
  axis(1, at=xx-1, labels=labels, las=2)
  axis(1, at=xx1-1, labels=format(time2[xx1]-3600*4, "%Y-%m-%d %H"), las=2)
  mtext(expression(paste(Log[10],"(",K[m],", ",m^2," ",s^{-1},")")), side=4, line=1)
  mtext("Time, EDT", side=1, line=7)
  
  imageTime <- as.POSIXct(paste(yyyy,ddd,hh), format="%Y %j %H", tz="UTC") 
  xxi <- which(time2 == imageTime)
#  abline(v=xxi-1, lwd=2, lty=2)
#  mtext("Image time", side=3, at=xxi-1 , line=0.25, cex=0.8)
  
  #  logo.draw( logo1, 0.01, 0.01, 5.0)
  logo.draw( logo1, 0.67, 0.005, 6.0)
  
  image.plot(rt
             ,legend.mar=6
             # ,legend.lab="Normalized concentration"
             ,legend.only=TRUE
             ,add=TRUE
             # ,main=""
             # ,xlim=c(200,245)
             # ,zlim=c(quantile(r2, probs=0.001),quantile(r2, probs=0.999))
             , zlim=zlim
  )
  
  graphics.off()
  
  
  ###################################################################
  
} # end ri loop over 1D stations

# Plot a map showing the 1D station locations

# plot surface temperature at the image time on the map plot
values <- get.var.ncdf(nc2, "temp", start=c(1,1,xxi-1), count=c(-1,1,1))

stas1d2 <- projLcc(stas1d[stas2plot,], prj_new)
xlim = c(xlimLcc[1], max(stas1d2$X)+20)
ylim = c(0, 81.477247)
temp <- values

require(akima)
plotres <- 0.5 # 0.2 resolution of interpolation grid, km
xlength <- (xlim[2]-xlim[1])/plotres
ylength <- (ylim[2]-ylim[1])/plotres
tol <- 10
xx <- which(coordsn$X <= xlim[2]+tol & coordsn$X >= xlim[1]-tol & coordsn$Y <= ylim[2]+tol & coordsn$Y >= ylim[1]-tol)
x=coordsn$X[xx]
y=coordsn$Y[xx]
z=values[xx]
data.interp <- interp(x,y,z, duplicate="mean"
                      ,xo=seq(min(x), max(x), length = xlength)
                      ,yo=seq(min(y), max(y), length = ylength)
                      , extrap=FALSE
                      , linear=TRUE
)
zlim <-  c(min(data.interp$z, na.rm=TRUE),max(data.interp$z, na.rm=TRUE))
zdiff <- zlim[2]-zlim[1]
if(zdiff < 10){
zdiff2 <- (10-zdiff)/2
zlim[1] <- zlim[1]-zdiff2
zlim[2] <- zlim[2]+zdiff2
}

fn <- paste(plotPath, "1D_stationMap.png",sep="")

setupPlot(fn)
cols <- tim.colors(256)

plotPolys(shoreline, projection=TRUE
          # ,plt = c(0.11, 0.98, 0.12, 0.88)
          ,plt = c(0.01, 1, 0.15, 0.9)
          # ,col=colorbar[1] # fill the background with the low-HAB color in case there are gaps
          ,xaxt="n",yaxt="n",bty="n",ylab="",xlab=""
          ,xlim=xlim, ylim=ylim
)
#box("figure")
datePlot <- format(imageTime-4*3600, "%Y-%m-%d %H:%M")
mtext(paste(datePlot,"EDT"), side=3, line=1, cex=2.0)

image(data.interp, add=TRUE,  col=cols)

addPolys(mask2, col="tan")
plotCities(cex.cities=0.8)
scalebar(188, 15, lwd=4)

# plot the 1D station locations and labels
cex1 <- 0.8
#xxb <- which(stas1d2$X == max(stas1d2$X))
xxb <- which(stas1d2$station %in% c("WE6","45176"))
xxnb <- which(!(1:nrow(stas1d2) %in% xxb))
x <- stas1d2$X
y <- stas1d2$Y
x[which(stas1d2$station == "WE25")] <- x[which(stas1d2$station == "WE25")] +1 # shift plotting position of WE25 slightly
points(x, y, pch=3)
text(stas1d2$X[xxnb], stas1d2$Y[xxnb], stas1d2$station[xxnb], font=2, pos=4, cex=cex1)
text(stas1d2$X[xxb], stas1d2$Y[xxb], stas1d2$station[xxb], font=2, pos=2, cex=cex1)

logo.draw( logo1, 0.02, 0.02, 11.0)

plot(c(0,10),c(0,10),type="n", yaxt="n", xaxt="n", xlab=""
     , ylab="", bty="n")

image.plot(zlim=zlim
           ,col=cols
           #,breaks= seq(zlim[1],zlim[2])
           #,lab.breaks=lab.breaks
           , legend.width = 6
           # ,legend.args=list(cex=2.0)
           #  ,legend.mar=18 #15.1
           ,legend.mar=25 #15.1
           , legend.only=T, add=F)

label <- expression(paste("Temperature, ",degree,"C"))
text(7.5,5,label,srt=90, cex=1.5)

par(def.par)
graphics.off()

etime <- (proc.time()[3]-ptm[3])/3600
print(paste("elapsed 1D plotting time =",round(etime,2),"hrs"))

} # end if(!dryrun)


  
  
