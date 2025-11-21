# Mark Rowe 1-24-2017
# This R script will read an ensemble of HAB Tracker output netcdf files, 
# and produce a single file with cyanobacterial chlorophyll concentration averaged over the ensemble

require(ncdf4)

topDir <- "/Users/rowe/atoz/hpc/hab/hab_tracker_2016/"
setwd(topDir)

archiveDirs <- Sys.glob(paste0(topDir,"hab_tracker/archive/*"))

# get some parameters from the habtracker.log file
massperpart <- 1e+10 # 1e+8   # microgram of chl mass per particle #1e+10 paper
fntrknc <- 'erienc-ptraj.nc'
fntrkfc <- 'eriefc-ptraj.nc'
conczres=1.0 #! vertical resolution for concentration output, meters, REAL

for(di in 1:length(archiveDirs)){ # loop over habtracker run IDs
#for(di in 147:length(archiveDirs)){ # loop over habtracker run IDs
  
  archiveDir <- paste0(archiveDirs[di],"/")
  pdirs <- Sys.glob(paste0(archiveDir, "parallel*"))
  npdirs <- length(pdirs)
  
  #Get file start time
  fn <- paste0(archiveDir,"timestamp.txt")
  # print(fn)
  date1 <- read.table(fn )
  year1 <-  date1[,1]
  day1 <- date1[,2]
  hour1 <- date1[,3]
  date <- as.POSIXct(paste(year1,day1,hour1), format="%Y %j %H", tz="GMT")
  print(paste("date = ",date))
  
  archiveShare <- paste0(topDir, "archiveShare/")
  if(file.exists(archiveShare) == FALSE){dir.create(archiveShare)}
  
  
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
      pdir <- paste0(pdirs[ri],"/") # paste0(habtrackerDir, pdirs[ri],"/")
      fntrk1 <- paste0(pdir,fntrk)
      if(file.exists(fntrk1)){
        print(paste("reading",fntrk1))
        trkList[[ri]] <- nc_open(fntrk1, readunlim=FALSE)
      }else{
        print(paste("nc file does not exist, fntrk1=",fntrk1))
        break
      }
    } # end ri
    
    # read 3D concentration from CONCOUT ptraj output
    # conc is output in particles/ m3
    # loop over output files from parallel runs and average the conc output arrays
    nc <- trkList[[1]] 
    vx <- ncvar_get(nc, "vx", start=1, count=-1)
    vy <- ncvar_get(nc, "vy", start=1, count=-1)
    dtout <- as.numeric(ncatt_get(nc, 0)$dtout) # output interval, s
    conc <- ncvar_get(nc, "conc", start=c(1,1,1), count=c(-1,-1,-1))*0.0 #get array size
    count <- 0
    for(ri in 1:npdirs){
      nc <- trkList[[ri]]
      c1 <- ncvar_get(nc, "conc", start=c(1,1,1), count=c(-1,-1,-1))*massperpart/1000
      if(sum(dim(c1)) == sum(dim(conc))){ # if one of the files doesn't have full number of records, skip it
        conc <- conc + c1 
        count <- count + 1
      }
      nc_close(nc)
    }
    conc <- conc/count
    rm(trkList)
    
    # define a new netcdf file
    node <- ncdim_def( "node", "FVCOM grid node number", 1:dim(conc)[1])
    zlay <- ncdim_def( "zlay", "z-layer number", 1:dim(conc)[2])
    nrec <- dim(conc)[3]
    t <- ncdim_def( "time", "time record number", 1:nrec, unlim=TRUE)
    dimnchar <- ncdim_def("nchar", "", 1:25)
    
    if(ic==1){
      dates <- date + seq(0, nrec-1)*dtout
      nrecnc <- nrec
    }else{
      dates <- date + seq(nrecnc-1, nrec+nrecnc-2)*dtout
    }
    
    # Make a variable with those dimensions.  Note order: time is LAST
    timestamp <- ncvar_def("timestamp", "GMT", list(dimnchar, t), NULL, prec="char")
    lon <- ncvar_def("lon", "degreesE", list(node), NA)
    lat <- ncvar_def("lat", "degreesN", list(node), NA)
    concvar <- ncvar_def("conc",    "microgram/L cyanobacterial chlorophyll",  list(node,zlay,t), NA )
    
    # Create a netCDF file with this variable
    runID <- strsplit(archiveDir,"/")[[1]]
    runID <- runID[length(runID)]
    archiveShare1 <- paste0(archiveShare,runID,"/")
    if(file.exists(archiveShare1) == FALSE){dir.create(archiveShare1)}
    if(ic==1){fn1 <- "habtracker_nowcast.nc"}
    if(ic==2){fn1 <- "habtracker_forecast.nc"}
    ncnew <- nc_create( paste0(archiveShare1,fn1), list(lon, lat, timestamp, concvar))
    
    ncatt_put(ncnew, varid="zlay", attname="z-layer_thickness_m", attval=conczres, prec="float" )
    ncatt_put(ncnew, varid=0, attname="title", attval="Experimental Lake Erie HAB Tracker", prec="char" )
    ncatt_put(ncnew, varid=0, attname="institution", attval="University of Michigan CIGLR, NOAA GLERL", prec="char" )
    ncatt_put(ncnew, varid=0, attname="source", attval="https://www.glerl.noaa.gov/res/HABs_and_Hypoxia/habTracker.html", prec="char" )
    
    datechar <- format(dates, "%Y-%m-%d %H:%M:%S GMT")
    ncvar_put(ncnew, timestamp, vals=datechar)
    ncvar_put(ncnew, lon, vals=vx)
    ncvar_put(ncnew, lat, vals=vy)
    ncvar_put(ncnew, concvar, vals=conc)
    
    nc_close(ncnew)
    rm(conc)
    
  } # end ic
  
  webDir <- paste0(archiveShare1, "webfiles/")
  if(file.exists(webDir) == FALSE){dir.create(webDir)}
  
  cpFiles <- Sys.glob(paste0(archiveDir, "webfiles/*"))
  cpFiles1 <- Sys.glob(paste0(archiveDir, "webfiles/e*"))
  cpFiles <- cpFiles[!(cpFiles %in% cpFiles1)]
  file.copy(cpFiles, webDir)
  
  cpFiles <- Sys.glob(paste0(archiveDir, "*.tif"))
  file.copy(cpFiles, archiveShare1)
  
} # end di loop

# ncnew <- nc_open(paste0(archiveShare1,fn1))
# ncvar_get(ncnew, "timestamp", start=c(1,1), count=c(-1,-1))
# c1 <- ncvar_get(ncnew, "conc", start=c(1,1,1), count=c(-1,-1,-1))
# lat <- ncvar_get(ncnew, "lat", start=c(1), count=c(-1))
# lon <- ncvar_get(ncnew, "lon", start=c(1), count=c(-1))



