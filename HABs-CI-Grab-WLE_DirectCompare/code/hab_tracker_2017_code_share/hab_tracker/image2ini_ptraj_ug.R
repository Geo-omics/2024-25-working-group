# Mark Rowe 12-16-2014
# Create ptraj initial particle position file for habtracker
# this "ug" version will use the FVCOM unstructured grid rather than the 1km grid and raster package

#.libPaths("~/zebra2/users/rowe/rlib")
#require(raster)
#require(fields) # only for plotting
require(PBSmapping)
require(ncdf)
def.par <- par(no.readonly = TRUE) # save default, for resetting...
ptm1 <- proc.time()

# run the next two lines to load up variables to run this script outside of habtracker.R
# dryrun <- TRUE
# source("habtracker.R")

# # limits for map plots
# xlimLL <- c(-83.47336, -82.01775)
# ylimLL <- c(41.36373, 42.08565)

# set parameters
#massperpart <- 1e+6 # 1e+8   # microgram of chl mass per particle
#concdep <- 0.1 # assume depth in m of grid cell (only used for doIDL)
pdepth <- 0.0 # depth to release particles  (only used for doIDL)
ttrelease <- 0.0 # time to release particles, hours elapsed from start
vres <- conczres # vertical resolution to calculate concentration, meters (only used if not doIDL and not COLDSTART)

#setwd("~/zebra2/users/rowe/hab/hab_tracker/")

#imagePath='image' #;mdr use relative path
#archivePath <- "../hab_text/"
#plotPath <- "../Rimages/"
#posPath <- "./old_pos/"
#fvcomPath="../fvcom_netcdf/"

print("image2ini, load common data and functions:  habtracker_functions.R")
source("habtracker_functions.R")

# read all image file names from directory
imageFiles <- Sys.glob(paste0(imagedir, "*.txt")) # list.files(imagePath, "*.txt", full.names=TRUE)
#imageFiles <- "~/zebra2/users/rowe/hab/hab_text/M2008255154800.txt"

fi <- 1
#for(fi in 1:length(imageFiles)){ # not sure why one would want to process more than one file; it would overwrite
#for(fi in 1:1){
  
  imageFile <- imageFiles[fi]
  
  print(paste("read image file, ",imageFile,", System time = ",Sys.time()))
  
  nchar1 <- nchar(imageFile)
  yyyyddd <- substr(imageFile, start=nchar1-16, stop=nchar1-10)
  hh <- substr(imageFile, start=nchar1-9, stop=nchar1-8)
  
  if(substr(imageFile, start=nchar1-17, stop=nchar1-17) == "D"){
    dummyFile <- TRUE
    print("imageFile is a dummy file")
    chlData <- data.frame(lon=NA, lat=NA, chl=NA, temp=NA )
  }else{
    dummyFile <- FALSE

  
  chlData <- read.table(imageFile, header=FALSE)
  names(chlData) <- c("lon","lat","chl","temp")
  
  chlData <- processChlData(chlData)
  
  # calculate some quality stats
  nwater <- length(chlData$chl)
  ngood <- length(which(is.na(chlData$chl)==FALSE))
  nbad <- length(which(is.na(chlData$chl)==TRUE))
  print(paste("nwater = ",nwater))
  print(paste("ngood = ",ngood))
  print(paste("nbad = ",nbad))
  print(paste(Sys.time()))
  
  # remove NA values  
  chlData <- chlData[is.na(chlData$chl)==FALSE,]
  
  chlData <- bufferFun(chlData, buffer_rad, buffer_rad2, fill=TRUE)
  
#             plotPolys(shorelineLL, projection=TRUE, bg="wheat", xlim=xlimLL, ylim=ylimLL)
#             cols <- getColors(chlData$chl)
#             #points(chlData$lon, chlData$lat, col=cols$cols, pch=16, cex=0.5)
#             points(chlData$lon, chlData$lat, pch=15, cex=0.2, col=cols$cols)
#             #points(chlData$lon, chlData$lat, col=fields::color.scale(chlData$chl), pch=16, cex=0.5)
#            # abline(h=lat1)
#           #  abline(v=lon1)
  
  } # end if not dummyFile

# MDR, 6-28-2016 we want to run the model even if there is no image or no new data  
#  if(ngood < 1 & !dummyFile){
#    print("No good satellite chlorophyll data, quitting")
#    print(paste(Sys.time()))
#  }else{
    
#         plotPolys(shoreline, projection=TRUE, bg="wheat")
#         cols <- getColors(chlData$chl)
#         points(chlData$X, chlData$Y, col=cols$cols, pch=16, cex=0.5)  
    
    if(!(doIDL)){
    if(!(dummyFile & coldstartOnly)){ 
      ####################################################################
      # assign MODIS chl values to FVCOM unstructured grid tracer control elements 
      # by inverse distance weigted or nearest neighbor interpolation
      print("Assign chl values to FVCOM TCEs")
      print(paste(Sys.time()))
      
      # remove NA values  
      chlData <- chlData[is.na(chlData$chl)==FALSE,]
      
      # calculate the mean chl in each FVCOM grid TCE
      if(dummyFile){
        chl2 <- c()
      }else{
        chl2 <- aggregate(chl ~ PID, data=chlData, mean)
      }

      nr <- nrow(coordsn)
      chl1 <- rep(NA, nr)    
      chl1[chl2$PID] <- chl2$chl
      
#                      cols <- getColors(chl1)
#                      plotPolys(psTces, col=cols$cols, border=cols$cols, projection=T
#                             #  , xlim=xlimLcc
#                             #  , ylim=ylimLcc
#                                )
      
      pts <- data.frame(
        X = coordsn$X#/1000
        ,Y = coordsn$Y#/1000
        ,chl = chl1
      )
      
      if(doIDW){
        
        pts <- getIDW(pts, IDWdist=IDWdist)  
        chl1 <- pts$chl
 
      }else{
        # if not doing IDW, fill in missing values with zeros  
        chl1[is.na(chl1)] <- 0.0
        
      } # end if doIDW
      
      #    
      ptsdf <- data.frame(X=pts[,1],Y=pts[,2])
      ptsdf$chl <- chl1

      
#                        cols <- getColors(chl1)
#                       # cols$cols[!(1:length(cols$cols) %in% chlNodes)] <- "gray"
#                       plotPolys(shoreline, projection=T
#                        #         ,xlim = xlimLcc
#                        #         ,ylim = ylimLcc
#                                 )      
#                       addPolys(psTces, col=cols$cols, border=cols$cols, projection=T )
#                       #plotPolys(psTces, col=fields::color.scale(chl1), border=fields::color.scale(chl1), projection=T)
#                 
      ####################################################################
      # initialize 3D particle positions by either 
      # coldstart or hotstart
      
      # read in ptraj output from previous run to initialize cloud-covered areas
      print(paste("Start 3D ini",Sys.time()))
      
      if(file.exists(paste(posPath,"timestamp.txt",sep=""))){
        ts <- read.table(paste(posPath,"timestamp.txt",sep=""))
      }else{
        ts <- data.frame(V1=0,V2=0,V3=0)
        #  ts <- data.frame(V1=2011,V2=211,V3=18)
      }
      
      y1 <- as.numeric(substr(yyyyddd, 1,4))
      y2 <- ts[1,1]
      dy <- ((y2-y1)*365)
      d2 <-  ts[1,2]                           # initial day of the archived model run
      d1 <-  as.numeric(substr(yyyyddd, 5,8))  # initial day of the current image
      dd <- ((d1-d2))
      
      print(paste("Forecast start, ",ts[1,1],ts[1,2],ts[1,3]))
      print(paste("Image date,     ",y1,d1,hh))
      
      if(dy == 0 & dd > 0 & dd <= fclength/24 & coldstartOnly == FALSE ){ # only use the model info if it was within 5 day forecast period
        coldstart <- FALSE
      }else{
        coldstart <- TRUE # force coldstart if not within the nowcast/forecast period
      }
      
      print(paste("coldstart = ",coldstart))
      print(paste(Sys.time()))
      
      # calculate the record number to read from the forecast ncdf file that corresponds 
      # to the time of the current image
      # record 1 of the forecast run would be at midnight the day after the timestamp.txt
      tind <- (d1-d2-1)*24 + as.numeric(hh) + 1
      
      #===================================
      #  Interpolate 1D concentration profiles to initialize 3D run
      #===================================
      if(ini2D==FALSE & do1Druns){
        print( 'Interpolating 1D conc profiles ...')
       # ptm <- proc.time() # track run time
        
        # find the 1D run record that corresponds to the image time 
#         nc <- open.ncdf(fvcom1d, readunlim=FALSE)
#         time1d <- get.var.ncdf(nc, "Times", start=c(1,1), count=c(-1,-1))
        tind1d <- 24 + as.numeric(hh) 
        
        # Loop through 1D run output files
        #   calculate the sigma layer at which conc is 1/2 of surface conc
        ri=1
        stas1d$nlev1 <- 0 # this will be the number of levels within SML at the station
        for(ri in 1:npdirs1d){

          station <- stas1d$station[ri]
          scenario <- station
         # node1 <- stas1d$node[ri]
          
          pdir <- paste0(habtrackerDir, pdirs1d[ri],"/")
          
          ncfile <- fnout
          
          nc <- open.ncdf(paste(pdir,ncfile,sep=""), readunlim=FALSE)
          
          att1 <- att.get.ncdf(nc, 0, "dtout")
          dtout1 <- as.numeric(att1$value)/3600 #10/3600 # output interval, hours
          tind1d <- ceiling(tind1d/dtout1)
          
          conc1 <- get.var.ncdf(nc, "conc", start=c(1,1,tind1d), count=c(-1,-1,1))
          conc1 <- apply(conc1, 2, sum)
          if(conc1[1] > 0){
            stas1d$nlev1[ri] <- length(which(conc1/conc1[1] > 0.5))
          }else{
            # assume uniform conc profile if surface conc <= 0
            stas1d$nlev1[ri] <- length(conc1) 
          }
          
#           plot(conc1/conc1[1], -1*(1:length(conc1)), type="l", main=station)
#           abline(h = -length(which(conc1/conc1[1] > 0.5)))
#           ri = ri + 1
          
        } # end ri
        
        # interpolate SML thickness at the 1D stations to all the nodes
        # variable to be interpolated
        stas1d <- projLcc(stas1d, prj_new)
        nlev1 <- stas1d$nlev1 
        X <- stas1d$X
        Y <- stas1d$Y
        
        # interpolation grid
        pts <- coordsn[,c("X","Y")]
        npts <- nrow(coordsn)
        smlthick <- rep(1, npts)
        
        nn <- 1 # number of nearest points to use
        idp <- 1 # IDW exponent
        
        # IDW interpolation of the SML thickness
        ro <- 1
        for(ro in 1:npts){
          x1 <- pts[ro,1]
          y1 <- pts[ro,2]
          dist <- ((x1-X)^2 + (y1-Y)^2)^0.5
          #  ind1 <- which(dist==min(dist))
          o <- rank(dist)
          ind1 <- which(o <= nn)
          dist1 <- dist[ind1]
          nlev11 <- nlev1[ind1]
          weights <- dist1^(-idp) # IDW weights
          smlthick[ro] <- sum(nlev11*weights/sum(weights))
        }
        
        # make a plot of the sml thickness interpolated from 1D runs
        fn <- paste(plotPath, "sml_thickness_",yyyy,ddd,hh,".png",sep="")
        setupPlot(fn)
                 cols <- fields::color.scale( smlthick, col=tim.colors(256))
                 plotPolys(psTces, col=cols, border=cols, projection=T, xlab="km", ylab="km"
                           ,main="SML thickness, m"
                           ,xlim=xlimLcc, ylim=ylimLcc
                           )
                 addPolys(shoreline, density=0)
                 text(stas1d$X, stas1d$Y, stas1d$nlev1*conczres, col="gray")
        graphics.off()
        
      } # end if ini2D==FALSE & do1Druns
      
      
      if(coldstart == FALSE){ # read in model data
        print(paste("coldstart = ",coldstart,"reading in model data from ",posPath))
        
        #         trks <- read.table(paste(posPath,"erienc-ptraj.end",sep=""),skip=1,header=FALSE)
        #         names(trks) <- c("PID","lon","lat","depth","ttrelease")
        #         trks$lon <- trks$lon-360
        #         trks <- projLcc(trks, prj_new)
        
        # instead of using the end file from the nowcast run, use the forecast run at the time
        # of the current satellite image so that the simulated vertical distribution will be for
        # the current conditions.
        
        if(tind < 1){
          tind <- as.numeric(hh) - ts[1,3] + 1
          print("Warning: time index earlier than forecast run, using nowcast")
          #  nc <- open.ncdf(paste(posPath,"erienc-ptraj.nc",sep=""), readunlim=FALSE)
          fntrk <- fntrknc
          if(tind < 1){ tind = 1}
          print(paste("tind=",tind,"hh=",hh,"ts=",paste(ts[1,1],ts[1,2],ts[1,3])))
        }else{
          print("Opening forecast ncdf file")
          #  nc <- open.ncdf(paste(posPath,"eriefc-ptraj.nc",sep=""), readunlim=FALSE)  
          fntrk <- fntrkfc
        }
        
        trkList <- list() # empty list to hold parallel nc files
        for(ri in 1:npdirs){
          pdir <- paste0(posPath, pdirs[ri],"/")
          fntrk1 <- paste0(pdir,fntrk)
          if(file.exists(fntrk1)){
            print(paste("reading",fntrk1))
            trkList[[ri]] <- open.ncdf(fntrk1, readunlim=FALSE)
          }else{
            print(paste("nc file does not exist, fntrk1=",fntrk1))
            break
          }
        } # end ri
        
        nc <- trkList[[1]]
        nrec <- length(get.var.ncdf(nc, "time", start=c(1), count=c(-1)))
        if(tind > nrec){
          print(paste("Warning: time index later than forecast run, tind=",tind))
          tind <- nrec
        }
        zvals <- get.var.ncdf(nc, "zvals", start=c(1), count=c(-1))
        
        # Even in the case of parallel runs (npdirs > 1), only use the first 
        # particle positions file. The number of particles will be updated to match
        # the observed satellite surface conc
        trks <- data.frame(
          lon = get.var.ncdf(nc, "x", start=c(1,tind), count=c(-1,1)) - 360.0
          ,lat = get.var.ncdf(nc, "y", start=c(1,tind), count=c(-1,1))
          ,z1 = get.var.ncdf(nc, "z", start=c(1,tind), count=c(-1,1))
        )
        trks$EID <- 1:nrow(trks)
        
        # calculate sigma position and add to trks data frame
        depth1 = get.var.ncdf(nc, "depth", start=c(1,tind), count=c(-1,1))
        elev1 = get.var.ncdf(nc, "elev", start=c(1,tind), count=c(-1,1))
        trks$sigma <- (trks$z1- elev1) / (elev1 + depth1)
        
        # add projected coordinates to data frame
        trks <- projLcc(trks, prj_new)
        
        # remove inactive particles
        trks <- trks[trks$z1 != 999,]
        
        trksu <- trks # trksu will be the updated 3D particle distribution
        
        #         # read forecast concentration on FVCOM nodes and z layers
        #         conc <- get.var.ncdf(nc, "conc", start=c(1,1,tind), count=c(-1,-1,1))*massperpart/1000
        
        # read 3D concentration from CONCOUT ptraj output
        # conc is output in particles/ m3
        # loop over output files from parallel runs and average the conc output arrays
        conc <- get.var.ncdf(nc, "conc", start=c(1,1,tind), count=c(-1,-1,1))*0.0 #get array size
        for(ri in 1:npdirs){
          nc1 <- trkList[[ri]]
          conc <- conc + get.var.ncdf(nc1, "conc", start=c(1,1,tind), count=c(-1,-1,1))*massperpart/1000
        }
        conc <- conc/ri
        
        #vres <- zvals[1] - zvals[2]
        
        # Limit the model concentration to a max value to avoid accumulating too many particles
        cmax <- 300
        conc[which(conc > cmax, arr.ind=TRUE)] <- cmax
        
        # read time step and output interval from output file
        biofn <- paste(posPath,pdirs[1],"/ptrajfc.out",sep="")
        for(ro in 0:20){
          line1 <- scan(biofn, what=character(), skip=ro, nlines=1, quiet=TRUE)
          if(length(line1)>2){
            if(line1[2] == "DTI="){dti <- as.numeric(line1[3])}
            if(line1[2] == "INSTP="){instp <- as.numeric(line1[3])}
          }
        }
        nsteps <- instp/dti # number of time steps averaged for conc output
        #         
        #                    cols <- getColors(conc[,1])
        #                   # cols$cols[!(1:length(cols$cols) %in% chlNodes)] <- "gray"
        #                    plotPolys(psTces, col=cols$cols, border=cols$cols, projection=T)
        #                   #plotPolys(psTces, col=fields::color.scale(chl1), border=fields::color.scale(chl1), projection=T)
        # 
        #         # plot a conc profile by clicking on the map
        #         xy <- locator(1)      
        #         pts <- matrix(c(coordsn$X, coordsn$Y), nrow=nrow(coordsn), ncol=2)
        #         pt <- c(xy$x*1000, xy$y*1000)
        #         dist <- spDistsN1(pts, pt, longlat = FALSE)
        #         node1 <- coordsn$node[which(dist == min(dist))]
        #         plot(conc[node1,],zvals-vres/2)
        
      }else{ # coldstart = TRUE, initialize from current MODIS image
        
        chl1[is.na(chl1)] <- 0.0
        chl2 <- data.frame(
          PID = 1:length(chl1)
          ,chl = chl1
        )
        
      } # end if coldstart == FALSE
      
#                             cols <- getColors(chl2$chl)
#                              plotPolys(shoreline, projection=T
#                                        ,xlim = xlimLcc
#                                        ,ylim = ylimLcc
#                                        )
#                             points(coordsn$X, coordsn$Y, col=cols$cols, pch=20)
      
      # need to know the surface elevation at each node
      # we can get this from ptraj output if particles are present
      # but if there are no particles at a node, we need to get it from FVCOM
      # ncCur <- open.ncdf(paste(fvcomPath,yyyyddd,"00.nc",sep=""), readunlim=FALSE)
      if(!dummyFile){
        print(paste("reading zeta from ",fvcomnc))
       ncCur <- open.ncdf(fvcomnc, readunlim=FALSE)
       tindCur <- as.numeric(hh) 
       zeta = get.var.ncdf(ncCur, "zeta", start=c(1,tindCur), count=c(-1,1))
      }else{
       zeta = rep(0.0, nrow(coordsn))
      }

      # Loop through FVCOM nodes that have current satellite data
      #     calculate concentration in each layer
      #     determine the depth at which conc differs from the surface
      #     calculate how many particles are needed as f(satellite conc, volume of water)
      #     delete excess particles or
      #     create new particles to update concentration to that depth from satellite
      # End loop through nodes
      # 
      
      # to speed up model simulations, only initialize particles west of specified X value
      Xcrit <- xlimLcc[2] # 120000
      xx <- which(coordsn$X < Xcrit)
      chl2 <- chl2[chl2$PID %in% xx,]
      if(coldstart==FALSE){
        xx <- which(coordsn$X >= Xcrit)
        conc[xx,] <- 0
        trks <- trks[trks$X < Xcrit,]
        trksu <- trks
      }
      
      chlNodes <- chl2$PID # nodes that have current satellite data
      
#                  cols <- getColors(chl2$chl)
#                  plotPolys(shoreline, projection=T
#                             ,xlim = xlimLcc
#                             ,ylim = ylimLcc
#                              )
#                  points(coordsn$X[chlNodes],coordsn$Y[chlNodes], col=cols$cols, pch=20)
#  
      
      if(length(chlNodes) > 0){
        # estimate the max number of new particles needed, to pre-allocate new particle array
        nparts=sum(ceiling(chl2$chl*art1[chlNodes]*(coordsn$depth[chlNodes]+zeta[chlNodes])/0.001/massperpart))
        newparts <- matrix(NA, nrow=nparts*1.5, ncol=4)
        nnn <- nrow(newparts)
        rowrite <- 1 # row to write in newparts matrix
        
        print("loop through nodes and assign initial particle locations")
        m <- chl2$PID[which(chl2$chl == max(chl2$chl))][1]
        for(m in chlNodes){
#                               mi <- 1
#                               mi <- mi + 10
#                               m <- chlNodes[mi]
#           
          # satellite surface chl at node m
          surfChl <- chl2[chl2$PID == m,"chl"]
          if(surfChl < minChl){surfChl <- 0}
          
          # find the points that are within TCE m
          poly1 <- psTces[psTces$PID==m,] # polygon of tracer control element m
          Xmin <- min(poly1$X)
          Xmax <- max(poly1$X)
          Ymin <- min(poly1$Y)
          Ymax <- max(poly1$Y)
          
          # vertical dimensions
          if(coldstart){
            dmin <- -coordsn$depth[m] #-0.1
            dmax <- zeta[m] #+ 0.1
            brks <- seq(dmax, dmin, -vres)
            brks[length(brks)] <- dmin
            brks <- unique(brks)
            nlay <- length(brks) - 1 # number of vertical z layers at this node
          }else{ # hotstart
            fp <- findPolys(events=trks, polys=poly1, maxRows = 1e+07)
            trks1 <- trks[trks$EID %in% fp$EID,]
            dmin <- min(c(-coordsn$depth[m], trks1$z1)) #-0.1
            dmax <- max(c(zeta[m], trks1$z1)) #+ 0.1
            brks <- seq(dmax, dmin, -vres)
            brks[length(brks)] <- dmin
            brks <- unique(brks)
            nlay <- min((length(brks) - 1),length(zvals)) # number of vertical z layers at this node
          }
          
          # sml is a vector with one value per z-layer
          # sml has a value of 1 if the layer is in the SML, 0 otherwise
          # layers with sml=1 will be updated to match the satellite surface conc
          sml <- rep(1, nlay)
          
          # determine depth of the surface mixed layer if hotstart
          # MODIS surface concentration will be applied over sml
          if(coldstart == FALSE & nlay > 1 ){          
            
            # remove the particles in TCE m from the updated particle list
            # while we consider their fate
            trksu <- trksu[!(trksu$EID %in% fp$EID),]  
            
#                                   plotPolys(shoreline, projection=TRUE, main=paste("Node =",m))
#                                   cols <- getColors(surfChl)
#                                   addPolys(poly1, col=cols$cols, border=cols$cols, projection=T)
#                                   
#                                   cols <- getColors(surfChl)
#                                   #   plotPolys(poly1, col=cols$cols, border=cols$cols, projection=T)
#                                   plotPolys(poly1, projection=T, main=paste("Node =",m))
#                                   cols <- fields::color.scale(abs(trks1$z1))
#                                   points(trks1$X, trks1$Y, col=cols, pch=16)
            
            # calculate number of forecast particles in each layer from particle positions (trks1)
            trks1$dcut <- cut(trks1$z1, breaks=brks)
            c1trks1 <- rev(as.numeric(table(trks1$dcut)))
            
            # or, calculate number of forecast particles in each layer from conc output
            c1 <- conc[m,]/massperpart*1000*vres*art1[m]*nsteps # number of particles counted per layer over output interval
            
            if(ini2D==FALSE){
              
              # use the sml thickness interpolated from 1D runs 
              sml[1:length(sml)] <- 0
              nlev11 <- min(c(smlthick[m],length(sml)))
              sml[1:nlev11] <- 1
              
#              # assign sml thickness from previous 3D run 
#               c1u <- (c1)^-0.5 * c1 # uncertainty of the particle count (Sawford 1985 )
#               c1u[c1 == 0] <- 2  # give a finite uncertainty to zero count
#               c1u <- c1u + 0.2*c1 # assign a minimum relative uncertainty to the concentration
#               
#               # xx are the levels differ significantly from the surface value
#               xx <- c(which(c1-c1u > c1[1]+c1u[1]), which(c1+c1u < c1[1]-c1u[1]))
#               if(length(xx)==0){xx=1}
#               
#               firstDiff <- xx[1] # first level that differs significantly from surface
#               if(length(firstDiff) > 0){
#                 sml[firstDiff[1]:nlay] <- 0
#               }
              
#                                     hist(trks1$z1, breaks=brks, xlab="Depth, m", ylab="Particle count", main=paste("Node =",m), freq=TRUE)
#                                     points(brks[1:length(zvals)]-vres/2, c1*dti/instp)
#                                     plot(c1, ylim = c(0, max(c1+c1u)), xlab="Layer", ylab="Particle count", main=paste("Node =",m))
#                                     points(c1-c1u, pch=6)
#                                     points(c1+c1u, pch=2)
#                                     abline(h=c1[1]-c1u[1])
#                                     abline(h=c1[1]+c1u[1])
#                                     abline(v=seq(1,nlay)*sml)
            } # end if ini2D==FALSE
            
          } # end if coldstart==FALSE
            
            # if the following line is commented, particles will be applied uniformly to all layers 
            if(ini2D & nlay > 1){
              sml[2:nlay] <- 0 # only apply particles to first layer 
            }
            
          ######### update particles in sml to match satellite concentration
          # place at random XY within TCE
          
          k=length(sml[sml == 1])

## MDR, 8-3-2015 don't carry particles forward at nodes with satellite data
## instead, wipe them out and initialize new particles.
            np1 <- 0
#           if(coldstart){
#             np1 <- 0
#           }else{ # hotstart
#            # pctDiff <- (conc[m,1]-surfChl)/(surfChl)*100
#            # if(pctDiff > 100 & ini2D == FALSE){k=nlay} #if model is much higher than satellite, apply satellite conc to whole water column
#             np1 <- sum(c1trks1[1:k]) # current number of particles in control volume
#           }
          
          d1 <- brks[1]
          d2 <- brks[k+1]
          
          # calculate number of particles to place in the TCE over sml to match 
          # satellite chl concentration
         # nparts=ceiling(surfChl*art1[m]*(d1-d2)/0.001/massperpart)
          nparts=round(surfChl*art1[m]*(d1-d2)/0.001/massperpart, 0)
          
          npDiff <- round((nparts - np1),0) # number of particles needed
                            #  npDiff

## MDR, 8-3-2015 don't carry particles forward at nodes with satellite data
#           if(coldstart == FALSE){
#             if(npDiff < 0){
#               # delete excess particles
#               xx <- which(trks1$z1 <= d1 & trks1$z1 > d2,)
#               npDiff1 <- min(abs(npDiff), length(xx))
#               xx <- xx[1:abs(npDiff1)]
#               trks1 <- trks1[!(1:nrow(trks1) %in% xx),]
#             } # end if npDiff < 0
#             
#             #if(ini2D){
#               # don't carry forward particles below the surface layer
#               xx <- which(trks1$z1 <= d2 ,)
#               trks1 <- trks1[!(1:nrow(trks1) %in% xx),]
#             #}
#             
#             # put the particles from the TCE back into the updated particle list
#             # less the ones that got deleted
#             trksu <- rbind(trksu, trks1[,names(trksu)])
#           } # end if coldstart == FALSE
          
          if(npDiff > 0){
            # create additional particles
            # randomly placed within TCE and layer
            
            # candidate locations
            np <- 100
            locations1 <- data.frame(
              X = runif(np, min=Xmin, max=Xmax)
              ,Y = runif(np, min=Ymin, max=Ymax)
              ,EID = 1:np               
            )  
            # select candidate locations within TCE
            fp <- findPolys(events=locations1, polys=poly1, maxRows = 1e+07)
            locations1 <- locations1[locations1$EID %in% fp$EID,]
            
            # place at random XY locations known to be within the TCE
            xx <- sample(1:nrow(locations1), npDiff, replace=TRUE )
            x1 = locations1$X[xx]
            y1 = locations1$Y[xx]
            
#             # place particles at random z within SML
#             z1 = runif(npDiff, min=d2, max=d1)
            
            # place particles at the center z of each layer, then randomly place the extras
            mids <- brks[1:k] - vres/2 # midpoints of the layers within the sml
            npDiff1 <- floor(npDiff/k) # number of particles in each layer
            z1 <- c(rep(mids, npDiff1), sample(mids, npDiff-(npDiff1*k)))
            
            xx <- seq(rowrite,rowrite+npDiff-1)
            
            if(max(xx) > nnn){ # need to add rows to the matrix
              additionalparts <- matrix(NA, nrow=length(xx)+1000, ncol=3)
              newparts <- rbind(newparts, additionalparts)
              nnn <- nrow(newparts)
            }
            
            newparts[xx,1] <- x1
            newparts[xx,2] <- y1
            newparts[xx,3] <- z1  # input in z coordinates
            newparts[xx,4] <- (z1-dmax)/(dmax-dmin) # input in sigma 
            rowrite <- max(xx) + 1
            
            #             plotPolys(poly1, projection=T, main=paste("Node =",m," Depth =",k,"Chl =",round(surfChl, 2)))
            #             cols <- fields::color.scale(abs(z1))
            #             points(x1, y1, col=cols, pch=16)
            
         #   print(paste("Node =",m,k,rowrite,nrow(newparts),Sys.time()))
        #    plot(rep(1,length(z1)), z1, ylim=c(min(c(dmin,min(z1))),max(c(dmax,max(z1)))))
        #    abline(h=c(dmin, dmax))
            
          } # end if npDiff
          
  
        } # end loop over chl nodes
        
        # finalize the list of 3D points
        newparts <- newparts[is.na(newparts[,1]) == FALSE,]
        ptsdf1 <- data.frame(
          X = newparts[,1]
          ,Y = newparts[,2]
          ,z1 = newparts[,3]
          ,sigma = newparts[,4]
        )
        ptsdf1 <- projLatLon(ptsdf1, prj_new)
        
        if(coldstart == FALSE){ # carry forward particles in the nodes that were not updated
          trksu <- trksu[is.na(trksu[,1]) == FALSE,]
          ptsdf1 <- rbind(ptsdf1, trksu[,names(ptsdf1)])
        }
        
      }else{ # if length(chlNodes)==0, this would only happen if coldstart==FALSE
        ptsdf1 <- trksu
      }    
      
      # only initialize particles west of the specified easternLimit
      ptsdf1 <- ptsdf1[ptsdf1$lon <= easternLimit,]
      
    }else{ #  (dummyFile & coldstartOnly)
    # there are no chlData and coldstart=TRUE, so initialize two particles within the model domain
    
    print("Dummy file, initialize two particles")
    nodes1 <- sample(coordsn$node[coordsn$lon < easternLimit], size=2)
    ptsdf1 <- data.frame(
      lon = coordsn$lon[nodes1]
      ,lat = coordsn$lat[nodes1]
      ,sigma = rep(0, length(nodes1))
    )

    } # end if !(dummyFile & coldstartOnly)
      
      # make sure there are at least two particles to start
      if(nrow(ptsdf1) < 2){
        print("No particles to start, initialize two particles")
        nodes1 <- sample(coordsn$node[coordsn$lon < easternLimit], size=2)
        ptsdf1 <- data.frame(
          lon = coordsn$lon[nodes1]
          ,lat = coordsn$lat[nodes1]
          ,sigma = rep(0, length(nodes1))
        )
        ptsdf1 <- projxy(ptsdf1, prj_new=prj_new)
        ptsdf1$X <- ptsdf1$X/1000
        ptsdf1$Y <- ptsdf1$Y/1000
        ptsdf1$z1 <- 0
      }

    } # end if !(doIDL)
    
    ####################################################################
    # make a few plots to visualize initial particle distribution
    if(!dummyFile){
      
      date <- as.POSIXct(paste(substr(yyyyddd,1,4),substr(yyyyddd,5,7),hh), format="%Y %j %H", tz="GMT")
      tzoffset <- -4*3600 # time zone offset from GMT for plot time stamp, seconds
      tzstring <- "EDT"   # name of the time zone for plot time stamp
      datePlot <- paste(format(date+tzoffset, "%Y-%m-%d %H:%M"), tzstring)
      
      ploti=1
      if(coldstart){
        nplots <- 3
      }else{
        nplots <- 4
      }
      for(ploti in 1:nplots){ 
        
        if(ploti == 1){
          print(" make a plot to show the initial condition surface concentration")
          fn <- paste(plotPath, "initial_condition.png",sep="")
          
          # calculate concentration for the plot
          ptsdf2 <- ptsdf1
          ptsdf2$EID <- 1:nrow(ptsdf2)
          
          fp <- findPolys(events=ptsdf2, polys=psTces, maxRows = 1e+07)
          ptsdf2 <- merge(ptsdf2, fp[,c("EID","PID")])
          if(nrow(ptsdf2) < 1){next}
          # PID is the FVCOM node associated with the satellite obs
          ptsdf2$zeta <- zeta[ptsdf2$PID]
          ptsdf2$z2 <- ptsdf2$z1-ptsdf2$zeta # adjust z so surface is at zero
          
          if(plot2D == FALSE){ # plot surface concentration, otherwise assume all particles in surface layer
            ptsdf2 <- ptsdf2[ptsdf2$z2 >= -conczres,]
            if(nrow(ptsdf2) < 1){next}
          }
          
          ptsdf2$count <- 1
          counts <- aggregate(count ~ PID, data=ptsdf2, sum)
          values <- rep(0,nrow(coordsn))
          xx <- counts$PID
          values[xx] <- counts[,"count"]*massperpart/1000/art1[xx]/conczres
          
        }
        
        if(ploti == 2){
          print(" make a plot of the satellite surface conc")
          fn <- paste(plotPath, "satellite_surface_conc.png",sep="")
          
          chlData <- read.table(imageFile, header=FALSE)
          names(chlData) <- c("lon","lat","chl","temp")
          
          easternLimit <- 0 # plot all the data
          chlData <- processChlData(chlData)
          values <- chlData$chl
          
        }
        
        if(ploti == 3){
          print(" make a plot of the initial particle distribution")
          fn <- paste(plotPath, "initial_particle_distribution.png",sep="")
        }
        
        if(ploti == 4){
          print(" make a plot of the model output carried forward")
          fn <- paste(plotPath, "previous_run_surface_conc.png",sep="")
          
          # calculate concentration for the plot
          ptsdf2 <- trks
          ptsdf2$EID <- 1:nrow(ptsdf2)
          
          fp <- findPolys(events=ptsdf2, polys=psTces, maxRows = 1e+07)
          ptsdf2 <- merge(ptsdf2, fp[,c("EID","PID")])
          # PID is the FVCOM node associated with the satellite obs
          ptsdf2$zeta <- zeta[ptsdf2$PID]
          ptsdf2$z2 <- ptsdf2$z1-ptsdf2$zeta # adjust z so surface is at zero
          
          if(plot2D == FALSE){ # plot surface concentration, otherwise assume all particles in surface layer
            ptsdf2 <- ptsdf2[ptsdf2$z2 >= -conczres,]
          }
          
          ptsdf2$count <- 1
          counts <- aggregate(count ~ PID, data=ptsdf2, sum)
          values <- rep(0,nrow(coordsn))
          xx <- counts$PID
          values[xx] <- counts[,"count"]*massperpart/1000/art1[xx]/conczres
          
        }
        
        width=4.5
        height=3.5
        png(filename = fn, width = width, height = height, units = "in", pointsize = 10, res=300)
        
        par(mar=c(2, 0, 2, 0) + 0.1)#c(5, 4, 4, 2) + 0.1 c(bottom, left, top, right)
        mat <- matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE)
        widths <- c(6,1.5)
        widths <- widths/sum(widths)
        heights <- c(1)
        heights <- heights/sum(heights)
        n <- layout(mat,widths=widths,heights=heights)
        #  layout.show(n)
        
        values[values < 0.1] = 0.1 # don't plot zero values as no-data 
        
        xlim = xlimLcc # c(min(xx$x),max(xx$x))
        ylim = ylimLcc # c(min(xx$y),max(xx$y))
        
        temp <- values
        
        cols <- getColors(temp)
        
        brk <- cols$brk
        brkm <- max(brk)
        label <- zlab
        brkLabs <- round(brk, 0)
        
        nlevel <- length(brk)-1
        zlim <- cols$zlim
        colorbar <- cols$colorbar
        
        plotPolys(shoreline, main=datePlot, # projection=TRUE
                  #   ,col=colorbar[1] # fill the background with the low-HAB color in case there are gaps
                  ,xaxt="n",yaxt="n",bty="n",ylab="",xlab=""
                  ,xlim=xlim, ylim=ylim
        )
        
        if(ploti==3){
          # just plot the points

          addPolys(psTces, density=0, lwd=0.5)
          points(ptsdf2$X, ptsdf2$Y, pch=20, cex=0.1, col="red")
          }
        
        if(ploti %in% c(1,4)){
          
          zoom=TRUE
          if(zoom ){ #zoom into an area by clicking on the plot
            # you have to make the plot first, and don't run par(def.par) after making the plot
            # xx <- locator(2) # click upper left, and lower right of zoom region
            psClip <- clipPolys(psTces, xlim=xlim, ylim=ylim)
            temp <- values[unique(psClip$PID)]   
          }else{
            xlim=NULL
            ylim=NULL
            temp <- values
            psClip <- psTces
          }
          cols <- getColors(temp)
          addPolys(psClip, col=cols$cols, border=cols$cols, bg="gray",lwd=0.8) #plot TCE polygons
          
#           require(akima)
#           plotres <- 0.2 # resolution of interpolation grid, km
#           xlength <- (xlim[2]-xlim[1])/plotres
#           ylength <- (ylim[2]-ylim[1])/plotres
#           xx <- which(coordsn$X <= xlim[2] & coordsn$X >= xlim[1] & coordsn$Y <= ylim[2] & coordsn$Y >= ylim[1])
#           x=coordsn$X[xx]
#           y=coordsn$Y[xx]
#           z=values[xx]
#           data.interp <- interp(x,y,z, duplicate="mean"
#                                 ,xo=seq(min(x), max(x), length = xlength)
#                                 ,yo=seq(min(y), max(y), length = ylength)
#                                 , extrap=FALSE
#                                 , linear=TRUE
#           )
#           brk1 <- brk
#           brk1[length(brk1)] <- 1e6
#           image(data.interp, add=TRUE,  col=cols$colorbar, breaks=brk1)
        }
        
        if(ploti==2){
          points(chlData$X, chlData$Y, pch=15, cex=0.35, col=cols$cols)
        }
        
        if(ploti == 1){ # plot empty squares at satellite data locations on initial condition plot
        #  points(chlData$X, chlData$Y, pch=22, cex=0.7)
          addPolys(psTces[psTces$PID %in% chlNodes,], density=0, lwd=0.5)
          
#           # identify a node by clicking on the map plot
#           xy <- locator(1)
#           dist <- sqrt((coordsn$X - xy$x)^2 + (coordsn$Y - xy$y)^2)
#           m <- which(dist==min(dist))
        }
        
        if(ploti != 3){
        addPolys(mask2, col="tan")
        plotCities(cex.cities=0.6)
        
        plot(c(0,10),c(0,10),type="n", yaxt="n", xaxt="n", xlab=""
             , ylab="", bty="n")
        #box()
        text(1,5,label,srt=90)
        
        brk[1] <- 0
        image.plot(zlim=zlim
                   ,col=colorbar 
                   ,breaks= seq(zlim[1],zlim[2])
                   ,lab.breaks=format(brkLabs)
                   , legend.width = 6
                   ,legend.mar=18 #15.1
                   , legend.only=T, add=F)
        } # end if ploti != 3
        
        par(def.par)
        graphics.off()
        
      } # end ploti loop
    } # end if making initial conditions plot
    
    if(FALSE){
      # plot an X-Z slice through highest chl conc. with raster
      require(raster)
      require(fields)
      graphics.off()
      
      fn <- paste( getwd(),"/p3d_plot/tmp/concProfile_",yyyyddd,".png",sep="")
      width=5
      height=5
      png(filename = fn, width = width, height = height, units = "in", pointsize = 10, res=300)
      
      
      thick <- 2 # thickness of slice, km
      long <- 20 # length of slice, km
      
      node1 <- chl2$PID[which(chl2$chl == max(chl2$chl))][1]
      if(node1 %in% coordsn$node){
        
        xc <- coordsn[node1,"X"]
        yc <- coordsn[node1,"Y"]
        xmin <- floor(xc - long/2)
        xmax <- ceiling(xc + long/2)
        ymin <- yc - thick/2
        ymax <- yc + thick/2
        plotpts <- ptsdf1[ptsdf1$Y > ymin & ptsdf1$Y < ymax,]
        plotpts <- plotpts[plotpts$X > xmin & plotpts$X < xmax,]
        zmin <- min(plotpts$z1)
        zmax <- max(plotpts$z1)
        dmin <- floor(zmin)
        dmax <- ceiling(zmax)
        
        times <- seq(xmin, xmax, 1)
        depths <- seq(dmin, dmax, 1)
        nt <- length(times)
        nd <- length(depths)
        
        r1 <- raster(nrows=nd, ncol=nt, xmn=xmin, xmx=xmax, ymn=dmin, ymx=dmax)
        xy <- cbind(plotpts$X, as.vector(t(plotpts$z1)))
        r2 <- rasterize(xy, r1, fun='count')
        cv <- (dmax-dmin)/nd * (xmax - xmin)/nt*1000 * thick*1000 # cell volum, m3
        r2 <- r2/cv*massperpart/1000 # convert particles per cell to ugchl/L
        # r2[is.na(r2)] <- 0.0 # set concentration =0 where there are no particles
        image.plot(r2
                   ,ylab="Depth, m"
                   ,xlab="X distance, km: X-Z section through initial max chl patch"
                   ,main=expression(paste("Chlorophyll-a ",mu,g," ",L^{-1}))
                   # ,xlim=c(200,245)
                   ,zlim=c(0,r2@data@max)
        )
        
        graphics.off()
      } # end if node1 in coordsn
    } # END IF FALSE
    
    
    ####################################################################
    if(doIDL){ # initialize as in image2ini.pro IDL script
      print("Initialize as in IDL script")
      #     ;=============================================================
      #       ;	Translate Chl ASCII to particles (*.ini file)
      #     ;=============================================================
      
      ptsdf <- chlData[is.na(chlData$chl) == FALSE & is.na(chlData$temp) == FALSE,]
      # plot(ptsdf$lon, ptsdf$lat, pch=16, col=getColors(ptsdf$chl)$cols)
      # plot(ptsdf$lon, ptsdf$lat, pch=16, col=getColors(ptsdf$temp)$cols)
      ptsdf <- ptsdf[ptsdf$chl >= chlCrit & ptsdf$temp >= tCrit,]
      nhab <- nrow(ptsdf)
      print(paste("nHAB = ",nhab))
      
      xres <- 835.59484 # x and y distance, meters, represented by the satellite point concentration estimate
      yres <- 1120.33
      
      nr <- nrow(ptsdf)
      
      # number of particles in each pixel
      ptsdf$nparts=ceiling(ptsdf$chl*xres*yres*concdep/0.001/massperpart)
      totparts <- sum(ptsdf$nparts)
      
      # randomly scatter points within grid cells
      x1 <- rep(ptsdf$X, times=ptsdf$nparts)
      x2 <- x1 + runif(totparts, min=-0.5, max=0.5)*xres # random numbers from uniform distribution
      y1 <- rep(ptsdf$Y, times=ptsdf$nparts)
      y2 <- y1 + runif(totparts, min=-0.5, max=0.5)*yres # random numbers from uniform distribution
      
      #   plot(x1, y1)
      #   xx <- locator(2) # click two points on the plot to zoom in
      #   xlim <- c(min(xx$x),max(xx$x))
      #   ylim <- c(min(xx$y),max(xx$y))
      #   plot(x1,y1, xlim=xlim, ylim=ylim)
      #   points(x2, y2, col="red", pch=3)
      #   points(x1, y1)
      
      # project grid to lat lon
      ptsdf1 <- data.frame(X=x2,Y=y2)
      ptsdf1 <- projLatLon(ptsdf1, prj_new)
      
    } # end if initialize as in image2ini.pro IDL script
    
    #####################################################################
    print("write particle xyz to ini file")
    print(paste(Sys.time()))
    
    # write *.ini file with new image data
    if(exists("inifile")==FALSE){ # assign file name here only if not being called by habtracker.R
      inifile <- "erienc-ptraj.ini"
    }
    file.remove(inifile) # delete the file if it already exists
    
    ptm <- proc.time()

    totparts <- nrow(ptsdf1)
    ptsmat <- matrix(0.0, nrow=totparts, ncol=6)
    ptsmat[,1] <- seq(1,totparts)
    ptsmat[,2] <- ptsdf1$lon
    ptsmat[,3] <- ptsdf1$lat
#    ptsmat[,4] <- -ptsdf1$z1 # input z depth is positive down from surface, input sigma is negative down from surface
    ptsmat[,4] <- ptsdf1$sigma # input z depth is positive down from surface, input sigma is negative down from surface
    ptsmat[,5] <- ttrelease
    
    # # assign vertical velocity based on diameters measured by flowcam and empirical relation from 
    # # Nakamura 1993
    # df <- read.csv("/Users/rowe/Documents/HABproject/ToxicToledoFlowcam04082014/WE12at20%dilution4Aug2015.csv",skip=3)
    df <- read.csv(paste0(habtrackerDir, "WE12at20%dilution4Aug2015.csv"),skip=3)
    dia <- sample(df$Diameter..ESD., size=totparts, replace = TRUE)
    #### buoyant velocity and colony diameter from Nakamura et al. 1993 Wat. Res. Vol. 27, No. 6, pp. 979-983, 1993 Fig. 3
    # # from DigitizeNakamuraFig.R :
    # wb1 <- 10^(log10(dia)*1.471-1.384) # 88*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 3-Aug-1990 sample
    # wb2 <- 10^(log10(dia)*1.1694-0.4058) # 36*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 18-Sep-1990 sample
    # revised regression 9-30-2015
     wb1 <- 10^(log10(dia)*1.4439 -0.9194) # 88*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 3-Aug-1990 sample
    #wb2 <- 10^(log10(dia)*1.471  -1.380 ) # 36*(df$Diameter..ESD.)^1.5 # Nakamura Fig3 18-Sep-1990 sample
    
    ptsmat[,6] <- wb1*1e-6  # 100.0e-6 #
    
    ptsmat1 <- sprintf(fmt="%12.0f %12.6f %12.6f %10.6f %10.6f %1.6E", ptsmat[,1], ptsmat[,2] , ptsmat[,3], ptsmat[,4], ptsmat[,5], ptsmat[,6] )
    #    ptsmat1 <- sprintf(fmt="%12.0f %12.6f %12.6f %10.2f %10.2f", ptsmat[,1], ptsmat[,2] , ptsmat[,3], ptsmat[,4], ptsmat[,5] )
    
    write.table(totparts, file = inifile, append = FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(ptsmat1, file = inifile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
    # write a file containing estimated mixed layer depth
    smlfile <- paste0("sml_thick_",yyyy,ddd,"00",mm,"00.csv") 
    write.table("nlev1 is the estimated mixing depth in meters of buoyant Microcystis colonies, see Rowe et al. 2016 J. Geophys. Res. doi:10.1002/2016JC011720", file = smlfile, append = FALSE, col.names=FALSE, row.names=FALSE,  quote=FALSE)
    write.table(stas1d, file = smlfile, append = TRUE, col.names=TRUE, row.names=FALSE,  quote=FALSE, sep=",")
    
    print(paste("elapsed write time =",proc.time()[3]-ptm[3]))
    
    print(paste("totparts = ",totparts))
    print(paste("elapsed overall time =",proc.time()[3]-ptm1[3]))
    print(paste(Sys.time()))
#  } # end if ngood > 1
  
#} # end fi loop

# for some reason the first three characters of the message will get cut off
system("echo '   image2ini check' >> image2ini.log")
