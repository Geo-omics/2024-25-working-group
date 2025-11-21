# MDR, 6-4-2015, modified from Eric Anderson's original habtracker.sh
# and hab_ptraj.sh
# re-written in R rather than bash shell
# use ptraj in place of p3d  and to run hindcast simulations
# MDR, 6-27-2016, merge latest version with changes that Eric made for 2015
# To do a coldstart on a specified date with no image, edit the timestamp in lastrun.txt
# for example, [rowe@bear hab_tracker]$ echo /mnt/hpc/hab/hab_tracker_2016/hab_tracker/archive/201617800/ > lastrun.txt
# and set coldstartOnly=TRUE. To do a hotstart with no image, there should be no image file
# in imagedir, set coldstartOnly=FALSE, and it will start a run one day after lastrun.txt, there
# also needs to be nowcast and forecast netcdf files for the appropriate date.
# MDR, 6-12-2017 update for 2017

# to run, use the following command on wings:
# [rowe@wings hab_tracker]$ pwd
# /zebra2/users/rowe/hab/hab_tracker
# [rowe@wings hab_tracker]$ R --no-save < habtracker.R > habtracker.log 2>&1
# [anderson@tiger hab_tracker]$ R --no-save < habtracker.R >& habtracker.log

#==============================================================================
# Use this script to run the nowcast/forecast bloom simulation:
#==============================================================================
# (1) Check for new satellite image
#  	- if present, then reinitialize model and run nowcast/forecast
#		- if not available, exit
# (2) Acquire GLCFS currents ( FVCOM)
#		- determine start time (based on image or end of last nowcast)
#		- end time = last hour of most recent nowcast
#		- concatenate the necessary files

# (3) Execute ptraj-nowcast
# (4) Copy nowcast end file to forecast.ini
# (5) Execute ptraj-forecast
# (6) Create plots/animations
#==============================================================================
def.par <- par(no.readonly = TRUE) # save default, for resetting...

ptm0 <- proc.time() # track run time
# set working directory
zebraDir_hpc <- "/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/"
zebraDir <- zebraDir_hpc
dontQuit <- FALSE
if(file.exists(zebraDir)==FALSE){
  zebraDir <- "/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/"
  dontQuit <- TRUE # when debugging on the mac, don't quit the R session for missing files
}
 habtrackerDir <- paste0(zebraDir, "/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/")
 habtrackerDir_hpc <- paste0(zebraDir_hpc, "/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/")
setwd("/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/")

# Currents directory (GLCFS - FVCOM)
# fvcomdirnc='/mnt/hpc/zebra/users/rowe/hab/fvcom_netcdf/' # nowcast
# fvcomdirfc='/mnt/hpc/zebra/users/rowe/hab/fvcom_netcdf/' # forecast
fvcomdirnc='/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/' # nowcast
fvcomdirfc='/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/' # forecast

# ptraj executable
exefile=paste0(zebraDir,"/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/") 

# file names for ptraj output
fntrknc <- 'erienc-ptraj.nc'
fntrkfc <- 'eriefc-ptraj.nc'

# image directory (*.txt processed file)
imagedir='./image/'
if(file.exists(imagedir) == FALSE){dir.create(imagedir)}

# temporary currents working directory
#tmpcursdir <- paste0(habtrackerDir,"tmp_curs/")
tmpcursdir <- paste0("/Users/pdenuyl/Documents/R/extract_CI/hab_tracker_2017_code_share/hab_tracker/")
if(file.exists(tmpcursdir) == FALSE){dir.create(tmpcursdir)}

# archive path 
archivePath <- paste0(habtrackerDir,"archive/")
if(file.exists(archivePath) == FALSE){dir.create(archivePath)}

# path containing output from the last run to be used in hotstart initialization
if(file.exists("lastrun.txt")){
  posPath <- scan("lastrun.txt", what=character(), nlines=1)
}else{
  print("Could not find lastrun.txt")
  posPath <- NA
}

# output path for plots
plotPath <- paste("./plots/", sep="")
if(file.exists(plotPath) == FALSE){dir.create(plotPath)}

# switch to use p3d in place of ptraj
# 6-17-2015 MDR: doP3D and doIDL probably won't function without some re-work
doP3D <- FALSE
doIDL <- FALSE   # initialize as in the image2ini.pro IDL script
concdep <- 0.1 # assume depth in m of grid cell (only used for doIDL)

# specify parameters to be used both in pre-processing and in post processing
massperpart <- 1e+10 # 1e+8   # microgram of chl mass per particle #1e+10 paper
minChl <- 12 # don't initialize any particles where chl < minChl
doSkill <- FALSE # calculate and output skill stats for hindcast simulation
ini2D <- FALSE    # apply particles to surface z layer only (conczres)
buffer_rad <- 1.5 # Buffer distance, km. Don't use satellite pixels within the specified buffer distance from land (applied only outside of SW basin)
buffer_rad2 <- 1.0 # Buffer distance, km.  (applied only INSIDE of SW basin)
doIDW <- TRUE # use nearest-neighbor interpolation to fill in missing data areas
IDWdist <-  2.0 # km, only fill in values by IDW if good obs are closer than specified distance
easternLimit <- -82.521 # don't initialize particles east of this longitude
coldstartOnly <- FALSE # if FALSE, previous run model info will be used to initialize no data areas if the new run is within forecast period of the old run
icOnly <- FALSE # if TRUE, exit habtracker after producing initial condition files and plots (don't run ptraj)
doWeb <- TRUE # copy files to webfiles directory and rename files for web display
machineName <-  "bear" #"bear" # "rhino" # will specify "ops" partition on bear, not on rhino
do1Druns <- TRUE # TRUE # Do 1D vertical mixing runs at stations
skipCur <- FALSE # if TRUE, skip concatenating netCDF files (for debugging)
doArchive <- TRUE # archive run (also used to carry model info forward to subsequent runs)

# output options
# limits for map plots
# xlimLL <- c(-83.47336, -82.01775) # zoom in on WE
# ylimLL <- c(41.36373, 42.08565)   # zoom in on WE
xlimLL <- c(-83.6, -80.9) # furthest Eastern extent of bloom on 2011-10-09
ylimLL <- c(41.3, 42.7)   # furthest Eastern extent of bloom on 2011-10-09
#source("habtracker_functions.R")
#plotPolys(shorelineLL, xlim=xlimLL, ylim=ylimLL)

doAnimation <- FALSE # assemble still images into a gif animation
plot2D <- FALSE  # calculate concentration for the map plots assuming all particles are in the surface layer

# run the script to initialize variables for the plotting script
# without running ptraj or other time consuming things
if(exists("dryrun")==FALSE){
  dryrun <- FALSE
}

#===================================
#  Set ptraj run parameters to be used
#  in both nowcast and forecast
#===================================

npdirs <- 10 # 10 # number of parallel runs #10 for paper
jobfile <- "job-1"  # file name for shell script to submit through qsub (habtracker.R will write this file)
inifile <- "erienc-ptraj.ini"
fclength <- 120 # 120 # forecast duration, hours

dti=600.0    # time step, s, REAL
instp=3600 # record interval in hydrodynamic file, INTEGER
dtout=3600    # output time step, s, INTEGER
sticky=0     # Flag to indicate particles stick to shores, INTEGER
f_depth="F"  #! Cartesian depth is kept constant if T.
out_sigma="F" #! Output particle depth in sigma if T.
backward="F" # run backward trajectory
dograwe="F"  # run Grawe test case
dhor=1.0     #! DHOR (horizontal diff coeff m^2/s, Chapra Fig 8.11: ~1m2/s for 1km grid), REAL
dtrw=-1.0 #1.0    #! time step (s) for RW (from the visser's criterion). mod(DTI/DTRW) should be 0 , REAL
wb=-1.0  #! buoyant velocity of particles, m/s, positive upward
rwmethod=3 #3 #1  # select random walk scheme
smooth=2    # select vertical diffusivity smoothing scheme
concout="T" #! output concentration if T
concnzlev=20 #! number of z-levels for concentration output, INTEGER
conczres=1.0 #! vertical resolution for concentration output, meters, REAL

system(paste("echo ",npdirs," > npdirs.txt"))
# set up directories to run a series of parallel runs
pdirs <- c()
for(ro in 1:npdirs){
  pdirs <- c(pdirs, paste0("parallel",ro))
}

if(dryrun){
  print("!!!!!!!!!!!!!!!!!!! Doing a dry run: no model execution !!!!!!!!!!!!!!!!!!!")
}

# MDR 7-6-2017 move remove old files until after check for netcdf file

#===================================
#  Check for new satellite image
#===================================

# MDR 7-10-2017 remove any old dummy files from image dir
#Clear old image files
rmFiles <- Sys.glob(paste0(imagedir,"D*.txt"))
file.remove(rmFiles)

# look for a geotiff file in the imagedir
tiffFile <- Sys.glob(paste0(imagedir,"*.tif"))[1] # if there is more than one, only look at the first one
# if there is a geotiff, write a habtracker txt file
if(file.exists(tiffFile)){
  print(paste("Convert geotiff to txt ",tiffFile))
  source("geotiff2txt.R")
  print("Get MODIS true color from Coastwatch and plot with time stamp")
  source("get_modis.R") # MDR 8-3-2017
}

imageFile <- Sys.glob(paste0(imagedir,"*.txt"))[1]
# MDR 8-3-2017  modisFile <- Sys.glob(paste0(imagedir,"*.jpg"))[1]

print('Check for new satellite image...')
print(paste(" image file =",imageFile))

if(file.exists(imageFile)){
  
  print( 'New satellite image found, reinitialize model')
  
  # get year month day hour and minute from satellite image file name
  nchar1 <- nchar(imageFile)
  yyyy <- substr(imageFile, nchar1-16, nchar1-13)
  ddd <- substr(imageFile, nchar1-12, nchar1-10)
  hh <- substr(imageFile, nchar1-9, nchar1-8)
  mm <- substr(imageFile, nchar1-7, nchar1-6)
  
  string1 <- paste("echo",yyyy,ddd,hh,"> timestamp.txt")
  system(string1)
  
}else{ # no image found
  #  print( 'No new image, quitting' )
  #  # mdr, for hindcast runs, do nothing if there is no image
  #  if(!dryrun){quit(save="no")}
  
#    print("No new image, cold start with latest nowcast")
#    ncfiles <- Sys.glob(paste0(fvcomdirnc,"*.nc"))
#    ncfiles <- ncfiles[rev(order(ncfiles))[1]]
#    
#    nchar1 <- nchar(ncfiles) - 2
#    yyyy <- substr(ncfiles,nchar1-9,nchar1-6)
#    ddd <- substr(ncfiles,nchar1-5,nchar1-3)
#    hh <- substr(ncfiles,nchar1-2,nchar1-1)
  
    print( 'No new image, initialize with lastrun.txt + 1 day ' )
    # to start on a specified date, edit the timestamp in lastrun.txt, and set coldstartOnly = TRUE
    nchar1 <- nchar(posPath)
    yyyy <- substr(posPath,nchar1-9,nchar1-6)
    ddd <- substr(posPath,nchar1-5,nchar1-3)
    ddd <- paste0(sprintf("%03.0f",as.numeric(ddd)+1))
    hh <- substr(posPath,nchar1-2,nchar1-1)
    mm="00"
    dummyImage=paste0(imagedir,"D",yyyy,ddd,"00",mm,"00.txt")
    
#   MDR 6-28-2016: the dummy file does not really need to contain any data
#    imageFile <- Sys.glob(paste0(imagedir,"*.txt"))[1]
#    system(paste("scp dummy.txt ",dummyImage))
    
    string1 <- paste("echo 'no data' > ",dummyImage)
    system(string1)
    
    hh="00"
    string1 <- paste("echo",yyyy,ddd,hh,"> timestamp.txt")
    system(string1)

} # end if image file exists

# don't initialize particles east of Pt Pelee -82.5 earlier than Aug 18
if(as.numeric(ddd) < 230){easternLimit <- -82.5 }


#===================================
#  Acquire GLCFS currents (FVCOM)
#===================================

  print( 'Acquire GLCFS (FVCOM) currents...')
  
  #Nowcast currents
  fvcomnc <- paste0(fvcomdirnc,yyyy,ddd,"00.nc")
  print(paste("Nowcast currents",fvcomnc))
  
  # quit if the nowcast file does not exist
  if(file.exists(fvcomnc)==FALSE){
    print( "Nowcast NetCDF file does not exist")
    print(fvcomnc)
    if(!dryrun & !icOnly  & !dontQuit){quit(save="no")}
  }
  
  # ptraj uses the FVCOM netcdf file, so it's not necessary to regrid the currents as for p3d
  
  #Forecast currents
#  #MDR, concatenate netcdf files for five day hindcast
#  
#  fvcomfcfiles <- c()
#  fclengthd <- ceiling(fclength/24)
#  for(fday in seq(1,fclengthd)){
#    dddf=as.numeric(ddd)+fday
#    fileChar <- as.character(dddf)
#    fileChar1 <- "000"
#    fileChar <- paste(substr(fileChar1,start=1,stop=3-nchar(fileChar)),fileChar,sep="")
#    fvcomfcfiles <- c(fvcomfcfiles, paste0(fvcomdirfc,yyyy,fileChar,'00.nc'))
#  }
#  fe <- min(file.exists(fvcomfcfiles))
#  # quit if any of the fvcomfcfiles do not exist
#  if(fe==0){
#    print( paste("quitting, Forecast NetCDF file does not exist",fvcomfcfiles[!(file.exists(fvcomfcfiles))]))
#    if(!dryrun & !icOnly  & !dontQuit){quit(save="no")}
#  }
#  fvcomfcfiles <- paste(fvcomfcfiles, collapse=" ")
#  
#  fvcomfc <- paste0(fvcomdirfc, "forecast.nc")
#  if(!dryrun & !skipCur & !icOnly){
#    file.remove(fvcomfc)
#    string1 <- paste("ncrcat -F", fvcomfcfiles, fvcomfc)
#    system(string1)
#  }else{
#    print("Skipping creation of the forecast currents file")
#  }
  
  # Edit 7/6/2015 by EJA: changed to look for LEOFS single forecast file (nowcast timestamp +1 day)
  fvcomfc <- paste0(fvcomdirfc,yyyy,sprintf("%03.0f",as.numeric(ddd)+1),"00.nc")
  print(paste("Forecast currents",fvcomfc))
  
  # quit if the forecast file does not exist
  if(file.exists(fvcomfc)==FALSE){
    print( "Forecast NetCDF file does not exist")
    print(fvcomfc)
    if(!dryrun & !icOnly  & !dontQuit){quit(save="no")}
  }

  # MDR 7-6-2017 move clear old files until after check for netcdf
  if(!dryrun){
    #Clear old files
    rmFiles <- Sys.glob("job-1.*")
    file.remove(rmFiles)
    
    rmFiles <- Sys.glob("sml_thick*")
    file.remove(rmFiles)
    
    rmFiles <- Sys.glob("*.ini")
    file.remove(rmFiles)
    
    rmFiles <- Sys.glob("*.out")
    file.remove(rmFiles)
    
    rmFiles <- Sys.glob("*ptraj.nc")
    file.remove(rmFiles)
    
    rmFiles <- Sys.glob("erie*c_run.dat")
    file.remove(rmFiles)
    
    #Clear old image files
    rmFiles <- Sys.glob(paste0(plotPath,"*.png"))
    file.remove(rmFiles)
  }  
  
#===================================
#  Do 1D runs
#===================================

if(do1Druns & !icOnly){
  # update ptraj input parameters specific to 1D runs
  p_sigma="T"  #! Input particle depth in sigma if T. 
  do1d="T"    # do 1D vertical diffusion only (u,v,w set to zero)
  irw=2        #! IRW (0-w/o rw; 1-hor rw; 2-vert rw; 3:hor+vert rw), INTEGER
  
  tdrift <- fclength + 48 # ptraj run duration, hr, for habtracker 2015, 168 = fclength + 48
  offset <- 0   # offset into flow file, hr
  runfile <- "erie1d_run.dat"
  fnout="erie1d-ptraj.nc"
  fnout1d <- fnout
  fnend="erie1d-ptraj.end"
  fnrunout="ptraj1d.out"
  casename="erie1d"
  
  # call the 1D run script
 # dryrun <- TRUE
  source("habtracker_1d.R")
#  dryrun <- FALSE
} # end if do1Druns

#===================================
#  Create the ini file for nowcast
#===================================

# run the script to create the initial particle position file for input to ptraj
# careful: scripts run with "source" share the same scope for variable names as this script
if(!dryrun){source("image2ini_ptraj_ug.R")}

# exit if doing icOnly
if(icOnly & !dryrun  & !dontQuit){
  print(paste("Stopping after initialization, icOnly = ",icOnly))
  quit(save="no")
}

#===================================
#  Edit ptraj input file for nowcast
#===================================

# update input values specific to nowcast

# set variables to start nowcast run at the time of the image and run to end of the day
ihrs <- 24-as.numeric(hh) #$(expr 24 - $hh)
tdrift <- ihrs # ptraj run duration, hr
offset <- as.numeric(hh)   # offset into flow file, hr

p_sigma= "T" #"F"  #! Input particle depth in sigma if T. 
do1d="F"    # do 1D vertical diffusion only (u,v,w set to zero)
irw= 2 #1 #2        #! IRW (0-w/o rw; 1-hor rw; 2-vert rw; 3:hor+vert rw), INTEGER
dtrw=-1.0 #-1.0

# use the absolute path for input files so that parallel runs can use the same input files
fncur=fvcomnc
paste("Nowcast currents",fvcomnc)
#MDR 8-4-2017  fnini=paste0(habtrackerDir,inifile) # this name matches the name in image2ini_ptraj.R script
fnini=paste0("../",inifile) # if the path name is too long, it will be truncated in ptraj

# use relative paths for output files so that parallel runs won't use the same files
runfile <- "erienc_run.dat"
fnout="erienc-ptraj.nc"
fnend="erienc-ptraj.end"
fnendnc=fnend
fnrunout="ptrajnc.out"
casename="erienc"

# write the run file
source(paste0(habtrackerDir, "write_ptraj_runfile.R"))

#===================================
#  Run ptraj nowcast
#===================================
print( 'Running ptraj-nowcast...')
ptm <- proc.time() # track run time

if(!dryrun){
  ri=1
  for(ri in 1:length(pdirs)){
    
    pdir <- paste0(habtrackerDir, pdirs[ri],"/")
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
    
    # copy run parameter file
    file.copy(from=paste0("../",runfile), to=paste0("./",runfile), overwrite = TRUE)

    setwd(habtrackerDir)
    
  } # end ri loop to setup parallel run dirs
    
  # write the job file

    write.table("#!/bin/bash", file = jobfile, append = FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(paste("#SBATCH -D",habtrackerDir_hpc), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    if(machineName == "bear"){write.table("#SBATCH --partition=ops", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)}
    write.table("#SBATCH -N 1", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table("#SBATCH --error=job-1.err", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table("#SBATCH --output=job-1.out", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)

    ri=1
    for(ri in 1:length(pdirs)){
      
      pdir <- paste0(habtrackerDir, pdirs[ri],"/")
      if(file.exists(pdir) == FALSE){dir.create(pdir)}

      write.table(paste0("cd ",paste0(habtrackerDir_hpc, pdirs[ri],"/")), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
      write.table(paste("/cm/shared/apps/slurm/current/bin/srun -n 1 -N 1", exefile, casename," >",fnrunout,"&"), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    } # end ri loop to start parallel runs
    write.table("wait", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)

  # start ptraj runs bear
  system("/cm/shared/apps/slurm/current/bin/sbatch job-1", wait=FALSE) 

  ###########################################################################
  # check to see if runs are finished, and don't proceed until all have finished.
  finished <- rep(0, npdirs)
  nfinished <- sum(finished)
  while(nfinished < npdirs){
    Sys.sleep(10)
    for(ri in 1:npdirs){
      pdir <- paste0(habtrackerDir, pdirs[ri],"/")
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
  print(paste("elapsed nowcast time =",round(etime,2),"hrs"))
  
  #===================================
  # forecast runs
  # continue each parallel run in its respective dir
  # update input values specific to forecast
  #===================================
  ptm <- proc.time() # track run time
  
  ri=1
  for(ri in 1:length(pdirs)){
    
    pdir <- paste0(habtrackerDir, pdirs[ri],"/")
    if(file.exists(pdir) == FALSE){dir.create(pdir)}
    
    setwd(pdir)
    
    # copy error.out file so it won't be overwritten
    file.copy(from=paste0("job-1.err"), to=paste0("job-1nc.err"), overwrite = TRUE)
    
    # edit the run parameter file
    # use the absolute path for input files so that parallel runs can use the same input files
    fncur=fvcomfc
	paste("Nowcast currents",fvcomnc)
    
    # use relative paths for output files so that parallel runs won't use the same files
    fnini="eriefc-ptraj.ini"
    file.remove(fnini)
    system(paste("ln -s ",fnendnc,fnini)) # use end file from nowcast to initialize forecast
    runfile <- "eriefc_run.dat"
    fnout="eriefc-ptraj.nc"
    fnend="eriefc-ptraj.end"
    fnrunout="ptrajfc.out"
    casename="eriefc"
    
    # update input values specific to forecast
    tdrift=fclength # ptraj run duration, hr
    offset=0   # offset into flow file, hr
    
    # write the run file
    source(paste0(habtrackerDir, "write_ptraj_runfile.R"))
    
    setwd(habtrackerDir)
  } # end ri loop to setup parallel dirs
    
# write the job file
  write.table("#!/bin/bash", file = jobfile, append = FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(paste("#SBATCH -D",habtrackerDir_hpc), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  if(machineName == "bear"){write.table("#SBATCH --partition=ops", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)}
  write.table("#SBATCH -N 1", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table("#SBATCH --error=job-1.err", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table("#SBATCH --output=job-1.out", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  ri=1
  for(ri in 1:length(pdirs)){
    
    pdir <- paste0(habtrackerDir, pdirs[ri],"/")
    if(file.exists(pdir) == FALSE){dir.create(pdir)}
    
    write.table(paste0("cd ",paste0(habtrackerDir_hpc, pdirs[ri],"/")), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(paste("/cm/shared/apps/slurm/current/bin/srun -n 1 -N 1", exefile, casename," >",fnrunout,"&"), file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  } # end ri loop to start parallel runs
  write.table("wait", file = jobfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
  
# start ptraj runs bear
system("/cm/shared/apps/slurm/current/bin/sbatch job-1", wait=FALSE) 

  
  ###########################################################################
  # check to see if runs are finished, and don't proceed until all have finished.
  finished <- rep(0, npdirs)
  nfinished <- sum(finished)
  while(nfinished < npdirs){
    Sys.sleep(10)
    for(ri in 1:npdirs){
      pdir <- paste0(habtrackerDir, pdirs[ri],"/")
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
} # end if not dryrun

etime <- (proc.time()[3]-ptm[3])/3600
print(paste("elapsed forecast time =",round(etime,2),"hrs"))

#===================================
#  	Create plots/animations
#===================================
ptm <- proc.time() # track run time
print( 'Generating plots...' )
#R --no-save < habplot_ptraj.R > habplot.log 2>&1
if(!dryrun){source( "habplot_ptraj_conc4.R" )}

etime <- (proc.time()[3]-ptm[3])/3600
print(paste("elapsed plotting time =",round(etime,2),"hrs"))

#===================================
#   Copy and rename files for web display
#===================================

if(doWeb){
  
  # set up archive dir
  webDir <- paste0(habtrackerDir, "webfiles/")
  print(paste(" Archiving files for the web in ",webDir))
  if(file.exists(webDir) == FALSE){dir.create(webDir)}
  
  # remove any old files from the web dir
  rmFiles <- Sys.glob(paste(webDir,"c*.png",sep=""))
  file.remove(rmFiles)
  rmFiles <- Sys.glob(paste(webDir,"d*.png",sep=""))
  file.remove(rmFiles)
  rmFiles <- Sys.glob(paste(webDir,"1D*.png",sep=""))
  file.remove(rmFiles)
  rmFiles <- Sys.glob(paste(webDir,"ew*.gif",sep=""))
  file.remove(rmFiles)
  
  # select the still files for the animation
  stillFiles <- Sys.glob(paste0(plotPath, "*.png"))
  nchar1 <- nchar(stillFiles)
  stillFiles1 <- substr(stillFiles, nchar1-16, nchar1-4)
  stillFiles <- stillFiles[which(is.na(as.numeric(stillFiles1))==FALSE)]
  
  nchar1 <- nchar(stillFiles)
  fnouts <- paste0(webDir, "cast", substr(stillFiles, nchar1-6, nchar1-4), ".png")
  
  file.copy(from=stillFiles, to=fnouts, overwrite=TRUE)
  
  # now for the displacement files
  stillFiles <- Sys.glob(paste0(plotPath, "*displacement*.png"))
  nchar1 <- nchar(stillFiles)
  fnouts <- paste0(webDir, substr(stillFiles, nchar1-16, nchar1))
  
  file.copy(from=stillFiles, to=fnouts, overwrite=TRUE)
  
  # 1D files
  stillFiles <- Sys.glob(paste0(plotPath, "1D*.png"))
  file.copy(from=stillFiles, to=webDir, overwrite=TRUE)
  
  # 3D station plots
  stillFiles <- Sys.glob(paste0(plotPath, "depth-time*.png"))
  file.copy(from=stillFiles, to=webDir, overwrite=TRUE)
  
  habanalysisFile <- Sys.glob(paste0(plotPath,"satellite_surface_conc.png"))[1]
  if(file.exists(habanalysisFile)){
    newhabanalysisFile <- paste0(webDir,"habanalysis.png")
    file.rename(from=habanalysisFile, to=newhabanalysisFile)
  }

  # MDR 8-3-2017  
#   if(file.exists(modisFile)){
#     newmodisFile <- paste0(webDir,"modis.jpg")
#     file.rename(from=modisFile, to=newmodisFile)
#   }
  file.copy(from=paste0(plotPath,"modis.jpg"), to=webDir, overwrite=TRUE) # MDR 8-3-2017  
  
  # copy GLCFS wind and wave plots to display on HAB Tracker page
  # corresponding to the hours of the HAB Tracker animation
  timestamps <- Sys.glob(paste0(plotPath, yyyy,"*.png"))
  ncount <- 0
  vcount <- 0
  ncdir <- "/mnt/archive/lang/GLCFS_Plots_Archive/2017/e/"
  fcdir <- "/mnt/archive/lang/GLCFS_Plots_Archive/forecast/e/"
  for(i in 1:length(timestamps)){
    time1 <- strsplit(timestamps[i], "/")[[1]]
    time1 <- time1[length(time1)]
    nchar <- nchar(time1)
    if(nchar == 17 ){
      time1 <- substr(time1, 1, 9)
      
      # for waves
      ncfile <- paste0(ncdir,"ewv",time1,".gif")
      fcfile <- paste0(fcdir,"ewv",time1,".gif")
      fileChar <- as.character(vcount)
      fileChar1 <- "0000"
      fileChar1 <- paste(substr(fileChar1,start=1,stop=4-nchar(fileChar)),fileChar,sep="")
      vfile <- paste0(webDir,"ewv",fileChar1,".gif")
      if(file.exists(ncfile)){
        file.copy(from=ncfile, to=vfile, overwrite=TRUE)
        vcount <- vcount +1
      }else if(file.exists(fcfile)){
        file.copy(from=fcfile, to=vfile, overwrite=TRUE)
        vcount <- vcount +1
      }
      
      # for wind
      ncfile <- paste0(ncdir,"ewn",time1,".gif")
      fcfile <- paste0(fcdir,"ewn",time1,".gif")
      fileChar <- as.character(ncount)
      fileChar1 <- "0000"
      fileChar1 <- paste(substr(fileChar1,start=1,stop=4-nchar(fileChar)),fileChar,sep="")
      nfile <- paste0(webDir,"ewn",fileChar1,".gif")
      if(file.exists(ncfile)){
        file.copy(from=ncfile, to=nfile, overwrite=TRUE)
        ncount <- ncount +1
      }else if(file.exists(fcfile)){
        file.copy(from=fcfile, to=nfile, overwrite=TRUE)
        ncount <- ncount +1
      }
      
    } # end if nchar=17
  } # for i
  

  
} # doWeb

#===================================
#    Copy model output to archive directory
#===================================

if(doArchive){
if(file.exists(archivePath) == FALSE){dir.create(archivePath)}
archivePath1 <- paste0(archivePath,yyyy,ddd,hh,"/")
if(file.exists(archivePath1) == FALSE){dir.create(archivePath1)}

# the file "lastrun.txt" will identify the most recent run in the 
# archive for use in hotstart initialization of a subsequent run
system(paste("echo ",archivePath1," > lastrun.txt"))

if(!dryrun){
  file.copy(from="timestamp.txt", to=archivePath1, overwrite=TRUE)
  file.copy(from="npdirs.txt", to=archivePath1, overwrite=TRUE)
#  file.copy(from="habtracker.log", to=archivePath1, overwrite=TRUE)
  stillFiles <- Sys.glob("*.log")
  file.copy(from=stillFiles, to=archivePath1, overwrite=TRUE)
  stillFiles <- Sys.glob("sml_thick_*.csv")
  file.copy(from=stillFiles, to=archivePath1, overwrite=TRUE)
  stillFiles <- Sys.glob(paste0(imagedir, "*"))
  file.copy(from=stillFiles, to=archivePath1, overwrite=TRUE)
  file.remove(stillFiles) # remove the image file from the imagedir
  print("Archiving the webfiles")
  file.copy(from="webfiles", to=archivePath1, overwrite=TRUE, recursive=TRUE)
  print("Archiving the plots")
  file.copy(from=plotPath, to=archivePath1, overwrite=TRUE, recursive=TRUE)
  
  ri=1
  print("Archiving the parallel run dirs")
  for(ri in 1:npdirs){
    
    pdir <- paste0(habtrackerDir, pdirs[ri],"/")
    pdira <- paste0(archivePath1, pdirs[ri],"/")
    if(file.exists(pdira) == FALSE){dir.create(pdira)}
    
    cpFiles <- Sys.glob(paste0(pdir, "*"))
    cpFiles1 <- Sys.glob(paste0(pdir, "*.end")) #no need to copy *.end files
    cpFiles <- cpFiles[!(cpFiles %in% cpFiles1)]
    file.copy(from=cpFiles, to=pdira, overwrite=TRUE)
  }
  
  if(do1Druns){
    print("Archiving the 1D run dirs")
    file.copy(from="stations_1d.csv", to=archivePath1, overwrite=TRUE)
    ri=1
    for(ri in 1:length(pdirs1d)){
      
      pdir <- paste0(habtrackerDir, pdirs1d[ri],"/")
      pdira <- paste0(archivePath1, pdirs1d[ri],"/")
      if(file.exists(pdira) == FALSE){dir.create(pdira)}
      
      cpFiles <- Sys.glob(paste0(pdir, "*"))
      cpFiles1 <- Sys.glob(paste0(pdir, "*.end")) #no need to copy *.end files
      cpFiles <- cpFiles[!(cpFiles %in% cpFiles1)]
      file.copy(from=cpFiles, to=pdira, overwrite=TRUE)
      
    }
  } # end if do1Druns
} # end if doArchive

print("plot microcystin map")
# MDR 7-26-2017
# plot microcystin map figure
  source( "plot_microcystin_map.R" )
  stillFiles <- Sys.glob(paste0(plotPath, "microcystin*.png"))
  file.copy(from=stillFiles, to=webDir, overwrite=TRUE)
  file.copy(from=paste0(plotPath,"microcystin_stationMap.png"), to=archivePath1, overwrite=TRUE, recursive=TRUE)

  
} # end if dryrun

#print( 'chmod -R 777 *')
#system("chmod -R 777 *")

print( 'HAB tracker complete!')
etime <- (proc.time()[3]-ptm0[3])/3600
print(paste("elapsed habtracker time =",round(etime,2),"hrs"))
print( '--------------------------------------------------')

# give time for output to write to the log file before quitting
Sys.sleep(30)


