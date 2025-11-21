# Mark Rowe 6-5-2015
# this R script will write the run parameter file for ptraj
# it is meant to be called from within habtracker.R, in which the variables will be defined

# !=========================================================================================!
#   ! INPUT FILE FOR PARAMETERS CONTROLLING EXECUTION OF OFF-LINE PARTICLE TRACKING           !
#   ! DESCRIPTION OF VARIABLES AND SUGGESTED PARAMETERS CAN BE FOUND AT BOTTOM                !
#   !                                                                                         !
#   !        FORMAT:    	                                                  !
#   !       1.) VARIABLE  = VALUE  (EQUAL SIGN MUST BE USED)                                  !
#   !       2.) FLOATING POINT VARIABLES MUST CONTAIN A PERIOD "." EX: 1.3, 2.,etc            !
#   !       3.) BLANK LINES ARE IGNORED AS ARE LINES BEGINNING WITH ! (F90 COMMENT)           !
#   !       4.) COMMENTS CAN FOLLOW VALUES IF MARKED BY !                                     !
#   !       5.) ORDER OF VARIABLES IS NOT IMPORTANT                                           !
#   !       6.) FOR MULTIPLE VALUE VARIABLES FIRST ENTRY IS NUMBER OF VARIABLES               !
#   !           TO FOLLOW (OR 0 IF NONE)                                                      !
#   !       7.) DO NOT USE COMMAS TO SEPARATE VARIABLES                                       !
#   !       8.) DO NOT EXCEED EIGHTY CHARACTERS PER LINE                                      !
#   !       9.) FOR LINE CONTINUATION ADD \\ TO END OF LINE TO FORCE CONTINUE                 !
#   !           TO NEXT LINE.  MAXIMUM 4 CONTINUATIONS                                        !
#   !       10.) TRUE = T, FALSE = F                                                          !
#   !                                                                                         !    
#   !=========================================================================================!

runpars <- c(
  paste("DTI      =",sprintf("%1.6E",dti)) #! Time step for scheme resolution (second). mod(INSTP/DTI) should be 0
  ,paste("INSTP    =",instp)  #! INPUT time step (second). 
  ,paste("DTOUT    =",dtout)   #  ! output time step in seconds (an integer multiple of DTI)
  ,paste("TDRIFT   =",tdrift)   #! Drift duration (hour)
  ,paste("OFFSET   =",offset)     #! Offset into flow file (hours) - skip before starting drift calculations
  ,paste("BACKWARD = F")    #!Flag to do backward trajectories, set OFFSET to starting record in ncdf file
  ,paste("STICKY =",sticky)      #! Flag to indicate particles stick to shores: 1=stick, other=no stick
  #         
  #         !=========Parameters Controlling Input/Output LOCATION=====================================
  #           
  ,paste("FNCUR   = ",fncur)
  ,paste0("FNINI   = ",fnini)
  ,paste("FNOUT   =",fnout)
  ,paste("FNEND   =",fnend)
  #         
  #         !=========Parameters Controlling SIGMA/CARTESIAN===========================================
  #           
  ,paste("F_DEPTH  =",f_depth)              #! Cartesian depth is kept constant if T.
  ,paste("P_SIGMA  =",p_sigma)            # ! Input particle depth in sigma if T. 
  ,paste("OUT_SIGMA  =",out_sigma)            #! Output particle depth in sigma if T.
  #         
  #         !=========Parameters Controlling Random Walk===============================================
  #           
  ,paste("IRW   =",irw  ) #! IRW (0-w/o rw; 1-hor rw; 2-vert rw; 3:hor+vert rw)
  ,paste("DHOR  =",sprintf("%1.6E",dhor)) #! DHOR (horizontal diff coeff m^2/s, miller 10^6m^2/day~11.57m^2/s)
  ,paste("DTRW  =", sprintf("%1.6E",dtrw)) #! time step (s) for RW (from the visser's criterion). mod(DTI/DTRW) should be 0 
  ,paste("WB1    = ",sprintf("%1.6E", wb))   #! buoyant velocity of particles, m/s, positive upward
  ,paste("RWMETHOD =",rwmethod) #!Select vert. random walk method: 1=Visser, 2=Milstein, 3=Grawe Milstein
  #! 4=Grawe Strong Second Order, 5=Grawe Platen Strong Second Order
  ,paste("SMOOTH =",smooth)  #!Select smoothing: 1=original 3pt, 2= 3-pt smoothing mirrored at surface and bottom,
  #!3 = no smooth
  ,paste("DOGRAWE =",dograwe) #! do a 1-D diffusion residence time test case
  ,paste("DO1D =",do1d) #! do a 1-D vertical diffusion simulation
  # ,paste0("KFILE=",kfilepath,KFILE)
  #                  
  #          !=========Parameters Controlling Concentration Output======================================
  #          
  ,paste("CONCOUT =",concout)  #! output concentration if T
  ,paste("CONCNZLEV =",concnzlev) #! number of z-levels for concentration output, INTEGER
  ,paste("CONCZRES =",sprintf("%1.6E",conczres)) #! vertical resolution for concentration output, meters, REAL
)

# write values to a file for input into ptraj
print(paste0("writing run parameter file ",getwd(),"/",runfile))
write.table("! run parameter file for ptraj ", file = runfile, append = FALSE, row.names=FALSE, col.names=FALSE, quote=FALSE)
for(ro in 1:length(runpars)){
  write.table(runpars[ro], file = runfile, append = TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)        
}