# Mark Rowe 8-2-2017
# get MODIS image from Coast Watch and plot with time stamp

print("load jpeg package")
require(jpeg)

# modis directory
modisdir=paste0(habtrackerDir,'coops_images/modis/')

# coops images directory
coopsdir=paste0(habtrackerDir,'coops_images/tmp/')

setwd(modisdir)
system("wget -nc -nd -nH -np -r --quiet --no-check-certificate -R '*.tif','index*' https://coastwatch.glerl.noaa.gov/modis/buf_img/")
setwd(habtrackerDir)

# get MODIS image corresponding to image used to initialize HAB Tracker
#string2 <- strsplit(tiffFile, "/")[[1]][3]
#string2 <- substr(string2, 1,3)

# MDR 8-25-2017 HAB Tracker may use MODIS or Sentinal 3 image, so get most recent MODIS CI image
fns1 <- Sys.glob(paste0(coopsdir,"terra*.tif"))
fns2 <- Sys.glob(paste0(coopsdir,"aqua*.tif"))
fns1 <- c(fns1, fns2)
fndates <- c()
for(i in 1:length(fns1)){
  string1 <- strsplit(fns1[i], "/")[[1]]
  string1 <- string1[length(string1)]
  string2 <- strsplit(string1, "\\.")[[1]]
  fndates <- c(fndates, paste0(string2[2]))
}
xx <- which(fndates==max(fndates, na.rm=TRUE))[1]
fn1 <- fns1[xx]
string1 <- strsplit(fn1, "/")[[1]]
string1 <- string1[length(string1)]
string2 <- substr(string1, 1,3)
imageDate <- as.POSIXct(paste0(fndates[xx],"1700"), format="%Y%j%H%M", tz="UTC")

if(string2 == "ter"){
  sensori <- "t1"
}else{
  sensori <- "a1"
}

yyyy <- format(imageDate, "%Y")
ddd <- format(imageDate, "%j")

yyddd <- substr(paste0(yyyy,ddd), 3, 7)
fn2 <- paste0(sensori,".",yyddd,"*LakeErie*.jpg")
fnmodis <- Sys.glob(paste0(modisdir,fn2))[1]
print(paste("Modis file = ",fnmodis))

string2 <- strsplit(fnmodis, "/")[[1]]
string2 <- string2[length(string2)]
modishr <- substr(string2, 10,13)

if(is.na(fnmodis)==FALSE){
  jpg <- readJPEG(fnmodis)
  
  res = dim(jpg)[2:1] # get the resolution, [x, y]

  jpeg(filename = paste0(plotPath,"modis.jpg"),
       width = res[1], height = res[2], units = "px", pointsize = 12,
       quality = 100)
  
  par(mar=c(0,0,0,0))
  plot(1,1,xlim=c(1,res[1]),ylim=c(1,res[2]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(jpg,1,1,res[1],res[2])
  
  dy <- 100
  polygon(x=c(1,res[1],res[1],1), y=c(res[2],res[2],res[2]-dy,res[2]-dy), col="white")
  
  plotDate <- as.POSIXct(paste(yyyy, ddd, modishr), format="%Y %j %H%M")
  
  text(grconvertX(0.50, from="npc"), grconvertY(0.95, from="npc"), format(plotDate-4*3600, "%Y-%m-%d %H:%M EDT"), col="black", cex=5)
  
  graphics.off()

} # end if fnmodis == false




