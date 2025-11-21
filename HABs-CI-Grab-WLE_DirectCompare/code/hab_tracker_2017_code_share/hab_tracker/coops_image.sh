#!/bin/sh

# Working directory
workdir='/mnt/projects/hpc/hab/hab_tracker_2017/hab_tracker/'
cd $workdir

# coops directory
coopsdir=$workdir'coops_images/tmp/'
cd $coopsdir

# modis directory
modisdir=$workdir'coops_images/modis/'

# Retrieve coops images
wget -nc -nd -nH -np -r --quiet -R 'index*' https://tidesandcurrents.noaa.gov/hab/lake_erie/data/
rm robots.txt
# MDR 8-23-2017   rm sentinel*

# Most recent file
frecent=$(ls -t *.tif | head -n1)

echo $frecent

cd ..

# Test for match
if [ -f $frecent ]; then
	echo "File exists - No new coops images"
else
	echo "Acquire new coops image..."
	cp $coopsdir$frecent .
	sensor="${frecent%%.*}"
	b="${frecent#*.}"
	yyyyddd="${b%%.*}"
	b1="${b#*.}"
	hhmm="${b1%%.*}"
	
	# Move file to image area for next HAB Tracker run
	cp $frecent $workdir'image/'$sensor'.'$yyyyddd'.'$hhmm'.L3.GL1.tif'
	
	# If lastrun is after coops image, edit lastrun.txt
	cd $workdir
	lastrun=$(<lastrun.txt)
	lastyyyyddd=${lastrun:59:7}
	echo $lastyyyyddd
	
	if [ "$lastyyyyddd" -gt "$yyyyddd" ]; then
		echo 'Need to edit lastrun.txt'
		
		# Find most recent archived run (before coops image time)
		cd $workdir'archive/'
		num1=1
		yyyydddm1=$(($yyyyddd-$num1))
		arch=$(find -maxdepth 1 -type d -name $yyyydddm1'*' | head -n1)
		arch=${arch:2:9}
		echo $arch
		
		# Delete archived folders beyond the coops image time
		for i in `seq 1 10`;
		do
			yyyydddnew=$(($yyyyddd+$i))
			fnew=$(find -maxdepth 1 -type d -name $yyyydddnew'*' | head -n1)
			flength=${#fnew}
			if [ $flength -gt 0 ]; then
				rm -R $fnew
			fi
		done
		
		# Edit the lastrun.txt
		cd $workdir
		echo '/mnt/projects/hpc/hab/hab_tracker_2017/hab_tracker/archive/'$arch'/' > lastrun.txt
	fi
	
# MDR 8-3-2017 replace with get_modis.R
#	# Get appropriate MODIS image
#	cd $modisdir
#	wget -nc -nd -nH -np -r --quiet --no-check-certificate -R '*.tif','index*' https://coastwatch.glerl.noaa.gov/modis/buf_img/
#	sensori=${sensor:0:1}
#	yyddd=${yyyyddd:2:5}
#	mfile=$(find -maxdepth 1 -type f -name $sensori'1.'$yyddd'*.LakeErie*.jpg' | head -n1)
#	echo $mfile
#	cp $mfile $workdir'webfiles/modis.jpg'
	
fi



