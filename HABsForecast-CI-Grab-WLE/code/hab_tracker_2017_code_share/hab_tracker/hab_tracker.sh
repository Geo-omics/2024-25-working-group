#!/bin/sh

# Working directory
workdir='/mnt/projects/hpc/hab/hab_tracker_2017/hab_tracker/'
cd $workdir

## echo 'Running HAB Tracker'
#date

# lastrun=$(<lastrun.txt)
# lastddd=${lastrun:63:3}
# nowddd=$(date -u +%j)
# num1=1
# nowdddm1=$(($nowddd-$num1))
# while [ $nowdddm1 -gt $lastddd ]; do
#	echo $nowdddm1
#	echo $lastddd 
#	R --no-save < habtracker.R >& habtracker.log
#	lastrun=$(<lastrun.txt)
#	lastddd=${lastrun:63:3}
# done
#
#R --no-save < habtracker.R >& habtracker.log
R --no-save < Bottom_temp_animation_cur.R >& bottomtemp.log


