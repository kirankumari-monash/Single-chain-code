#!/bin/bash
#PBS -P g16
#PBS -q normal 
#PBS -l walltime=10:00:00
#PBS -l mem=650MB
#PBS -l ncpus=1
#PBS -l wd
#PBS -N h25_modified
# directory containing actual executable, and its name
progdir=`pwd`
progname=sens
exe=$progdir/$progname
# mrun="prun -n $PBS_NCPUS"
# mrun=prun

echo Nodes assigned
cat -b $PBS_NODEFILE

 
date=`date +"%H:%M:%S on %d %b %Y"`
echo
echo "============================================================="
echo "Timing: Commenced at $date "
# execute the program
# $mrun $exe
$exe
sleep 5
date=`date +"%H:%M:%S on %d %b %Y"`
echo "Timing: Finished at $date "
echo "============================================================="

# remove empty error files
efiles=`find . -empty`
if [ "$efiles" != "" ]; then /bin/rm  $efiles; fi

