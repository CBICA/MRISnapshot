#! /bin/bash

indir=`pwd`/../Input/Scans
outdir=`pwd`/../Output/Test1

# Make out dir
mkdir -pv $outdir

# Create list and configuration file
cmd="mrisnapshot_prep_data -i $indir -s _T1_DS.nii.gz -d $outdir --mask _T1_ICVMASK.nii.gz --olay _T1_ICVMASK.nii.gz 
--olay2 _T1_BRAINMASK.nii.gz"
$cmd

# echo $cmd
# echo lllllllll


# Run MRISnaphot report
cmd="mrisnapshot_create_report -d $outdir"
$cmd

