#! /bin/bash

indir=`pwd`/../Input/Scans
outdir=`pwd`/../Output/Test1

# Make out dir
mkdir -pv $outdir

# Create list and configuration file
cmd="prep_dataset -i $indir --ulay _T1_DS.nii.gz --olay _T1_ICVMASK.nii.gz --olay2 _T1_BRAINMASK.nii.gz -d $outdir --mask _T1_ICVMASK.nii.gz"
$cmd

## Run MRISnaphot report
cmd="create_report -d $outdir"
$cmd

