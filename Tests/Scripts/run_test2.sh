#! /bin/bash

indir=`pwd`/../Input/Scans
outdir=`pwd`/../Output/Test2

# Make out dir
mkdir -pv $outdir

# Create list and configuration file
cmd="prep_dataset -i $indir --ulay _FL_DS.nii.gz --olay _FL_WMLMASK.nii.gz -d $outdir"
echo "About to run: $cmd"
$cmd

## Run MRISnaphot report
cmd="create_report -d $outdir"
echo "About to run: $cmd"
$cmd

