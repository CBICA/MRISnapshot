#! /bin/bash

indir=`pwd`/../Input/Scans
indir='/home/guraylab/AIBIL/Projects/MRISnapshot/Pipelines/20230215_PrepTestData/Test'
outdir=`pwd`/../Output/Test5


# Make out dir
mkdir -pv $outdir

# Create list and configuration file
cmd="prep_dataset -i $indir --ulay _T1.nii.gz --olay _T1_ROIMASK.nii.gz --mask _T1_BRAINMASK.nii.gz -d $outdir"
# echo "About to run: $cmd"
$cmd

## Run MRISnaphot report
cmd="create_report -d $outdir"
# echo "About to run: $cmd"
$cmd

