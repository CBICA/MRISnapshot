#! /bin/bash

indir=`pwd`/../Input/Scans
indir='/home/guraylab/AIBIL/Projects/MRISnapshot/Pipelines/20230215_PrepTestData/Test'
indir='/home/guray/Github/Test/MRISnapshotTestData/Test'
outdir=`pwd`/../Output/Test3

# Make out dir
mkdir -pv $outdir

# Create list and configuration file
# cmd="prep_dataset -i $indir --ulay _FL.nii.gz --olay _FL_WMLMASK.nii.gz --mask _FL.nii.gz -d $outdir"

cmd="prep_dataset -i $indir --ulay _FL.nii.gz --olay _FL_WMLMASK.nii.gz -d $outdir"

echo "About to run: $cmd"
$cmd

## Run MRISnaphot report
cmd="create_report -d $outdir"
echo "About to run: $cmd"
$cmd

