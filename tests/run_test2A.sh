#! /bin/bash

config=`pwd`/Input/Test1/mriqc_config_ROI.csv
# indir=`pwd`/Input
indir='/home/guraylab/AIBIL/Github/Test/TestData/MRISnapshot/Test2_Init'
outdir=`pwd`/Output/Test2A

# Make out dir
mkdir -pv $outdir

cmd="prep_dataset -i $indir --ulay _T1.nii.gz --olay _MUSE.nii.gz -d $outdir --mask _DLICV.nii.gz"
# cmd="prep_dataset -i $indir --ulay _T1.nii.gz --olay _BMASK.nii.gz -d $outdir --mask _DLICV.nii.gz"
echo "About to run: $cmd"
$cmd

## Run test
cmd="create_report -d $outdir"
echo "About to run: $cmd"
$cmd

