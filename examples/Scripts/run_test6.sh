#! /bin/bash

indir='/home/guraylab/AIBIL/Projects/MRISnapshot/Pipelines/20230220_Test-MUSE/Scripts/Yuhan_replicate/Scans'
outdir='/home/guraylab/AIBIL/Projects/MRISnapshot/Pipelines/20230220_Test-MUSE/Scripts/Yuhan_replicate/MRISnapshots'

# Make out dir
mkdir -pv $outdir

# Create list and configuration file
cmd="prep_dataset -i $indir --ulay _T1_LPS.nii.gz --olay _T1_LPS_dlicv_fastbc_muse.nii.gz --olay2 _T1_LPS_dlicv_dlmuse.nii.gz --mask _T1_LPS_dlicvmask.nii.gz -d $outdir"
# echo "About to run: $cmd"
$cmd

## Run MRISnaphot report
cmd="create_report -d $outdir"
# echo "About to run: $cmd"
$cmd

