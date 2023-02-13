#! /bin/bash

indir='/cbica/projects/LookAhead/Pipelines/LookAhead_3D_2021/Protocols'
outdir=`pwd`/Output/Test_vA

mkdir $outdir

cmd="prep_dataset -i $indir -d $outdir --ulay _T1_LPS.nii.gz --olay _T1_LPS_dlicvmask.nii.gz --mask _T1_LPS_dlicvmask.nii.gz"
echo "About to run: $cmd"
$cmd

config=`pwd`/Input/Test1/mriqc_config_ROI.csv

# ## Run test
# cmd="create_report -d $outdir2"
# echo "About to run: $cmd"
# $cmd

