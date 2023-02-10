#! /bin/bash

list=`pwd`/Input/Test1/listScan.csv
config=`pwd`/Input/Test1/mriqc_config_ROI.csv
outdir=`pwd`/Output/Test1

indir=`pwd`/Input
outdir2=`pwd`/Output/Test2

## Make out dir
mkdir -pv $outdir

# ## Run test
# cmd="create_report -l $list -c $config -d $outdir"
# echo "About to run: $cmd"
# $cmd

cmd="prep_dataset -i $indir --s_ulay _T1.nii.gz --s_olay _ROI.nii.gz -d $outdir2 --s_mask _bmask.nii.gz"
echo "About to run: $cmd"
$cmd

