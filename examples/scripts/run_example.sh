#! /bin/bash
indir='../scans'
outdir='../output'

## Run the command for data preparation.
cmd="mrisnapshot_prep_data -i $indir -s _T1.nii.gz -d $outdir --mask _T1_ICVMASK.nii.gz --olay _T1_ICVMASK.nii.gz --olay2 _T1_BRAINMASK.nii.gz"
echo "Running: $cmd"
$cmd

## Run the command for report creation.
cmd="mrisnapshot_create_report -d $outdir"
echo "Running: $cmd"
$cmd

