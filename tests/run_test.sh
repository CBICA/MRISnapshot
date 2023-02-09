#! /bin/bash

list=`pwd`/Input/Test1/listScan.csv
config=`pwd`/Input/Test1/mriqc_config_ROI.csv
outdir=`pwd`/Output/Test1

## Make out dir
mkdir -pv $outdir

## Run test
cmd="create_report -l $list -c $config -d $outdir"
echo "About to run: $cmd"
$cmd

