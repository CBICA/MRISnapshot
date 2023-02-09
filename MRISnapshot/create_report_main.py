#!/usr/bin/env python

### Import modules
import getopt
import nibabel as nib
import numpy as np
import scipy as sp
import pandas as pd
import pickle
from PIL import Image, ImageFilter
import os, sys, time, shutil
from scipy.ndimage.interpolation import zoom
import MRISnapshot.utils.img_overlays as imolay
import MRISnapshot.utils.html_utils as html

def create_report(param_filelist, param_config, param_outdir):

    ### Parse input

    ### Read input file list
    try:
        dfFiles = pd.read_csv(param_filelist)
    except:
        sys.exit("Could not read list file (" +  param_filelist + "). Aborting operations !!!");
    fileTypes = dfFiles.columns.values
    labelID = fileTypes[0]

    ### Read config file
    try:
        dfConf = pd.read_csv(param_config, comment='#', index_col='PARAM_NAME').fillna('')
    except:
        sys.exit("Could not read config file (" +  param_config + "). Aborting operations !!!");
    params = dfConf[dfConf.columns[0]].T    
    dfTmp = pd.to_numeric(params, errors='coerce')
    params = dfTmp.combine_first(params)
