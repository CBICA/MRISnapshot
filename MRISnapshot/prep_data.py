#!/usr/bin/env python
import argparse
from argparse import ArgumentParser, SUPPRESS
import numpy as np
import pandas as pd
import glob
import os
import sys

import MRISnapshot.utils.mylogger as mylogger
logger = mylogger.logger

def add_img_names(df, in_dir, suffix, col_name):
    """Finds in the input directory images with a specific suffix  
    (recursively, including nested folders), and adds them to the 
    input dataframe as a new column.
    
    :param df: Input dataframe that keeps the image names
    :param in_dir: Input directory
    :param suffix: Image suffix
    :param col_name: Column name
    :return: Output dataframe
    """

    if suffix != None:
        ## Find images with given suffix
        list_img = glob.glob(in_dir + os.sep + '**' + os.sep + '*' + suffix, recursive = True)
        df_new = pd.DataFrame(data = list_img, columns = [col_name])
        num_all = df_new.shape[0]
        logger.info('  Found ' + str(num_all) + ' ' + col_name + ' images in the input folder')

        ## Read IDs and eliminate duplicates
        df_new['ScanID'] = df_new[col_name].apply(lambda x : os.path.basename(x).replace(suffix, ''))
        df_new = df_new.drop_duplicates(subset = 'ScanID')        
        num_nodupl = df_new.shape[0]
        logger.info('  Found ' + str(num_nodupl) + ' ' + col_name + ' images with unique IDs in the input folder')

        ## Merge with list of underlay scans
        df = df.merge(df_new, how = 'left', left_on = 'ScanID', right_on = 'ScanID')
        num_match = df[df[col_name].isna() == False].shape[0]
        if num_match == 0:
            logger.warning('Found ' + str(num_match) + ' ' + col_name + ' images matching ' + 
                           'the underlay image IDs in the input folder')
        else:
            logger.info('  Found ' + str(num_match) + ' ' + col_name + ' images matching ' + 
                           'the underlay image IDs in the input folder')
        
    return df

def prep_dataset(params):
    """Function to prepare input files required for QC report creation. 
    Creates two output files in the output directory: 
    
    * list_images.csv: a list with full names of all input image files,
    
    * config.csv: a list of all configuration parameters and their default values.

    :param params: Input parameters
    """
    logger = mylogger.logger

    ## Create a list with underlay and overlay images
    out_list = os.path.join(params.out_dir, 'list_images.csv')
    if os.path.exists(out_list):
        logger.warning('  File ' + out_list + ' already exists, skipping.')
    else:
        ## Find underlay images
        list_ulay = glob.glob(params.in_dir + os.sep + '**' + os.sep + '*' + params.s_ulay, recursive = True)
        df = pd.DataFrame(data = list_ulay, columns = ['UnderlayImg'])
        num_all = df.shape[0]
        logger.info('  Found ' + str(num_all) + ' underlay images in the input folder')
        
        ## Read IDs and eliminate duplicates
        df['ScanID'] = df.UnderlayImg.apply(lambda x : os.path.basename(x).replace(params.s_ulay, ''))
        df = df.drop_duplicates(subset = 'ScanID')
        df = df[['ScanID', 'UnderlayImg']]

        ### Eliminate scan if ID is empty
        if df[df.ScanID == ''].shape[0]:
            logger.warning('  Removing image with empty ScanID from list: ' + 
                           df[df.ScanID==''].UnderlayImg.values[0])
            df = df[df.ScanID != '']

        if df.shape[0] == 0:
            logger.warning('  No underlay images were found in the input folder. Output image list is empty!')
        else:
            logger.info('  Found ' + str(df.shape[0]) + ' underlay images with unique IDs in the input folder')

            ## Add mask and overlay images
            df = add_img_names(df, params.in_dir, params.s_mask, 'MaskImg')
            df = add_img_names(df, params.in_dir, params.s_olay, 'OverlayImg')
            df = add_img_names(df, params.in_dir, params.s_olay2, 'OverlayImg2')

        ## Create output folder
        if os.path.exists(params.out_dir) == False:
            os.makedirs(params.out_dir)
        
        ## Create output list
        df.to_csv(out_list, index=False)
        logger.info('  Created output list: ' + out_list)
        
    
    ## Create a default configuration file
    out_config = os.path.join(params.out_dir, 'config.csv')
    if os.path.exists(out_config):
        logger.warning('  File ' + out_config + ' already exists, skipping.')
    else:
        dict_default = {'id_col' : 'ScanID', 
                        'ulay_col' : 'UnderlayImg', 'mask_col' : 'MaskImg', 
                        'olay_col' : 'OverlayImg', 'olay_col2' : 'OverlayImg2', 
                        'sel_vals_olay' : '', 'sel_vals_olay2' : '',
                        'view_plane' : 'A+S+C',
                        'num_slice' : 5, 'step_size_slice' : '',
                        'min_vox' : 1,
                        'crop_to_mask' : 0,
                        'crop_to_olay' : 0,
                        'padding_ratio' : 0,
                        'bin_olay' : 0,
                        'segment_olay' : 0,
                        'num_classes_olay' : 0,
                        'is_edge' : 1,
                        'alpha_olay' : 1, 
                        'perc_high' : 100, 'perc_low' : 0, 
                        'is_out_single' : 0, 'is_out_noqc' : 0, 
                        'img_width' : 300,
                        'label_checkbox1' : 'PASS',
                        'label_checkbox2' : 'FAIL',
                        'label_editbox' : 'Notes'}
        
        ## Reset default values for missing image types
        if params.s_mask == None:
            dict_default['mask_col'] = ''
        if params.s_olay == None:
            dict_default['olay_col'] = ''
        if params.s_olay2 == None:
            dict_default['olay_col2'] = ''
        
        df_config = pd.DataFrame.from_dict(dict_default, orient='index', columns=['ParamValue']).reset_index()
        df_config.columns = ['ParamName', 'ParamValue']

        ## Create output config file
        df_config.to_csv(out_config, index=False)
        logger.info('  Created config file: ' + out_config)
    
def main():
    """Helper script to prepare input files required for QC report creation.

    :param indir: Input image directory (full or relative path) (required)
    :param s_ulay: Suffix of underlay images (required)
    :param outdir: Input image directory (full or relative path) (required)
    :param s_mask: Suffix of mask images (optional)
    :param s_olay: Suffix of underlay images (optional)
    :param s_olay2: Suffix of underlay images (optional)
    """

    descr = 'Helper script to prepare input files (image list and configuration file) ' \
            'required for QC report creation. After running this script, the user can ' \
            'manually edit the files created in OUTDIR, and run the command ' \
            '"snap_create_report [OUTDIR]" to create the QC report for their dataset ' \
            'and selected configuration.\n\n' \
            'See https://cbica.github.io/MRISnapshot for usage and examples.\n\n'

    ## Create parser
    parser = argparse.ArgumentParser(add_help=False,
                                     prog="mrisnap_prep_data",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description = descr,
                                     epilog = 'Contact: guray.erus@pennmedicine.upenn.edu')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Add back help 
    optional.add_argument(
        '-h',
        '--help',
        action = 'help',
        default = SUPPRESS,
        help = 'show this help message and exit'
    )
    
    required.add_argument("-i", dest="indir", type=str, required=True, help="Input image directory")
    required.add_argument("-s", dest="s_ulay", type=str, required=True, help="Suffix of underlay images")
    required.add_argument("-d", dest="outdir", type=str, required=True, help="Output directory")
    optional.add_argument("--mask", dest="s_mask", type=str, help = "Suffix of mask images")
    optional.add_argument("--olay", dest="s_olay", type=str, help = "Suffix of overlay images")
    optional.add_argument("--olay2", dest="s_olay2", type=str, help="Suffix of second overlay images")

    # Parse input params
    params = parser.parse_args()

    logger.info('-----------------------------------------')    
    logger.info('Running : ' + ' '.join(sys.argv))
    logger.info('-----------------------------------------')    

    ## Update params
    params.in_dir = os.path.abspath(params.indir)
    params.out_dir = os.path.abspath(params.outdir)
    
    ## Prepare data
    logger.info('  Preparing image list and configuration file ...')
    prep_dataset(params)
    
    
    

    
    
    
