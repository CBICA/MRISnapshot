#!/usr/bin/env python
import numpy as np
import pandas as pd
import glob
import os

## Set logger  ## FIXME to be updated
import logging
format='%(levelname)-8s [%(filename)s : %(lineno)d - %(funcName)20s()] %(message)s'
format='%(levelname)-8s %(message)s'
logging.basicConfig(level=logging.DEBUG, format = '\n' + format, datefmt='%Y-%m-%d:%H:%M:%S')
logger = logging.getLogger(__name__)

##logger.setLevel(logging.DEBUG)      ## While debugging
logger.setLevel(logging.INFO)    ## FIXME Debug comments will be removed in release version
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"

def add_img_names(df, in_dir, suffix, col_name):
    if suffix != None:
        list_img = glob.glob(in_dir + os.sep + '**' + os.sep + '*' + suffix, recursive = True)
        df_new = pd.DataFrame(data = list_img, columns = [col_name])
        df_new['ScanID'] = df_new[col_name].apply(lambda x : os.path.basename(x).replace(suffix, ''))
        df_new = df_new.drop_duplicates(subset = 'ScanID')        
        df = df.merge(df_new, how = 'left', left_on = 'ScanID', right_on = 'ScanID')
    return df

def prep_dataset(params):
    '''Create image list and copy default configuration file
    '''

    ## Create a dataframe with all underlay images
    list_ulay = glob.glob(params.in_dir + os.sep + '**' + os.sep + '*' + params.s_ulay, recursive = True)
    df = pd.DataFrame(data = list_ulay, columns = ['UnderlayImg'])
    df['ScanID'] = df.UnderlayImg.apply(lambda x : os.path.basename(x).replace(params.s_ulay, ''))
    df = df.drop_duplicates(subset = 'ScanID')
    df = df[['ScanID', 'UnderlayImg']]

    ## Find and add mask and overlay images
    df = add_img_names(df, params.in_dir, params.s_mask, 'MaskImg')
    df = add_img_names(df, params.in_dir, params.s_olay, 'OverlayImg')
    df = add_img_names(df, params.in_dir, params.s_olay2, 'OverlayImg2')

    ## Create output folder
    if os.path.exists(params.out_dir) == False:
        os.makedirs(params.out_dir)
        
    ## Write output list
    out_list = os.path.join(params.out_dir, 'list_images.csv')
    if os.path.exists(out_list) == False:
        df.to_csv(out_list, index=False)
    else:
        logger.warning('Output list file exists!')
    
    ## Create default configuration file
    dict_default = {'id_col' : 'ScanID', 'ulay_col' : 'UnderlayImg', 
                    'mask_col' : 'MaskImg', 'olay_col' : 'OverlayImg', 
                    'olay_col2' : 'OverlayImg2', 'sel_vals_olay' : '', 
                    'sel_vals_olay2' : '', 'num_slice' : 5, 'view_plane' : 'A+S+C', 
                    'is_edge' : '', 'bin_olay' : 0, 'min_vox' : 0, 
                    'is_edge' : 0, 'transp' : 1, 'perc_high' : 100, 
                    'perc_low' : 0, 'is_out_single' : 0, 'is_out_noqc' : 0, 
                    'img_width' : 100
                   }
    df_config = pd.DataFrame.from_dict(dict_default, orient='index', columns=['ParamValue']).reset_index()
    df_config.columns = ['ParamName', 'ParamValue']
    
    ## Write configuration
    out_config = os.path.join(params.out_dir, 'config.csv')
    if os.path.exists(out_config    ) == False:
        df_config.to_csv(out_config, index=False)
    else:
        logger.warning('Output configuration file exists!')
    
    
    

    
    
    
