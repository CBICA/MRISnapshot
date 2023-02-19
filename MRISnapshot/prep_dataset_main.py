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

    ## Create a list with underlay and overlay images
    out_list = os.path.join(params.out_dir, 'list_images.csv')
    if os.path.exists(out_list):
        logger.warning('  Output image list exists. To overwrite it, delete the image list and rerun: ' + out_list)
    else:
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
        
        ## Create output list
        df.to_csv(out_list, index=False)
        logger.info('  Created output list: ' + out_list)
    
    ## Create a default configuration file
    out_config = os.path.join(params.out_dir, 'config.csv')
    if os.path.exists(out_config):
        logger.warning('  Output config file exists. To overwrite it, delete the config file and rerun: ' + out_config)
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
                        'is_edge' : 1, 
                        'alpha_olay' : 1, 
                        'perc_high' : 100, 'perc_low' : 0, 
                        'is_out_single' : 0, 'is_out_noqc' : 0, 
                        'img_width' : 300
                    }
        
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
    
    
    
    

    
    
    
