#!/usr/bin/env python

### Import modules
import os, sys, time, shutil
import getopt
import numpy as np
import scipy as sp
import pandas as pd
import pickle
from PIL import Image, ImageFilter
#from scipy.ndimage.interpolation import zoom
import MRISnapshot.utils.img_overlays as imolay
import MRISnapshot.utils.html_utils as html

import nibabel as nib
import nibabel.processing as nibp
from nibabel.orientations import axcodes2ornt, ornt_transform, inv_ornt_aff

import MRISnapshot.utils.mylogger as mylogger
logger = mylogger.logger

def parse_config(df_conf, list_col_names):
    '''Read config list and check params
    '''
        
    #### Check integrity of the config dataframe
    
    ## Check column names of the config file
    if 'ParamName' not in df_conf:
        sys.exit("\nERROR: Missing column in config file: ParamName" + '\n')
    if 'ParamValue' not in df_conf:
        sys.exit("\nERROR: Missing column in config file: ParamValue" + '\n')
    
    ## Set index
    df_conf = df_conf.set_index('ParamName')

    ## Check required config values
    cols_required = ['id_col', 'ulay_col']
    for tmp_col in cols_required:
        if tmp_col not in df_conf.index:
            sys.exit("\nERROR: Missing required value in config file ParamName column: " + tmp_col + '\n')
    
    ## Create dataframe with default values
    cols_default = ['id_col', 'ulay_col', 'mask_col', 'olay_col', 'olay_col2', 'sel_vals_olay', 
                    'sel_vals_olay2', 'view_plane', 'num_slice', 'step_size_slice',
                    'min_vox', 'crop_to_mask', 'crop_to_olay', 'padding_ratio', 'bin_olay', 
                    'is_edge', 'alpha_olay', 'perc_high', 'perc_low', 
                    'is_out_single', 'is_out_noqc', 'img_width',
                    'label_checkbox1', 'label_checkbox2', 'label_editbox']
    vals_default = ['', '', '', '', '', '', 
                    '', 'A', '4', '',
                    '1', '0', '0', '', '1', 
                    '1', '1', '99', '1', 
                    '1', '0', '300',
                    'PASS', 'FAIL', 'Notes']
    if 'step_size_slice' in df_conf.index:              ## If step_size_slice is set default num_slices will be empty
        if df_conf.loc['step_size_slice'].values[0] != '':
            vals_default[cols_default.index('num_slice')] = ''
    df_default = pd.DataFrame(index = cols_default, data = vals_default, columns = ['ParamValue'])

    ## Update user config values with default ones
    df_default.update(df_conf)
    
    ## Create parameters variable
    params = df_default[df_default.columns[0]].T
    
    #### Verify that image parameters are in image list
    
    if params.id_col not in list_col_names:
        sys.exit("\nERROR: ID column missing in image list " + params.id_col + '\n')
    if params.ulay_col not in list_col_names:
        sys.exit("\nERROR: Underlay column missing in image list: " + params.ulay_col + '\n')
    if (params.mask_col != '') & (params.mask_col not in list_col_names):
            sys.exit("\nERROR: Mask column missing in image list: " + params.mask_col + '\n')
    if (params.olay_col != '') & (params.olay_col not in list_col_names):
            sys.exit("\nERROR: Overlay column missing in image list: " + params.olay_col + '\n')
    if (params.olay_col2 != '') & (params.olay_col2 not in list_col_names):
            sys.exit("\nERROR: Overlay column 2 missing in image list: " + params.olay_col2 + '\n')

    #### Check parameters
    
    ### Find number of overlay and mask images
    params['num_olay'] = int(params.olay_col != '') + int(params.olay_col2 != '')
    params['num_mask'] = int(params.mask_col != '')

    ### Set additional params
    params['img_width_single'] = 1000  ### FIXME

    ### Convert args with multiple values to lists
    params.view_plane = [n for n in params.view_plane.split('+') if n != '']
    params.sel_vals_olay = [int(n) for n in params.sel_vals_olay.split('+') if n != '']
    params.sel_vals_olay2 = [int(n) for n in params.sel_vals_olay2.split('+') if n != '']

    #logger.info(params)
    #input('p')

    ### Convert numeric args from str to int or float
    for tmp_arg in ['num_slice', 'step_size_slice', 'min_vox', 'crop_to_mask', 'crop_to_olay', 'bin_olay', 
                    'is_edge', 'is_out_single', 'is_out_noqc', 'img_width']:
        if params[tmp_arg] != '':
            params[tmp_arg] = int(params[tmp_arg])
        
    for tmp_arg in ['alpha_olay', 'perc_high', 'perc_low', 'padding_ratio']:
        if params[tmp_arg] != '':
            params[tmp_arg] = float(params[tmp_arg])
    
    ### Correct inconsistent parameters
    if (params.is_edge == 1) & (params.alpha_olay != 1):
        params.alpha_olay = 1
        logger.warning('    is_edge = 1 was selected, alpha_olay is set to 1')

    return params

def create_dir(dir_name):
    '''Create output directories
    '''
    try:        
        if not os.path.exists(dir_name):
                os.makedirs(dir_name)
    except:
        sys.exit("\nERROR: Could not create out folder !!!" + '\n');

def create_log_files(outdir):
    '''Create log files
    '''
    #logFile = outdir + os.sep + 'log_' + EXEC_NAME + '_' + startTimePretty + '.stdout'
    #errFile = outdir + os.sep + 'log_' + EXEC_NAME + '_' + startTimePretty + '.stderr'
    writeLog(logFile, '''------------------------''')

def copy_edited_js(path_utils, out_dir, report_hdr_txt):
    '''Copy js scripts
    '''
    fin = os.path.join(path_utils, 'save_qcform.js')
    fout = os.path.join(out_dir, 'save_qcform.js')
    new_text = 'textHdr = "' + report_hdr_txt + '"\n'

    open(fout, 'w').write(new_text + open(fin).read())

def copy_js(path_utils, out_dir):
    '''Copy js scripts
    '''
    shutil.copy(os.path.join(path_utils, 'misc_func.js'), out_dir)
    shutil.copy(os.path.join(path_utils, 'shortcut.js'), out_dir)
    shutil.copy(os.path.join(path_utils, 'load_back.js'), out_dir)

def get_nifti(df_images, sub_index, col_name, orient = 'LPS'):
    ''' Read nifti image in image list and reorient 
    '''
    fname = ''
    nii = None
    if col_name in df_images:
        fname = df_images.loc[sub_index][col_name]    
        nii = nib.load(fname)
        orig_ornt = nib.io_orientation(nii.affine)
        targ_ornt = axcodes2ornt(orient)
        if np.all(orig_ornt == targ_ornt) == False:
            transform = ornt_transform(orig_ornt, targ_ornt)
            nii = nii.as_reoriented(transform)
    return nii, fname

def get_img_mat(nii, orient = 'LPS'):
    ''' Reorient nifti and get data matrix
    '''
    img_mat = None
    if nii != None:
        orig_ornt = nib.io_orientation(nii.affine)
        targ_ornt = axcodes2ornt(orient)
        if np.all(orig_ornt == targ_ornt) == False:
            transform = ornt_transform(orig_ornt, targ_ornt)
            nii = nii.as_reoriented(transform)
        img_mat = nii.get_fdata()
    return img_mat

def calc_sel_slices_tmp():
    '''Select slices that will be used to create snapshots
    '''
    slices = [50, 70, 90, 110, 130, 150]      ## FIXME
    return slices

def calc_sel_slices(img_ulay, img_mask, img_olay, img_olay2, params, sub_index, sub_id):
    '''Select slices that will be used to create snapshots
       Slice selection algorithm:
        - Detect eligible voxels:
            - Mask image non-zero voxels if there is a mask
            - Overlay image non-zero voxels if there is an overlay            
            - Underlay image non-zero voxels otherwise
        - Detect eligible slices
            - Those where the number of non-zero voxels in a slice is larger than user threshold
        - Select from them based on user parameters
            - num_slices (n): Select n slices with equal step sizes
            - step_size_slice (s): Select all slices with step size s
    '''
    ## Create non-zero mask
    if img_mask is not None:
        nz_mask = (img_mask > 0).astype(int)
    else:
        if img_olay is not None:
            if img_olay2 is not None:
                img_olay = img_olay + img_olay2
            nz_mask = (img_olay > 0).astype(int)
        else:
            nz_mask = (img_ulay > 0).astype(int)

    ## Detect indices of non-zero slices
    ind_nz = np.where(np.sum(nz_mask, axis = (0, 1)) >= params.min_vox)[0]
    num_nz = ind_nz.size
    
    logger.info(ind_nz)

    sl_sel = ind_nz
    if num_nz <= 1:       ## 0 or 1 slice to select
        return sl_sel
        
    ## Select slices based on params.step_size_slice
    if params.step_size_slice != '':
        hstep = float(num_nz % params.step_size_slice) /  2
        sl_sel = np.unique(np.round(np.arange(hstep, num_nz-1-hstep, params.step_size_slice)))

    ## Select slices based on params.num_slice
    if params.num_slice != '':
        hstep = (float(num_nz - 1) / params.num_slice) / 2
        sl_sel = np.unique(np.round(np.linspace(hstep, num_nz-1-hstep, params.num_slice)))
        if sl_sel.shape[0] < params.num_slice:
            sl_sel = np.unique(np.round(np.linspace(0, num_nz-1, params.num_slice)))

    logger.info(sl_sel)
    logger.info(ind_nz)
    #input()

    return ind_nz[sl_sel.astype(int)]

def scale_img_contrast(nii_img, nii_mask, perc_low, perc_high):
    '''Change contrast of the image using percentile values
    '''
    img = nii_img.get_fdata()
    mask = img > 0
    if nii_mask != None:
        mask = nii_mask.get_fdata()
        mask = (mask > 0)
    scale_low, scale_high = np.percentile(img[mask], [perc_low, perc_high])
    img[img > scale_high] = scale_high
    img = img - scale_low
    nii_out = nib.Nifti1Image(img, nii_img.affine, nii_img.header)
    
    logger.info('      Underlay image contrast adjusted to ([min, max]): [' + str(scale_low) + ', ' + str(scale_high) + 
']')
    return nii_out

def extract_snapshot(img_ulay, img_olay, img_olay2, params, curr_view, curr_slice, slice_index, 
                     sub_id, dir_snapshots_full, list_sel_slices):
    
    # Get underlay slice
    img2d_ulay = img_ulay[:,:,curr_slice].astype(float)
    
    # Scale underlay image between 0 and 1
    img2d_ulay = (img2d_ulay - img2d_ulay.min()) / (img2d_ulay.max() - img2d_ulay.min())

    ## Resize underlay slice
    #img2d_ulay = zoom(img2d_ulay, (scX,scY), order=1)
    
    ### params.num_olay = 1  ## FIXME - just for tmp tests
    
    # Create final images and save
    snapshot_name = sub_id + '_' + curr_view + '_' + str(slice_index)
    
    if params.num_olay == 0:
        pil_under = imolay.singleImage(img2d_ulay)
        pil_under.convert('RGB').save(os.path.join(dir_snapshots_full, snapshot_name + '.png'))
    
    if params.num_olay == 1:
        img2d_olay = img_olay[:,:,curr_slice].astype(float)
        pil_under, pil_fused = imolay.overlayImage(img2d_ulay, img2d_olay,
                                                   params.alpha_olay, params.is_edge)
        pil_under.convert('RGB').save(os.path.join(dir_snapshots_full,snapshot_name + '.png'))
        pil_fused.convert('RGB').save(os.path.join(dir_snapshots_full,snapshot_name + '_olay.png'))

    if params.num_olay == 2:
        img2d_olay = img_olay[:,:,curr_slice].astype(float)
        img2d_olay2 = img_olay2[:,:,curr_slice].astype(float)

        pil_under, pil_fused = imolay.overlayImageDouble(img2d_ulay, img2d_olay, img2d_olay2, 
                                                         params.alpha_olay, params.is_edge)
        pil_under.convert('RGB').save(os.path.join(dir_snapshots_full,snapshot_name + '.png'))
        pil_fused.convert('RGB').save(os.path.join(dir_snapshots_full,snapshot_name + '_olay.png'))

    # Keep image information for later creation of html files
    snapshot_caption = 'Slice: ' + curr_view + '_' + str(list_sel_slices[slice_index] + 1)
    
    return [snapshot_caption, snapshot_name]

def crop_nifti(nii_mask, nii_arr, padding_ratio  = 0.1):
    ''' Crops a set of nifti images to the bounding box of the mask
    '''
    
    ## Min size in each dimension for the cropped image
    CROP_MIN_SIZE = 30
    nii_vox_size = np.array(nii_mask.header.get_zooms())[0:3]
    crop_min_xyz = np.ceil(float(CROP_MIN_SIZE) / nii_vox_size).astype(int)
    
    ## Get mask image
    tmp_img = nii_mask.get_fdata()
    img_dim = tmp_img.shape
    
    ## Calculate init cropping boundaries in 3 orientations
    b_coors = np.zeros([3, 2]).astype(int)
    b_sizes = np.zeros([3])
    for i, tmp_axis in enumerate([(1, 2), (0, 2), (0, 1)]):
        ind_nz = np.where(np.any(tmp_img, axis = tmp_axis))[0]
        
        ## If empty mask with all values are zero; no cropping
        if ind_nz.shape[0] == 0:
            return nii_arr
        
        ## Get non-zero boundaries
        b_min = ind_nz[0]                  ## Left boundary is the first non-zero index
        if ind_nz.shape[0] > 1:           
            b_max = ind_nz[-1]           ## Right boundary is the last non-zero index
        else:
            b_max = b_min                 ## Mask has a single non-zero slice; special case

        b_coors[i, 0] = int(b_min)
        b_coors[i, 1] = int(b_max)
        b_sizes[i] = (b_max - b_min) * nii_vox_size[i]
    
    ## Calculate padded crop size (crop to size of the view with max size) 
    b_max_size = b_sizes.max()
    b_padded_size = b_max_size + b_max_size * padding_ratio * 2
    b_padded_size = int(np.max([b_padded_size, CROP_MIN_SIZE]))
    b_padded_dims = np.ceil(float(b_padded_size) / nii_vox_size).astype(int)

    ## Calculate final cropping boundaries in 3 orientations
    for i, tmp_axis in enumerate([0, 1, 2]):
        b_center = float(b_coors[i, 0] + b_coors[i, 1]) / 2 
        b_half_size = float(b_padded_dims[i]) / 2
        b_min = int(np.ceil(b_center - b_half_size))
        b_max = int(np.ceil(b_center + b_half_size))
        b_coors[i, 0] = np.max([0, b_min])           # Correct if start of crop boundary is smaller than 0
        b_coors[i, 1] = np.min([img_dim[i] - 1, b_max])  # Correct if end of crop boundary is larger than img size
        
    ## Crop and reshape all images
    out_arr = []
    for i, tmp_nii in enumerate(nii_arr):
        if tmp_nii is None:
            out_arr.append(tmp_nii)
        else:
            ## Crop images
            cropped_nii = tmp_nii.slicer[b_coors[0,0]:b_coors[0,1], b_coors[1,0]:b_coors[1,1], 
                                         b_coors[2,0]:b_coors[2,1]]
            out_arr.append(cropped_nii)

    logger.info('      Images cropped to (x, y, z): ' + str(b_coors[0]) + 
                ', ' + str(b_coors[1]) + ', ' + str(b_coors[2]))

    return out_arr


def resize_nifti(nii_arr, interp_order):
    '''Resize to 1x1x1 mm with size max_size x max_size x max_size
    '''

    nii_vox_size = np.array(nii_arr[0].header.get_zooms())[0:3]
    nii_img_dims = np.array(nii_arr[0].get_fdata().shape)[0:3]
    max_size = int(np.ceil(np.max(nii_vox_size * nii_img_dims)))

    out_arr = []
    for i, tmp_nii in enumerate(nii_arr):
        if tmp_nii is None:
            out_arr.append(tmp_nii)
        else:
            ## Resize images
            tmp_nii = nibp.conform(tmp_nii, out_shape = [max_size, max_size, max_size],
                                       order = interp_order[i], orientation = 'LPS')
            out_arr.append(tmp_nii)

    logger.info('      Images resized to: ' + str([max_size, max_size, max_size]))

    return out_arr



def sel_vals_nifti(in_nii, sel_vals):
    '''Select a set of values from the input image
    '''
    
    ## Read img data
    in_img = in_nii.get_fdata()
    
    ## Set values not included in the selected list to zero
    in_img[np.isin(in_img, sel_vals) == False] = 0
    
    ## Create updated nifti
    #out_nii = nib.Nifti1Image(in_img.astype(np.float), nii_original_scan.affine)
    out_nii = nib.Nifti1Image(in_img, in_nii.affine, in_nii.header)

    ## Return out nifti    
    return out_nii

def create_snapshots(params, df_images, dir_snapshots_full):
    '''Create snapshots and meta data about them
    '''
    
    # Dictionary with img orientation for different views
    d_orient = {'A':'PLS', 'S':'IPR', 'C':'IRP'}  
    
    num_images = df_images.shape[0]
    
    ### Extract and save snapshots with metadata
    if not os.path.isfile(dir_snapshots_full + os.sep + 'img_info_all.pickle'):

        img_info_all = [];
        for sub_index, sub_id in enumerate(df_images[params.id_col]):

            logger.info('    Reading images for subject ' + str(sub_index + 1) + ' / ' + str(num_images))
            nii_ulay, fname_ulay  = get_nifti(df_images, sub_index, params.ulay_col)
            nii_mask, fname_mask  = get_nifti(df_images, sub_index, params.mask_col)
            nii_olay, fname_olay  = get_nifti(df_images, sub_index, params.olay_col)
            nii_olay2, fname_olay2  = get_nifti(df_images, sub_index, params.olay_col2)

            ## Select values on overlay images
            if len(params.sel_vals_olay) > 0:
                nii_olay = sel_vals_nifti(nii_olay, params.sel_vals_olay)
            if len(params.sel_vals_olay2) > 0:
                nii_olay2 = sel_vals_nifti(nii_olay2, params.sel_vals_olay2)

            ## Crop input images
            if params.crop_to_mask == 1:
                [nii_ulay, nii_mask, nii_olay, nii_olay2] = crop_nifti(nii_mask, 
                                                                       [nii_ulay, nii_mask, 
                                                                       nii_olay, nii_olay2],
                                                                       params.padding_ratio)
            if params.crop_to_olay == 1:
                [nii_ulay, nii_mask, nii_olay, nii_olay2] = crop_nifti(nii_olay, 
                                                                       [nii_ulay, nii_mask, 
                                                                       nii_olay, nii_olay2],
                                                                       params.padding_ratio)

            ## Resize input images
            interp_order = [1, 0, 0, 0]
            [nii_ulay, nii_mask, nii_olay, nii_olay2] = resize_nifti([nii_ulay, nii_mask, 
                                                                      nii_olay, nii_olay2], interp_order)

            # Initialize containers to keep image info
            snapshot_name_all = []
            snapshot_caption_all = []
            list_sel_slices_all = []

            ## Scale ulay image intensities
            nii_ulay = scale_img_contrast(nii_ulay, nii_mask, params.perc_low, params.perc_high)

            ### Create snapshots for each orientation
            for view_index, curr_view in enumerate(params.view_plane):

                ## Get data in selected orientation
                img3d_ulay = get_img_mat(nii_ulay, d_orient[curr_view])
                img3d_mask = get_img_mat(nii_mask, d_orient[curr_view]) 
                img3d_olay = get_img_mat(nii_olay, d_orient[curr_view]) 
                img3d_olay2 = get_img_mat(nii_olay2, d_orient[curr_view]) 
                       
                ### Select slices to show
                list_sel_slices = calc_sel_slices(img3d_ulay, img3d_mask, img3d_olay, 
                                                  img3d_olay2, params, sub_index, sub_id)
                logger.info('      Selected slices in view ' + str(curr_view) + ' : ' + str(list_sel_slices))
                
                for slice_index, curr_slice in enumerate(list_sel_slices):
                    info_snapshot = extract_snapshot(img3d_ulay, img3d_olay, img3d_olay2, params, 
                                                     curr_view, curr_slice, slice_index, sub_id, 
                                                     dir_snapshots_full, list_sel_slices)
                    snapshot_caption_all.append(info_snapshot[0])
                    snapshot_name_all.append(info_snapshot[1])

                list_sel_slices_all.append(list_sel_slices)
                    
            ### Keep image information for later creation of html files
            snapshot_info = {'sub_index' : sub_index, 'sub_id' : sub_id, 
                              'fname_ulay' : fname_ulay, 'fname_olay' : fname_olay, 
                              'fname_olay2' : fname_olay2, 'view_plane' : params.view_plane, 
                              'list_sel_slices_all' : list_sel_slices_all, 
                              'snapshot_name_all' : snapshot_name_all, 
                              'snapshot_caption_all' : snapshot_caption_all}
                
            img_info_all.append(snapshot_info)

        pickle.dump( img_info_all, open( dir_snapshots_full + os.sep + 'img_info_all.pickle', "wb" ) )
        
    ### If slices were already extracted, use saved meta data
    else:
        img_info_all = pickle.load( open( dir_snapshots_full + os.sep + 'img_info_all.pickle', "rb" ) )

    return img_info_all


def create_html_report(params, out_dir, dir_subjects_full, dir_snapshots, dir_snapshots_full, dir_subjects, 
img_info_all, out_report):
    ###################################################################
    ### CREATE HTML REPORTS ###########################################

    ### Create html files for: 
    ###        - single snapshots with/without overlay for each selected slice, 
    ###        - all snapshots from a subject together
    ###        - all subjects together

    ## Write the stylesheet
    HTML_stylesheet = html.html_stylesheet(params.img_width);
    ofp = open(os.path.join(dir_subjects_full, 'scripts', 'pagestyle.css'),'w')
    ofp.write(HTML_stylesheet)
    ofp.close

    ## Prepare text for final report
    HTML_multi_subj = ''
    numItem = len(img_info_all)
    fout_final = os.path.join(dir_subjects_full, 'finalPage.html')

    for i,item in enumerate(img_info_all):

        ### Read meta-data about snapshot
        sub_id = item['sub_id']
        prev_sub_id = img_info_all[np.max([0,i-1])]['sub_id']
        next_sub_id = img_info_all[np.min([numItem-1,i+1])]['sub_id']

        fname_ulay = item['fname_ulay']
        fname_olay = item['fname_olay']
        fname_olay2 = item['fname_olay2']
        snapshot_name_all = item['snapshot_name_all']
        snapshot_caption_all = item['snapshot_caption_all']

        # Header for file with multiple snapshots
        HTML_multi_snapshot = html.htmlSubjectPrefix(sub_id, i+1, numItem, fname_ulay, 
                                                     fname_olay, fname_olay2)

        ## Add QC form
        HTML_multi_snapshot = HTML_multi_snapshot + '\n' + '<p style="clear: both;"><br>'        
        if params.is_out_noqc == 0:
            # Add qc form to html page for the subject
            HTML_multi_snapshot = HTML_multi_snapshot + '\n' + \
                    html.htmlQCForm(params.label_checkbox1, params.label_checkbox2, params.label_editbox)
            
        for j, snapshot_name in enumerate(snapshot_name_all):
            snapshot_caption = snapshot_caption_all[j]
            
            # Write html page for each snapshot without overlay
            tmp_img = snapshot_name + '.png'
            tmp_txt = 'Fig: ' + str(j) + snapshot_caption
            tmp_html = snapshot_name + '_olay.html'            
            HTML_single_snapshot = html.htmlSnapshot(tmp_img, tmp_txt, tmp_html, params.img_width_single) 

            ofp = open(os.path.join(dir_snapshots_full, snapshot_name + '.html'), 'w')
            ofp.write(HTML_single_snapshot)
            ofp.close

            # Write html page for each snapshot
            if params.num_olay == 0:
                
                # Append html code for all snapshots from a subject
                tmp_img = os.path.join(dir_snapshots, snapshot_name + '.png')
                tmp_html = os.path.join(dir_snapshots, snapshot_name + '.html')

                HTML_tmp = html.htmlSubjectAddImage(tmp_img, snapshot_caption, tmp_html)

            else:
                tmp_img = snapshot_name + '_olay.png'
                tmp_txt = 'Fig: ' + str(j) + snapshot_caption
                tmp_html = snapshot_name + '.html'
                HTML_single_snapshot = html.htmlSnapshot(tmp_img, tmp_txt, tmp_html, params.img_width_single)
                
                tmp_html = os.path.join(dir_snapshots_full, snapshot_name + '_olay.html')
                
                ofp = open(tmp_html, 'w')
                ofp.write(HTML_single_snapshot)
                ofp.close

                # Append html code for all snapshots from a subject
                tmp_img = os.path.join(dir_snapshots, snapshot_name + '_olay.png')
                tmp_html = os.path.join(dir_snapshots, snapshot_name + '_olay.html')
                HTML_tmp = html.htmlSubjectAddImage(tmp_img, snapshot_caption, tmp_html)
            
            HTML_multi_snapshot = HTML_multi_snapshot + '\n' + HTML_tmp

        # Prepare html page for the subject (with all snapshots)
        html_sub = sub_id + '_mriqc.html'

        if params.is_out_single == 1:
            # Prepare html page for all subjects together
            HTML_multi_subj = HTML_multi_subj + '\n' + HTML_multi_snapshot
            
        else:
            # Add navigation to html page for the subject
            html_next_sub = next_sub_id + '_mriqc.html'
            html_prev_sub = prev_sub_id + '_mriqc.html'

            HTML_multi_snapshot = html.html_Navig(html_prev_sub, html_next_sub) + \
                                  HTML_multi_snapshot + '\n' + \
                                  html.html_Navig(html_prev_sub, html_next_sub)

        # Write html page for the subject
        ofp = open(os.path.join(dir_subjects_full, html_sub), 'w')
        ofp.write(HTML_multi_snapshot)
        ofp.close()
        

    ### Write the output report
    if params.is_out_single == 1:
        # Write html page for all subjects together (final mriqc report)
        html_sub = dir_subjects_full + os.sep + 'ALLSUBJ_mriqc.html'

        ofp = open(html_sub, 'w')
        ofp.write(HTML_multi_subj)
        ofp.close()
        
        file_report = dir_subjects + os.sep + 'ALLSUBJ_mriqc.html'

    else:
        # Create link to the report for the first subject
        first_sub_id = img_info_all[0]['sub_id']
        file_report = dir_subjects + os.sep + first_sub_id + '_mriqc.html'

    TEXT_HTML_MAINPAGE = html.htmlMainPage(file_report)
    ofp = open(out_report, 'w')
    ofp.write(TEXT_HTML_MAINPAGE)
    ofp.close()


def create_report(list_file, config_file, out_dir):

    ## Get the path for utils
    path_root = os.path.abspath(os.path.dirname(__file__))
    path_templates = os.path.join(path_root, 'templates')

    ### Check output file
    out_report = os.path.join(out_dir, 'qcreport.html')
    if os.path.exists(out_report):
        logger.warning('Output report exists. Aborting. To overwrite, delete the output report and rerun: '  \
            + out_report)
        sys.exit();

    ## Read input file list
    logger.info('  Reading image list from: ' + list_file)
    try:
        df_images = pd.read_csv(list_file)
        list_col_names = df_images.columns.values
    except:
        sys.exit("\nERROR: Could not read image list file: " + list_file + '\n')
    logger.info('    Image columns: [' + ', '.join(list_col_names) + ']')

    ## Read config file
    logger.info('  Reading config file: ' + config_file)
    try:
        df_conf = pd.read_csv(config_file, dtype = 'str').fillna('')
    except:
        sys.exit("\nERROR: Could not read config file: " +  config_file + '\n');
    logger.info('     User config data: ' + '\n' + str(df_conf))

    ## Extract and check params
    logger.info('  Checking config parameters ...')    
    params = parse_config(df_conf, list_col_names)
    
    logger.info('    Parameters for QC report (missing values are updated with default ones): ' )
    logger.info(params)

    ### Create out dirs
    dir_subjects = 'subjects'
    dir_subjects_full = os.path.join(out_dir, dir_subjects)
    dir_snapshots = 'snapshots'
    dir_snapshots_full = os.path.join(dir_subjects_full, dir_snapshots)
    dir_scripts = os.path.join(dir_subjects_full, 'scripts')

    create_dir(out_dir)
    create_dir(dir_subjects_full)
    create_dir(dir_snapshots_full)
    create_dir(dir_scripts)
    
    ## Copy js files
    report_hdr_txt = 'ID,' + params.label_checkbox1 + ',' + params.label_checkbox2 + ',' + \
                             params.label_editbox
    
    copy_js(path_templates, dir_scripts)
    copy_edited_js(path_templates, dir_scripts, report_hdr_txt)
    
    ### Create snapshots
    logger.info('  Creating snapshots ...' )
    img_info_all = create_snapshots(params, df_images, dir_snapshots_full)
    
    ### Create report
    logger.info('  Creating html report ...' )
    create_html_report(params, out_dir, dir_subjects_full, dir_snapshots, 
                       dir_snapshots_full, dir_subjects, img_info_all, out_report)
    
    
    
