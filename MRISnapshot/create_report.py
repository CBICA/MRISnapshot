#!/usr/bin/env python

### Import modules
import argparse
from argparse import ArgumentParser, SUPPRESS

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

    :param df_conf: Dataframe with the list of configuration parameters
    :param list_col_names: List of column names in image list

    :return params: List of all input parameters
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
                    'segment_olay', 'num_classes_olay',
                    'is_edge', 'alpha_olay', 'perc_high', 'perc_low', 
                    'is_out_single', 'is_out_noqc', 'img_width',
                    'label_checkbox1', 'label_checkbox2', 'label_editbox']
    vals_default = ['', '', '', '', '', '', 
                    '', 'A', '4', '',
                    '1', '0', '0', '', '1',
                    '0', '0',
                    '1', '1', '99', '1', 
                    '1', '0', '300',
                    'PASS', 'FAIL', 'Notes']
    
    ## If step_size_slice is set, default num_slices will be empty
    if 'step_size_slice' in df_conf.index:              
        if df_conf.loc['step_size_slice'].values[0] != '':
            vals_default[cols_default.index('num_slice')] = ''
    df_default = pd.DataFrame(index = cols_default, data = vals_default, columns = ['ParamValue'])

    ## Update user config values with default ones
    df_default.update(df_conf)
    
    ## Create parameters variable
    params = df_default[df_default.columns[0]].T
    
    #### Verify that image parameters match with the contents of image list
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
    params['img_width_single'] = 1000  ### FIXME This can be a user parameter in future

    ### Convert args with multiple values to lists
    params.view_plane = [n for n in params.view_plane.split('+') if n != '']
    params.sel_vals_olay = [int(n) for n in params.sel_vals_olay.split('+') if n != '']
    params.sel_vals_olay2 = [int(n) for n in params.sel_vals_olay2.split('+') if n != '']

    ### Convert numeric args from str to int or float
    for tmp_arg in ['num_slice', 'step_size_slice', 'min_vox', 'crop_to_mask', 'crop_to_olay', 'bin_olay',
                    'segment_olay', 'num_classes_olay',
                    'is_edge', 'is_out_single', 'is_out_noqc', 'img_width']:
        if params[tmp_arg] != '':
            params[tmp_arg] = int(params[tmp_arg])
        
    for tmp_arg in ['alpha_olay', 'perc_high', 'perc_low', 'padding_ratio']:
        if params[tmp_arg] != '':
            params[tmp_arg] = float(params[tmp_arg])
    
    ### Correct inconsistent parameters
    if (params.is_edge == 1) & (params.alpha_olay != 1):
        params.alpha_olay = 1
        logger.warning('    Parameter alpha_olay is set to 1, because is_edge was set to 1')

    ### Correct inconsistent parameters
    if (params.segment_olay == 1):
        if params.num_classes_olay > 10:
            params.num_classes_olay = 10
            logger.warning('    Parameter num_classes_olay is set to 10 (max value)')
        if params.num_classes_olay < 2:
            params.num_classes_olay = 2
            logger.warning('    Parameter num_classes_olay is set to 2 (min value)')
        if params.is_edge == 0:
            params.is_edge = 1
            logger.warning('    Parameter is_edge is set to 1, because segment_olay was set to 1')
        params.alpha_olay = 1
        logger.warning('    Parameter alpha_olay is set to 1, because is_edge was set to 1')

    return params

def create_dir(dir_name):
    '''Create output directories

    :param dir_name: Output folder
    '''
    try:        
        if not os.path.exists(dir_name):
                os.makedirs(dir_name)
    except:
        sys.exit("\nERROR: Could not create out folder !!!" + '\n');

def copy_edited_js(path_utils, out_dir, report_hdr_txt):
    '''Copy js scripts to QC report folder, while adding hdr text.
    This function was required to modify QC form columns in js
    template file dynamically based on user input.
    
    :param path_utils: Path to js template files
    :param out_dir: Output folder
    :param report_hdr_txt: Header text to add to the js file
    '''
    
    try:       
        fin = os.path.join(path_utils, 'save_qcform.js')
        fout = os.path.join(out_dir, 'save_qcform.js')
        new_text = 'textHdr = "' + report_hdr_txt + '"\n'

        open(fout, 'w').write(new_text + open(fin).read())
    except:
        sys.exit("\nERROR: Could not copy edited template .js file to output directory " + out_dir + '\n')

def copy_js(path_utils, out_dir):
    '''Copy js scripts to QC report folder
    
    :param path_utils: Path to js template files
    :param out_dir: Output folder
    '''
    try:        
        shutil.copy(os.path.join(path_utils, 'misc_func.js'), out_dir)
        shutil.copy(os.path.join(path_utils, 'shortcut.js'), out_dir)
        shutil.copy(os.path.join(path_utils, 'load_back.js'), out_dir)
    except:
        sys.exit("\nERROR: Could not copy template .js files to output directory " + out_dir + '\n');

def read_and_check_images(df_images, params, sub_index, orient = 'LPS'):
    ''' Read underlay, mask and overlay image names for a subject from a dataframe, 
    read images as Nifti, reorient them, and check that dimensions are consistent
    
    :param df_images: Input dataframe with image names
    :param sub_index: Index of the current subject
    :param orient: Desired image orientation
    
    :return qc_ok_flag: Binary flag that indicates if images are OK and consistent
    :return qc_msg: Text message to report QC issues
    :return nii: Output nifti images
    :return img_mat: Image file names
    '''
    nii_out = []
    fnames_out = []
    col_names = [params.ulay_col, params.mask_col, params.olay_col, params.olay_col2]
    qc_ok_flag = 1
    qc_msg = 'PASS'
    ref_affine = []         ## Check that all images have the same affine (in the same space)
    
    ## Read and check each image 
    for i, col_name in enumerate(col_names):
        if col_name == '':          ## If the specific image column is not used in the report, skip it
            fname = ''
            nii = None
        else:
            fname = df_images.loc[sub_index][col_name]
            if fname == '':                                     ## Image missing, return with error flag
                qc_ok_flag = 0
                qc_msg = 'Missing ' + col_name + ' image'
                logger.warning('   ' + qc_msg + ', subject discarded!')
                return qc_ok_flag, qc_msg, nii_out, fnames_out
                
            else:
                try:
                    nii = nib.load(fname)
                except:
                    nii = None
                    qc_ok_flag = 0
                    qc_msg = 'Could not read ' + col_name + ' image ' + fname
                    logger.warning('   ' + qc_msg + ', subject discarded!')
                    return qc_ok_flag, qc_msg, nii_out, fnames_out

                try:
                    orig_ornt = nib.io_orientation(nii.affine)
                    targ_ornt = axcodes2ornt(orient)
                    if np.all(orig_ornt == targ_ornt) == False:
                        transform = ornt_transform(orig_ornt, targ_ornt)
                        nii = nii.as_reoriented(transform)
                    if ref_affine == []:
                        ref_affine = nii.affine
                    else:

                        ## We compare the affine matrix for underlay and overlay images
                        ##   - we use 3 digit precision
                        ##     otherwise it may be too sensitive to numeric differences in headers
                        if (np.abs(ref_affine - nii.affine).sum()) > 0.001:
                            qc_ok_flag = 0
                            qc_msg = 'Inconsistent image vs. overlay: affine matrices in headers not consistent'
                            logger.warning('   ' + qc_msg + ', subject discarded!')
                            return qc_ok_flag, qc_msg, nii_out, fnames_out
                except:
                    nii = None
                    qc_ok_flag = 0
                    qc_msg = 'Could not reorient ' + col_name + ' image ' + fname
                    logger.warning('   ' + qc_msg + ', subject discarded!')
                    return qc_ok_flag, qc_msg, nii_out, fnames_out
                    
                        
        nii_out.append(nii)
        fnames_out.append(fname)
        
    return qc_ok_flag, qc_msg, nii_out, fnames_out

def get_img_mat(nii, orient = 'LPS'):
    ''' Reorient nifti and get data matrix
    
    :param nii: Input image
    :param orient: Desired image orientation
    
    :return img_mat: Image voxels in a numpy matrix    
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

#def calc_sel_slices_tmp():
    #'''Select slices that will be used to create snapshots
    #'''
    #slices = [50, 70, 90, 110, 130, 150]      ## FIXME
    #return slices

def calc_sel_slices(img_ulay, img_mask, img_olay, img_olay2, params, sub_index, sub_id):
    '''Select slices that will be used to create snapshots
    
    :param img_ulay: Underlay image
    :param img_mask: Mask image
    :param img_olay: Overlay image
    :param img_olay2: Second overlay image
    :param params: Input parameters
    :param sub_index: Index of the subject
    :param sub_id: Id of the subject
    
    :return ind_sel: A list of indices for selected slices
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

    ind_sel = ind_nz[sl_sel.astype(int)]
    return ind_sel

def scale_img_contrast(nii_img, nii_mask, perc_low, perc_high):
    '''Change contrast of the image using percentile values. Image intensities 
    will be mapped to values calculated from the min and max percentiles.
    
    :param nii_img: Input nifti image
    :param nii_mask: Input nifti mask
    :param perc_low: Lower percentile value to determine new min intensity value
    :param perc_high: Higher percentile value to determine new max intensity value
    
    :return nii_out: Output nifti image
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

def digitize_olay(nii, num_classes, perc_low = 5, perc_high = 95):
    '''Simple quick segmentation of overlay image by quantizing intensities.
    The purpose is to make sure that edge detection on overlay will be smooth. 
    Bins for digitization will be estimated from robust intensity ranges
    (using percentile values to exclude outliers)
    
    :param nii: Input nifti image
    :param num_classes: Number of classes (digits) in output image
    
    :return nii_out: Output nifti image
    '''
    ## Check input img
    if nii == None:
        return nii
    
    ## Get img data
    tmp_img = nii.get_fdata()
    
    ## Find min and max intensities
    pmin, pmax = np.percentile(tmp_img.flatten(), [perc_low, perc_high])
    
    ## Digitize image
    tmp_img = np.digitize(tmp_img, np.arange(pmin, pmax, (pmax-pmin)/num_classes))
    
    ## Create out nifti
    nii_out = nib.Nifti1Image(tmp_img, nii.affine, nii.header)

    ## Return out nifti
    logger.info('      Overlay image segmented to ' + str(num_classes) + ' classes')
    return nii_out

def extract_snapshot(img_ulay, img_olay, img_olay2, params, curr_view, curr_slice, slice_index, 
                     sub_id, dir_snapshots_full, list_sel_slices):
    ''' Extracts an image snapshot based on input parameters, and writes it to 
    output folder. Returns the name and caption of the extracted snapshot.
    
    :param img_ulay: Underlay image
    :param img_olay: Overlay image
    :param img_olay2: Second overlay image
    :param curr_view: View plane to extract the snapshot
    :param curr_slice: Slice number to extract the snapshot
    :param slice_index: Index of the extracted slice
    :param sub_id: Id of the current subject
    :param dir_snapshots_full: Output directory for snapshots (full path)
    :param list_sel_slices: List of selected slices

    :return snapshot_caption: Caption of the extracted snapshot
    :return snapshot_name: Name of the extracted snapshot
    '''
    
    # Get underlay slice
    img2d_ulay = img_ulay[:,:,curr_slice].astype(float)
    
    # Scale underlay image between 0 and 1
    if img2d_ulay.max() - img2d_ulay.min() > 0:
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
    ''' Crops a set of nifti images to the bounding box of the mask. Images are also
    padded in each direction by the padding ratio
    
    :param nii_mask: Input nifti mask  image
    :param nii_arr: An array of input nifti images
    :param padding_ratio: Padding ratio
    :return out_arr: An array of output nifti images
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
    '''Resize image to 1x1x1 mm with size max_size x max_size x max_size
    
    :param nii_arr: An array of input nifti images
    :param interp_order: Order of interpolation
    :return out_arr: An array of output nifti images
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
    '''Selects a set of values from the input image, and resets all other values to 0
    
    :param in_nii: Input nifti image
    :param sel_vals: List of selected values
    :return out_nii: Output nifti image
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

def check_foreground_mask(params, nii_mask, nii_olay, nii_olay2):
    '''Check if mask and overlay images have non-zero foreground voxels for slice selection
    
    :param in_nii: Input nifti image
    :param sel_vals: List of selected values
    :return out_nii: Output nifti image
    '''

    qc_ok_flag = 1
    qc_msg = 'PASS'

    if params.num_mask == 1:
        num_nz_mask = (nii_mask.get_fdata() > 0).sum()
        if num_nz_mask == 0:
            qc_ok_flag = 0
            qc_msg = 'Mask image has 0 foreground voxels'
            logger.warning('   ' + qc_msg + ', subject discarded!')

    else:
        if params.num_olay == 1:
            num_nz_olay1 = (nii_olay.get_fdata() > 0).sum()
            if num_nz_olay1 == 0:
                qc_ok_flag = 0
                qc_msg = 'Overlay image has 0 foreground voxels'
                logger.warning('   ' + qc_msg + ', subject discarded!')
        if params.num_olay == 2:
            num_nz_olay1 = (nii_olay.get_fdata() > 0).sum()
            num_nz_olay2 = (nii_olay2.get_fdata() > 0).sum()
            if num_nz_olay1 + num_nz_olay2 == 0:
                qc_ok_flag = 0
                qc_msg = 'Overlay images have 0 foreground voxels in total'
                logger.warning('   ' + qc_msg + ', subject discarded!')

    return qc_ok_flag, qc_msg 

def create_snapshots(params, df_images, dir_snapshots_full, out_dir):
    '''Creates and returns image snapshots and meta-data about snapshots.
    Also creates a QC dataframe to keep QC PASS/FAIL information for initial 
    scans and saves it in the output directory
    
    :param params: Input parameters
    :param df_images: Dataframe with image names
    :param dir_snapshots_full: Output directory for snapshots (full path)
    
    :return img_info_all: A dictionary with meta-data about extracted snapshots
    '''
    
    df_qc_images = None
    
    # Dictionary with img orientation for different views
    d_orient = {'A':'PLS', 'S':'IPR', 'C':'IRP'}  
    
    num_images = df_images.shape[0]
    
    ### Snapshots were already extracted, use saved snapshots and meta-data
    fname_img_info_all = os.path.join(dir_snapshots_full, 'img_info_all.pickle')
    if os.path.isfile(fname_img_info_all):
        logger.info('  File ' + fname_img_info_all + ' already exists, ' +
                    ' skipping extraction of snapshots')
        logger.warning('  Pre-saved data may be incomplete or incorrect. Please delete ' + 
                       'the complete QCReport output folder and rerun to re-create it')
        try:
            img_info_all = pickle.load(open(fname_img_info_all, "rb"))
        except:
            sys.exit("\nERROR: Could not read pre-saved data:" + fname_img_info_all + '\n');

    ### Extract snapshots and metadata
    else:
        
        qc_ok_flag_all = []
        qc_msg_all = []

        img_info_all = [];
        for sub_index, sub_id in enumerate(df_images[params.id_col]):

            ## Read images
            logger.info('    Reading images for subject ' + str(sub_index + 1) + ' / ' + str(num_images))
            qc_ok_flag, qc_msg, nii_all, fname_all = read_and_check_images(df_images, params, 
                                                                           sub_index, orient = 'LPS')
            
            ## Save image QC info to record QC fail cases, and report them in a QC csv output file
            if qc_ok_flag == 0:
                qc_ok_flag_all.append(qc_ok_flag)
                qc_msg_all.append(qc_msg)
            
            ## Create snapshots if images passed the QC
            else:
                [nii_ulay, nii_mask, nii_olay, nii_olay2] = nii_all
                [fname_ulay, fname_mask, fname_olay, fname_olay2] = fname_all

                ## Select values on overlay images
                if len(params.sel_vals_olay) > 0:
                    nii_olay = sel_vals_nifti(nii_olay, params.sel_vals_olay)
                if len(params.sel_vals_olay2) > 0:
                    nii_olay2 = sel_vals_nifti(nii_olay2, params.sel_vals_olay2)

                ## Check that mask or overlay images have non-zero voxels to select slices
                qc_ok_flag, qc_msg = check_foreground_mask(params, nii_mask, nii_olay, nii_olay2)

                ## Save image QC info
                qc_ok_flag_all.append(qc_ok_flag)
                qc_msg_all.append(qc_msg)
                
                ## Create snapshots if images passed the QC
                if qc_ok_flag == 1:
                    
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

                    ## Digitize olay image intensities
                    if params.segment_olay == 1:
                        nii_olay = digitize_olay(nii_olay, params.num_classes_olay)
                        nii_olay2 = digitize_olay(nii_olay2, params.num_classes_olay)

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

        ## Save meta-data for the extracted snapshots to a pickle file
        pickle.dump(img_info_all, open( fname_img_info_all, "wb" ) )

        logger.info('  Created snapshots for ' + str(len(img_info_all)) + ' subjects')

        ## Create and save output QC dataframe for image information
        df_qc_images = pd.DataFrame(data = {'ScanID' : df_images.ScanID.tolist(), 
                                            'qc_ok_flag' : qc_ok_flag_all,
                                            'qc_message' : qc_msg_all}, columns=['ScanID', 'qc_ok_flag', 'qc_message'])
        
        #df_qc_images = 
        #df_images[['ScanID']].copy()
        #df_qc_images['qc_ok_flag'] = qc_ok_flag_all
        #df_qc_images['qc_message'] = qc_msg_all
        try:
            fname_out = os.path.join(out_dir, 'log_qc_images_all.csv')
            df_qc_images.to_csv(fname_out, index = False)
            logger.info('  QC log saved to: ' + fname_out)
            
            fname_out = os.path.join(out_dir, 'log_qc_images_fail.csv')
            df_qc_images[df_qc_images.qc_ok_flag == 0].to_csv(fname_out, index = False)
            logger.info('  QC log for failed subjects saved to: ' + fname_out)

        except:
            logger.warning('  Could not save QC log files for img reading: ' + fname_out)

    return img_info_all


def create_html_report(params, out_dir, dir_subjects_full, dir_snapshots, dir_snapshots_full, dir_subjects, 
img_info_all, out_report):
    '''Function to create html report from pre-saved image snapshots and meta-data.
    
    :param params: All parameters
    :param out_dir: Output directory
    :param dir_subjects_full: Output directory for subjects (full path)
    :param dir_snapshots: Output directory for snapshots
    :param dir_snapshots_full: Output directory for snapshots (full path)
    :param dir_subjects: Output directory for subjects
    :param img_info_all: Meta data about snapshots
    :param out_report: Output report file name
    '''

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

    logger.info('  Report created: ' + out_report)


def create_report(list_file, config_file, out_dir):
    """Function to create the QC report.
    
    :param list_file: Image list file
    :param config_file: Configuration file
    :param out_dir: Output directory
    """

    ## Get the path for utils
    path_root = os.path.abspath(os.path.dirname(__file__))
    path_templates = os.path.join(path_root, 'js_templates')

    ### Check output file
    out_report = os.path.join(out_dir, 'qcreport.html')
    if os.path.exists(out_report):
        logger.warning('  Output report ' + out_report + ' already exists, skipping.')
        sys.exit();

    ## Read input file list
    logger.info('  Reading image list from: ' + list_file)
    try:
        df_images = pd.read_csv(list_file)
        df_images = df_images.fillna('') 
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
    #logger.info('     User config data: ' + '\n' + str(df_conf))

    ## Extract and check params
    logger.info('  Checking config parameters ...')    
    params = parse_config(df_conf, list_col_names)
    
    logger.info('    Parameters for QC report (missing values are updated with default ones): ' )
    logger.info('\n' + params.to_string())

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
    img_info_all = create_snapshots(params, df_images, dir_snapshots_full, out_dir)
    
    ## No image that passes the QC for creating snapshots
    if len(img_info_all) == 0:
        logger.info('  No images to display, skipping report creation ...')
    
    else:
        ### Create report
        logger.info('  Creating html report ...')
        create_html_report(params, out_dir, dir_subjects_full, dir_snapshots, 
                        dir_snapshots_full, dir_subjects, img_info_all, out_report)
    
#if __name__ == "__main__":
def main():
    """Script to generate the QC report for an input image dataset.

    :param outdir: Output report directory (full or relative path)
    """
    
    descr = 'Script to generate the QC report for an input image dataset.\n\n' \
            'The output QC report is an html file that displays snapshots of underlay ' \
            'and overlay (optional) images. \n\n' \
            'Before running the script, users should create an output folder (OUTDIR), ' \
            'and create two files in it: \n' \
            '- list_images.csv: List of underlay (and optionally mask and overlay) image files. \n' \
            '- config.csv: a configuration file that includes user parameters and their values.\n\n' \
            'Users can use the script "mrisnapshot_prep_data" to prepare these two files, and edit ' \
            'them for their specific dataset and parameter selection.\n\n' \
            'See https://cbica.github.io/MRISnapshot for usage and examples.\n\n'

    ## Create parser
    parser = argparse.ArgumentParser(add_help=False,
                                     prog="mrisnap_create_report",
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

    required.add_argument("-d", dest="outdir", type=str, required=True, help="Output directory")

    ## Parse input args
    args = parser.parse_args()
    
    logger.info('-----------------------------------------')    
    logger.info('Running : ' + ' '.join(sys.argv))    
    logger.info('-----------------------------------------')    

    ## Derive args
    outdir = os.path.abspath(args.outdir)    
    list_file = os.path.join(outdir, 'list_images.csv')
    config_file = os.path.join(outdir, 'config.csv')
    report_dir = os.path.join(outdir, 'QCReport')
    
    ## Make out directory
    if os.path.exists(report_dir) == False:
        os.makedirs(report_dir)

    ## Create report
    create_report(list_file, config_file, report_dir)
    
    
