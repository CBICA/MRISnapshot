#!/usr/bin/env python

### Import modules
import getopt
import numpy as np
import scipy as sp
import pandas as pd
import pickle
from PIL import Image, ImageFilter
import os, sys, time, shutil
from scipy.ndimage.interpolation import zoom
import MRISnapshot.utils.img_overlays as imolay
import MRISnapshot.utils.html_utils as html

import nibabel as nib
from nibabel.orientations import axcodes2ornt, ornt_transform, inv_ornt_aff


## Set logger  ## FIXME to be updated
import logging
format='%(levelname)-8s [%(filename)s : %(lineno)d - %(funcName)20s()] %(message)s'
format='%(levelname)-8s %(message)s'
logging.basicConfig(level=logging.DEBUG, format = '\n' + format, datefmt='%Y-%m-%d:%H:%M:%S')
logger = logging.getLogger(__name__)

##logger.setLevel(logging.DEBUG)      ## While debugging
logger.setLevel(logging.INFO)    ## FIXME Debug comments will be removed in release version
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s"


def check_params(params, list_col_names):
    '''Verify input arguments
    '''
    
    ### Check underlay column in file list
    if params.ulay_col not in list_col_names:
        sys.exit("Underlay column missing in image list")

    ### Check overlay columns
    params['num_olay'] = 0
    if params.olay_col in list_col_names:
        params['num_olay'] = 1
    if params.olay_col2 in list_col_names:
        params['num_olay'] = 2

    ### Check mask column
    params['num_mask'] = 0
    if params.mask_col in list_col_names:
        params['num_mask'] = 1

    ### Set additional params
    params['img_width_single'] = 1000  ### FIXME

    ### Convert values for view plane to a list
    if (params.view_plane == ''):
        params.view_plane = []
    else:
        params.view_plane = [n for n in params.view_plane.split('+')]

    ### Convert selected values for overlay to a list
    if (params.sel_vals_olay == ''):
        params.sel_vals_olay = []
    else:
        params.sel_vals_olay = [int(n) for n in params.sel_vals_olay.split('+')]
        
    ### Convert selected values for overlay2 to a list
    if params.sel_vals_olay2 == '':
        params.sel_vals_olay = []
    else:
        params.sel_vals_olay2 = [int(n) for n in params.sel_vals_olay2.split('+')]

    ### Update few params
    if params.is_edge == 1:
        params.is_transparent = 1
        logger.warning('is_edge selected for overlay, transparency reset to 1')

    return params

def create_dir(dir_name):
    '''Create output directories
    '''
    
    try:        
        if not os.path.exists(dir_name):
                os.makedirs(dir_name)
    except:
        sys.exit("Could not create out folder !!!");


def create_log_files(outdir):
    '''Create log files
    '''
    #logFile = outdir + os.sep + 'log_' + EXEC_NAME + '_' + startTimePretty + '.stdout'
    #errFile = outdir + os.sep + 'log_' + EXEC_NAME + '_' + startTimePretty + '.stderr'
    writeLog(logFile, '''------------------------''')

def copy_js(path_utils, out_dir):
    '''Copy js scripts
    '''
    
    shutil.copy(os.path.join(path_utils, 'save_qcform.js'), out_dir)
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

def calc_sel_slices():
    '''Select slices that will be used to create snapshots
    '''
    slices = [50, 70, 90, 110, 130, 150]      ## FIXME
    return slices

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
    
    logger.info('Img scaled to : ' + str(scale_low) + '  ' + str(scale_high))
    return nii_out

def extract_snapshot(img_ulay, img_olay, img_olay2, params, curr_view, curr_slice, slice_index, 
                     sub_id, dir_snapshots_full, list_sel_slices):
    
    # Get underlay slice
    img2d_ulay = img_ulay[:,:,curr_slice].astype(float)
    
    # Scale underlay image between 0 and 1
    img2d_ulay = (img2d_ulay - img2d_ulay.min()) / (img2d_ulay.max() - img2d_ulay.min())

    ## Resize underlay slice
    #img2d_ulay = zoom(img2d_ulay, (scX,scY), order=1)
    
    params.num_olay = 1
    
    # Create final images and save
    if params.num_olay == 0:
        pil_under = imolay.singleImage(img2d_ulay)
        outimg_noolay_suffix = sub_id + '_orient_' + curr_view + '_slice_' + str(slice_index)
        pil_under.convert('RGB').save(dir_snapshots_full + os.sep + outimg_noolay_suffix + '.png')
        outimg_witholay_suffix=''

    if params.num_olay == 1:
        
        img2d_olay = img_olay[:,:,curr_slice].astype(float)
        pil_under , pil_fused = imolay.overlayImage(img2d_ulay, img2d_olay, 
                                                    params.is_transparent, params.is_edge)

        outimg_noolay_suffix = sub_id + '_orient_' + curr_view + '_slice_' + str(slice_index)
        pil_under.convert('RGB').save(dir_snapshots_full + os.sep + outimg_noolay_suffix + '.png')

        outimg_witholay_suffix = sub_id + '_orient_' + curr_view + '_slice_' + str(slice_index) + '_witholay'
        pil_fused.convert('RGB').save(dir_snapshots_full + os.sep + outimg_witholay_suffix + '.png')


    # Keep image information for later creation of html files
    outimg_caption = 'Slice: ' + curr_view + '_' + str(list_sel_slices[slice_index] + 1)
    
    return [outimg_caption, outimg_noolay_suffix, outimg_witholay_suffix]

def create_snapshots(params, df_images, dir_snapshots_full):
    '''Create snapshots and meta data about them
    '''
    
    # Dictionary with img orientation for different views
    d_orient = {'A':'PLS', 'S':'IPR', 'C':'IRP'}  
    
    ### Extract and save snapshots with metadata
    if not os.path.isfile(dir_snapshots_full + os.sep + 'img_info_all.pickle'):

        img_info_all = [];
        for sub_index, sub_id in enumerate(df_images[params.id_col]):

            ### Read input images
            nii_ulay, fname_ulay  = get_nifti(df_images, sub_index, params.ulay_col)
            nii_mask, fname_mask  = get_nifti(df_images, sub_index, params.mask_col)
            nii_olay, fname_olay  = get_nifti(df_images, sub_index, params.olay_col)
            nii_olay2, fname_olay2  = get_nifti(df_images, sub_index, params.olay_col2)

            # Initialize containers to keep image info
            outimg_noolay_suffix_all = []
            outimg_witholay_suffix_all = []
            outimg_caption_all = []
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
                list_sel_slices = calc_sel_slices()
                for slice_index, curr_slice in enumerate(list_sel_slices):
                    info_snapshot = extract_snapshot(img3d_ulay, img3d_olay, img3d_olay2, params, 
                                                     curr_view, curr_slice, slice_index, sub_id, 
                                                     dir_snapshots_full, list_sel_slices)
                    outimg_caption_all.append(info_snapshot[0])
                    outimg_noolay_suffix_all.append(info_snapshot[1])
                    outimg_witholay_suffix_all.append(info_snapshot[2])

                list_sel_slices_all.append(list_sel_slices)
                    
            ### Keep image information for later creation of html files
            snapshot_info = {'sub_index' : sub_index, 'sub_id' : sub_id, 
                              'fname_ulay' : fname_ulay, 'fname_olay' : fname_olay, 
                              'fname_olay2' : fname_olay2, 'views' : params.view_plane, 
                              'list_sel_slices_all' : list_sel_slices_all, 
                              'outimg_noolay_suffix_all' : outimg_noolay_suffix_all, 
                              'outimg_witholay_suffix_all' : outimg_witholay_suffix_all, 
                              'outimg_caption_all' : outimg_caption_all}
                
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
        outimg_noolay_suffix_all = item['outimg_noolay_suffix_all']
        outimg_witholay_suffix_all = item['outimg_witholay_suffix_all']
        outimg_caption_all = item['outimg_caption_all']

        # Header for file with multiple snapshots
        HTML_multi_snapshot = html.htmlSubjectPrefix(sub_id, i+1, numItem, fname_ulay, 
                                                     fname_olay, fname_olay2)
            
        for j, outimg_noolay_suffix in enumerate(outimg_noolay_suffix_all):
            outimg_witholay_suffix = outimg_witholay_suffix_all[j]
            outimg_caption = outimg_caption_all[j]
            
            # Write html page for each snapshot without overlay
            tmp_img = outimg_noolay_suffix + '.png'
            tmp_txt = 'Fig: ' + str(j) + outimg_caption
            tmp_html = outimg_witholay_suffix + '.html'            
            HTML_single_snapshot = html.htmlSnapshot(tmp_img, tmp_txt, tmp_html, params.img_width_single) 

            ofp = open(os.path.join(dir_snapshots_full, outimg_noolay_suffix + '.html'), 'w')
            ofp.write(HTML_single_snapshot)
            ofp.close

            # Write html page for each snapshot
            if params.num_olay == 0:
                
                # Append html code for all snapshots from a subject
                tmp_img = os.path.join(dir_snapshots, outimg_noolay_suffix + '.png')
                tmp_html = os.path.join(dir_snapshots, outimg_noolay_suffix + '.html')

                HTML_tmp = html.htmlSubjectAddImage(tmp_img, outimg_caption, tmp_html)

            else:
                tmp_img = outimg_witholay_suffix + '.png'
                tmp_txt = 'Fig: ' + str(j) + outimg_caption
                tmp_html = outimg_noolay_suffix + '.html'
                HTML_single_snapshot = html.htmlSnapshot(tmp_img, tmp_txt, tmp_html, params.img_width_single)
                
                tmp_html = os.path.join(dir_snapshots_full, outimg_witholay_suffix + '.html')
                
                ofp = open(tmp_html, 'w')
                ofp.write(HTML_single_snapshot)
                ofp.close

                # Append html code for all snapshots from a subject
                tmp_img = os.path.join(dir_snapshots, outimg_witholay_suffix + '.png')
                tmp_html = os.path.join(dir_snapshots, outimg_witholay_suffix + '.html')
                HTML_tmp = html.htmlSubjectAddImage(tmp_img, outimg_caption, tmp_html)
            
            HTML_multi_snapshot = HTML_multi_snapshot + '\n' + HTML_tmp

        # Prepare html page for the subject (with all snapshots)
        html_sub = sub_id + '_mriqc.html'
        HTML_multi_snapshot = HTML_multi_snapshot + '\n' + '<p style="clear: both;"><br>'
        
        if params.is_out_noqc == 0:
            # Add qc form to html page for the subject
            HTML_multi_snapshot = HTML_multi_snapshot + '\n' + html.htmlQCForm('EXC-Img', 'EXC-Proc', 'Notes',)

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

    ## Read input file list
    try:
        df_images = pd.read_csv(list_file)
        list_col_names = df_images.columns.values
    except:
        sys.exit("Could not read list file: " + list_file);

    ## Read params from config file
    try:
        df_conf = pd.read_csv(config_file).fillna('')
    except:
        sys.exit("Could not read config file: " +  config_file);
    params = df_conf.set_index('ParamName').ParamValue

    ## Convert numeric params
    df_tmp = pd.to_numeric(params, errors='coerce')
    params = df_tmp.combine_first(params)

    ## Verify params
    params = check_params(params, list_col_names)

    ### Check output file
    out_report = os.path.join(out_dir, 'qcreport.html')
    if os.path.exists(out_report):
        sys.exit("Output report exists: " + out_report);

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
    copy_js(path_templates, dir_scripts)
    
    ### Create snapshots
    img_info_all = create_snapshots(params, df_images, dir_snapshots_full)
    
    ### Create report
    create_html_report(params, out_dir, dir_subjects_full, dir_snapshots, 
                       dir_snapshots_full, dir_subjects, img_info_all, out_report)
    
    
    
