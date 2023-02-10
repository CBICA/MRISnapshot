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


def check_params(params, file_types):
        
    ### Check params
    if params.Underlay not in file_types:
        sys.exit("The underlay image type (" + params.Underlay + ") is not found in the file list (" + param_filelist + "); please check the headers in the file")

    params['NumOverlay'] = 0
    if params.Overlay  != '':
        params.NumOverlay = 1
        if params.Overlay not in file_types:
            sys.exit("The overlay image type (" + params.Overlay + ") is not found in the file list (" + param_filelist + "); please check the headers in the file")

    if params.Overlay2  != '':
        params.NumOverlay = 2
        if params.Overlay2 not in file_types:
            sys.exit("The overlay2 image type (" + params.Overlay2 + ") is not found in the file list (" + param_filelist + "); please check the headers in the file")

    params['IsMask'] = 0
    if params.Mask != '':
        params.IsMask = 1
        if params.Mask not in file_types:
            sys.exit("The mask image type (" + params.Mask + ") is not found in the file list (" + param_filelist + "); please check the headers in the file")

    ### Set additional params
    params['ImgWidthSingle'] = 1000  ### FIXME

    params['ressize'] = 1 ### FIXME

    if params.SelValsOverlay == '':
        params.SelValsOverlay = []
    else:
        params.SelValsOverlay = [int(n) for n in params.SelValsOverlay.str.split('+')]
        
    if params.SelValsOverlay2 == '':
        params.SelValsOverlay = []
    else:
        params.SelValsOverlay2 = [int(n) for n in params.SelValsOverlay2.str.split('+')]

    if params.IsEdge == 1:
        params.Transp = 1
        print("Warning: IsEdge selected for overlay; transparency reset to 1")

    logger.info(params)

    return params

def create_dir(dir_name):
    ### Create output directories
    try:        
        if not os.path.exists(dir_name):
                os.makedirs(dir_name)
    except:
        sys.exit("Could not create out folder !!!");


def create_log_files(outdir):
    ### Create log files
    #logFile = outdir + os.sep + 'log_' + EXEC_NAME + '_' + startTimePretty + '.stdout'
    #errFile = outdir + os.sep + 'log_' + EXEC_NAME + '_' + startTimePretty + '.stderr'
    writeLog(logFile, '''------------------------''')

def copy_js(outdir):
    ### Copy js scripts
    print('TODO')
    #shutil.copy(EXEC_DIR + os.sep + 'utils' + os.sep + 'saveqcform.js', dir_scripts)
    #shutil.copy(EXEC_DIR + os.sep + 'utils' + os.sep + 'miscfunc.js', dir_scripts)
    #shutil.copy(EXEC_DIR + os.sep + 'utils' + os.sep + 'shortcut.js', dir_scripts)
    #shutil.copy(EXEC_DIR + os.sep + 'utils' + os.sep + 'loadback.js', dir_scripts)

def get_nifti(fname, orient = 'LPS'):
    ''' Read nifti image and reorient 
    '''
    nii = nib.load(fname)
    orig_ornt = nib.io_orientation(nii.affine)
    targ_ornt = axcodes2ornt(orient)
    if np.all(orig_ornt == targ_ornt) == False:
        transform = ornt_transform(orig_ornt, targ_ornt)
        nii = nii.as_reoriented(transform)
    return nii

def get_img_mat(nii, orient = 'LPS'):
    ''' Reorient nifti and get data matrix
    '''
    orig_ornt = nib.io_orientation(nii.affine)
    targ_ornt = axcodes2ornt(orient)
    if np.all(orig_ornt == targ_ornt) == False:
        transform = ornt_transform(orig_ornt, targ_ornt)
        nii = nii.as_reoriented(transform)
    return nii.get_fdata()

def create_snapshots(params, df_files, dir_snapshots_full, labelID):
    
    # Dictionary with img orientation for different views
    d_orient = {'A':'LPS', 'S':'SLP', 'C':'PSL'}  
    
    params.NumOverlay = 0 ## FIXME
    
    
    ### Extract and save snapshots with metadata
    if not os.path.isfile(dir_snapshots_full + os.sep + 'mriqc_img_info_all.pickle'):

        mriqc_img_info_all = [];
        for sub_index, sub_id in enumerate(df_files[labelID]):

            ### Read input images
            fname_under = df_files.loc[sub_index][params.Underlay]
            nii_ulay = get_nifti(fname_under)

            fname_over = df_files.loc[sub_index][params.Overlay]
            nii_olay = get_nifti(fname_over)

            fname_over2 = 'None'    ## FIXME
            
            fname_mask = df_files.loc[sub_index][params.Mask]
            nii_mask = get_nifti(fname_mask)

            # Initialize containers to keep image info
            outimg_noolay_suffix_all = []
            outimg_witholay_suffix_all = []
            outimg_caption_all = []
            sel_slices_all = []

            ### Create snapshots for each orientation
            for view_index, view in enumerate(params.View):

                ### Get data in selected orientation
                img3d_ulay = get_img_mat(nii_ulay, d_orient[view]) 
                img3d_olay = get_img_mat(nii_olay, d_orient[view]) 
                img3d_mask = get_img_mat(nii_mask, d_orient[view]) 
            
                ### Select slices to show
                sel_slices = [50, 70, 90, 110]      ## FIXME
                for slice_index, curr_slice in enumerate(sel_slices):

                    # Get underlay slice
                    img2d_ulay = img3d_ulay[:,:,curr_slice].astype(float)
                    
                    # Scale underlay image between 0 and 1
                    img2d_ulay = (img2d_ulay - img2d_ulay.min()) / (img2d_ulay.max() - img2d_ulay.min())

                    ## Resize underlay slice
                    #img2d_ulay = zoom(img2d_ulay, (scX,scY), order=1)
                    
                    # Create final images and save
                    if params.NumOverlay == 0:
                        pil_under = imolay.singleImage(img2d_ulay)
                        outimg_noolay_suffix = sub_id + '_orient_' + view + '_slice_' + str(slice_index)
                        pil_under.convert('RGB').save(dir_snapshots_full + os.sep + outimg_noolay_suffix + '.png')
                        outimg_witholay_suffix=''

                    # Keep image information for later creation of html files
                    outimg_caption = 'Slice: ' + view + '_' + str(sel_slices[slice_index] + 1)
                    outimg_caption_all.append(outimg_caption)
                    
                    outimg_noolay_suffix_all.append(outimg_noolay_suffix)
                    outimg_witholay_suffix_all.append(outimg_witholay_suffix)

                sel_slices_all.append(sel_slices)
                    
            ### Keep image information for later creation of html files
            mriqc_img_info = {'sub_index' : sub_index, 'sub_id' : sub_id, 'fname_under' : fname_under, 'fname_over' : fname_over, 'fname_over2' : fname_over2, 'views' : params.View, 'sel_slices_all' : sel_slices_all, 'outimg_noolay_suffix_all' : outimg_noolay_suffix_all, 'outimg_witholay_suffix_all' : outimg_witholay_suffix_all, 'outimg_caption_all' : outimg_caption_all}
                
            mriqc_img_info_all.append(mriqc_img_info)

        pickle.dump( mriqc_img_info_all, open( dir_snapshots_full + os.sep + 'mriqc_img_info_all.pickle', "wb" ) )
        
    ### If slices were already extracted, use saved meta data
    else:
        mriqc_img_info_all = pickle.load( open( dir_snapshots_full + os.sep + 'mriqc_img_info_all.pickle', "rb" ) )

    return mriqc_img_info_all


def create_html_report(params, dir_subjects_full, dir_snapshots, dir_snapshots_full, dir_subjects, mriqc_img_info_all, outReportName):
    ###################################################################
    ### CREATE HTML REPORTS ###########################################

    ### Create html files for: 
    ###        - single snapshots with/without overlay for each selected slice, 
    ###        - all snapshots from a subject together
    ###        - all subjects together

    # Write the stylesheet
    TEXT_HTML_STR_STYLESHEET = html.html_stylesheet(params.ImgWidth);
    ofp = open(dir_subjects_full + os.sep + 'scripts' + os.sep + 'pagestyle.css','w')
    ofp.write(TEXT_HTML_STR_STYLESHEET)
    ofp.close

    # Prepare text for final report
    TEXT_HTML_STR_MULTI_SUBJ = ''

    numItem = len(mriqc_img_info_all)

    html_final = dir_subjects_full + os.sep + 'finalPage.html'

    for itemindex,item in enumerate(mriqc_img_info_all):

        ### Read meta-data about snapshots
        sub_id = item['sub_id']
        prev_sub_id = mriqc_img_info_all[np.max([0,itemindex-1])]['sub_id']
        next_sub_id = mriqc_img_info_all[np.min([numItem-1,itemindex+1])]['sub_id']

        fname_under = item['fname_under']
        fname_over = item['fname_over']
        fname_over2 = item['fname_over2']
        outimg_noolay_suffix_all = item['outimg_noolay_suffix_all']
        outimg_witholay_suffix_all = item['outimg_witholay_suffix_all']
        outimg_caption_all = item['outimg_caption_all']

        # Header for file with multiple snapshots
        TEXT_HTML_STR_MULTI_IMG = html.htmlSubjectPrefix(sub_id, itemindex+1, numItem, fname_under, fname_over, fname_over2)
            
        for imgindex, outimg_noolay_suffix in enumerate(outimg_noolay_suffix_all):
            
            outimg_witholay_suffix = outimg_witholay_suffix_all[imgindex]
            
            outimg_caption = outimg_caption_all[imgindex]
            
            # Write html page for each snapshot without overlay
            TEXT_HTML_STR_SINGLE_IMG = html.htmlSnapshot(outimg_noolay_suffix + '.png', 'Fig: ' + str(imgindex) + outimg_caption, outimg_witholay_suffix + '.html', params.ImgWidthSingle);
            ofp = open(dir_snapshots_full + os.sep + outimg_noolay_suffix + '.html','w')
            ofp.write(TEXT_HTML_STR_SINGLE_IMG)
            ofp.close

            # Write html page for each snapshot
            if params.NumOverlay == 0:
                # Append html code for all snapshots from a subject
                TEXT_HTML_STR_TMP = html.htmlSubjectAddImage(dir_snapshots + os.sep + outimg_noolay_suffix + '.png', outimg_caption, dir_snapshots + os.sep + outimg_noolay_suffix + '.html')

            else:
                TEXT_HTML_STR_SINGLE_IMG = html.htmlSnapshot(outimg_witholay_suffix + '.png', 'Fig: ' + str(imgindex) + outimg_caption, outimg_noolay_suffix + '.html', params.ImgWidthSingle);
                ofp = open(dir_snapshots_full + os.sep + outimg_witholay_suffix + '.html','w')
                ofp.write(TEXT_HTML_STR_SINGLE_IMG)
                ofp.close

                # Append html code for all snapshots from a subject
                TEXT_HTML_STR_TMP = html.htmlSubjectAddImage(dir_snapshots + os.sep + outimg_witholay_suffix + '.png', outimg_caption, dir_snapshots + os.sep + outimg_witholay_suffix + '.html')
            
            TEXT_HTML_STR_MULTI_IMG = TEXT_HTML_STR_MULTI_IMG + '\n' + TEXT_HTML_STR_TMP

        # Prepare html page for the subject (with all snapshots)
        html_sub = sub_id + '_mriqc.html'
        TEXT_HTML_STR_MULTI_IMG = TEXT_HTML_STR_MULTI_IMG + '\n' + '<p style="clear: both;"><br>'
        if params.IsOutNoQC == 0:
            # Add qc form to html page for the subject
            TEXT_HTML_STR_MULTI_IMG = TEXT_HTML_STR_MULTI_IMG + '\n' + html.htmlQCForm('EXC-Img', 'EXC-Proc', 'Notes',)

        if params.IsOutSingle == 1:
            # Prepare html page for all subjects together
            TEXT_HTML_STR_MULTI_SUBJ = TEXT_HTML_STR_MULTI_SUBJ + '\n' + TEXT_HTML_STR_MULTI_IMG
        else:
            # Add navigation to html page for the subject
            html_next_sub = next_sub_id + '_mriqc.html'
            html_prev_sub = prev_sub_id + '_mriqc.html'

            TEXT_HTML_STR_MULTI_IMG = html.html_Navig(html_prev_sub, html_next_sub) + TEXT_HTML_STR_MULTI_IMG + '\n' + html.html_Navig(html_prev_sub, html_next_sub)

        # Write html page for the subject
        ofp = open(dir_subjects_full + os.sep + html_sub,'w')
        ofp.write(TEXT_HTML_STR_MULTI_IMG)
        ofp.close()
        

    ### Write the output report
    if params.IsOutSingle == 1:
        # Write html page for all subjects together (final mriqc report)
        html_sub = dir_subjects_full + os.sep + 'ALLSUBJ_mriqc.html'

        ofp = open(html_sub, 'w')
        ofp.write(TEXT_HTML_STR_MULTI_SUBJ)
        ofp.close()
        
        file_report = dir_subjects + os.sep + 'ALLSUBJ_mriqc.html'

    else:
        # Create link to the report for the first subject
        first_sub_id = mriqc_img_info_all[0]['sub_id']
        file_report = dir_subjects + os.sep + first_sub_id + '_mriqc.html'

    TEXT_HTML_MAINPAGE = html.htmlMainPage(file_report)
    ofp = open(outReportName, 'w')
    ofp.write(TEXT_HTML_MAINPAGE)
    ofp.close()


def create_report(param_filelist, param_config, param_outdir):

    ### Read input file list
    try:
        df_files = pd.read_csv(param_filelist)
    except:
        sys.exit("Could not read list file (" +  param_filelist + "). Aborting operations !!!");
    file_types = df_files.columns.values
    labelID = file_types[0]

    ### Read params from config file
    try:
        dfConf = pd.read_csv(param_config, comment='#', index_col='PARAM_NAME').fillna('')
    except:
        sys.exit("Could not read config file (" +  param_config + "). Aborting operations !!!");
    params = dfConf[dfConf.columns[0]].T
    
    
    df_tmp = pd.to_numeric(params, errors='coerce')
    params = df_tmp.combine_first(params)

    #logger.info(params)
    #input()

    params = check_params(params, file_types)


    ### Check output file
    outReportName = param_outdir + os.sep + 'mriqc_report.html'
    if os.path.exists(outReportName):
        sys.exit("Output report exists, delete it first and rerun: " + outReportName);

    ### Create out dirs
    create_dir(param_outdir)    
    dir_subjects = 'subjects'
    dir_subjects_full = param_outdir + os.sep + dir_subjects

    dir_snapshots = 'snapshots'
    dir_snapshots_full = dir_subjects_full + os.sep + dir_snapshots

    create_dir(dir_subjects_full)
    create_dir(dir_snapshots_full)

    dir_scripts = dir_subjects_full + os.sep + 'scripts'
    create_dir(dir_scripts)
    
    ### Create snapshots
    mriqc_img_info_all = create_snapshots(params, df_files, dir_snapshots_full, labelID)
    
    ### Create report
    create_html_report(params, dir_subjects_full, dir_snapshots, dir_snapshots_full, dir_subjects, mriqc_img_info_all, outReportName)
    
    
    
