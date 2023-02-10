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


def check_params(params, file_types):
    ### Set values for optional params
    if 'Overlay' not in params:
        params.Overlay = ''

    if 'Overlay2' not in params:
        params.Overlay2 = ''

    if 'Mask' not in params:
        params.Mask = ''

    if 'SelValsOverlay' not in params:
        params.SelValsOverlay = ''

    if 'SelValsOverlay2' not in params:
        params.SelValsOverlay2 = ''

    if 'isEdge' not in params:
        params.isEdge = 0

    if 'BinOverlay' not in params:
        params.BinOverlay = 0

    if 'MinVox' not in params:
        params.MinVox = 0
        
    if 'MinVox' not in params:
        params.MinVox = 0

    if 'IsEdge' not in params:
        params.IsEdge = 0
        
    if 'Transp' not in params:
        params.Transp = 1
        
    if 'PercHigh' not in params:
        params.PercHigh = 100

    if 'PercLow' not in params:
        params.PercLow = 0
        
    if 'IsOutSingle' not in params:
        params.IsOutSingle = 0
        
    if 'IsOutNoQC' not in params:
        params.IsOutNoQC = 0
        
    if 'ImgWidth' not in params:
        params.ImgWidth = 100
        
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

    return params

def create_out_dirs(outdir):
    ### Create output directories
    try:        
        dir_subjects = 'subjects'
        dir_subjects_full = outdir + os.sep + dir_subjects

        dir_snapshots = 'snapshots'
        dir_snapshots_full = dir_subjects_full + os.sep + dir_snapshots

        if outdir and not os.path.exists(outdir):
                os.makedirs(outdir)

        if dir_subjects_full and not os.path.exists(dir_subjects_full):
                os.makedirs(dir_subjects_full)

        if dir_snapshots_full and not os.path.exists(dir_snapshots_full):
                os.makedirs(dir_snapshots_full)

        dir_scripts = dir_subjects_full + os.sep + 'scripts'
        if dir_scripts and not os.path.exists(dir_scripts):
                os.makedirs(dir_scripts)
    except:
        sys.exit("Could not create output folders. Aborting operations !!!");


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


def create_snapshots():
    ### If slices were not already extracted, extract them and save meta data
    if not os.path.isfile(dir_snapshots_full + os.sep + 'mriqc_img_info_all.pickle'):

        writeLog(logFile,'''------------------------ 
    Extracting snapshots ... (''' + time.strftime("%Y-%m-%d_%H-%M-%S") + ')') 

        mriqc_img_info_all = [];

        for sub_index, sub_id in enumerate(dfFiles[labelID]):

            writeLog(logFile,"        Sub " + str(sub_index) + ": " + sub_id) 

            ### Read the input images
            try:
                fname_under = dfFiles.loc[sub_index][params.Underlay]
                nii_under = nib.load( os.path.join( fname_under ) )
                img3d_under = nii_under.get_fdata()
                hdr_under = nii_under.header
                
                d1,d2,d3 = img3d_under.shape
                r1,r2,r3 = nii_under.header.get_zooms()
                sc1 = r1 / params.ressize
                sc2 = r2 / params.ressize
                sc3 = r3 / params.ressize
                
                print(img3d_under.shape)
                print(nii_under.header.get_zooms())
                
            except:
                msg="Could not read underlay image for subject " + str(sub_index) + ":" + sub_id
                writeLog(errFile, msg)
                continue

            fname_over = ''
            fname_over2 = ''
            if params.NumOverlay > 0:
                try:
                    fname_over = dfFiles.loc[sub_index][params.Overlay]
                    nii_over = nib.load( os.path.join( fname_over ) )
                    img3d_over = nii_over.get_fdata()
                    hdr_over = nii_over.header
                except:
                    msg="Could not read overlay image for subject " + str(sub_index) + ":" + sub_id
                    writeLog(errFile, msg)
                
            if params.NumOverlay > 1:
                try:
                    fname_over2 = dfFiles.loc[sub_index][params.Overlay2]
                    nii_over2 = nib.load( os.path.join( fname_over2 ) )
                    img3d_over2 = nii_over2.get_fdata()
                    hdr_over2 = nii_over2.header
                except:
                    msg="Could not read overlay image for subject " + str(sub_index) + ":" + sub_id
                    writeLog(errFile, msg)

            img3d_mask = img3d_under * 0 + 1    ## Mask with all voxels set to 1
            if params.IsMask == 1:
                try:
                    fname_mask = dfFiles.loc[sub_index][params.Mask]
                    nii_mask = nib.load( os.path.join( fname_mask ) )
                    img3d_mask = nii_mask.get_fdata()
                    img3d_mask = (img3d_mask>0).astype(int)
                    hdr_mask = nii_mask.header
                except:
                    msg="Could not read mask image for subject " + str(sub_index) + ":" + sub_id
                    writeLog(errFile, msg)
                    continue
                
            img3d_under_masked = img3d_under * img3d_mask


            print(fname_under)
            print(fname_over)

            # Scale underlay image between min and max percentile values
            scaleCutoffs = np.percentile(img3d_under_masked[np.nonzero(img3d_under_masked)],[params.PercLow, params.PercHigh])
            img3d_under[np.where(img3d_under<scaleCutoffs[0])] = scaleCutoffs[0]
            img3d_under[np.where(img3d_under>scaleCutoffs[1])] = scaleCutoffs[1]
            img3d_under = img3d_under - scaleCutoffs[0]

            # Select values in the overlay image
            if len(params.SelValsOverlay) > 0:
                    img3d_over[np.isin(img3d_over, params.SelValsOverlay)==False]=0

            # Select values in the overlay image
            if len(params.SelValsOverlay2) > 0:
                    img3d_over2[np.isin(img3d_over2, params.SelValsOverlay2)==False]=0

            # Binarize overlay image
            if params.BinOverlay == 1:
                if params.NumOverlay > 0:
                    img3d_over = (img3d_over>0).astype(int)
                if params.NumOverlay > 1:
                    img3d_over2 = (img3d_over2>0).astype(int)

            # Initialize containers to keep image info
            outimg_noolay_suffix_all = []
            outimg_witholay_suffix_all = []
            outimg_caption_all = []
            sl_sel_all = []

            ### Create snapshots
            for view_index, view in enumerate(params.View):

                if view != 'A' and view != 'S' and view != 'C':
                    continue

                scX=sc1
                scY=sc2

                # Process ulay
                img3d_under_reshaped = img3d_under.copy();
                img3d_under_masked_reshaped = img3d_under_masked.copy();
                if view == 'S':                # Rotate image for sagittal or coronal views
                    scX=sc3
                    scY=sc2
                    img3d_under_tmp = np.transpose(img3d_under,(1,2,0))
                    img3d_under_reshaped = img3d_under_tmp[:,::-1,:]

                    img3d_under_tmp = np.transpose(img3d_under_masked,(1,2,0))
                    img3d_under_masked_reshaped = img3d_under_tmp[:,::-1,:]
                    
                if view == 'C':
                    scX=sc3
                    scY=sc1
                    img3d_under_tmp = np.transpose(img3d_under,(0,2,1))
                    img3d_under_reshaped = img3d_under_tmp[:,::-1,:]

                    img3d_under_tmp = np.transpose(img3d_under_masked,(0,2,1))
                    img3d_under_masked_reshaped = img3d_under_tmp[:,::-1,:]

                # Process olay
                
                if params.NumOverlay > 0:
                    img3d_over_reshaped = img3d_over.copy()
                    if view == 'S':                # Rotate image for sagittal or coronal views
                        img3d_over_tmp = np.transpose(img3d_over,(1,2,0))
                        img3d_over_reshaped = img3d_over_tmp[:,::-1,:]
                        
                    if view == 'C':
                        img3d_over_tmp = np.transpose(img3d_over,(0,2,1))
                        img3d_over_reshaped = img3d_over_tmp[:,::-1,:]

                # Process olay2
                if params.NumOverlay > 1:
                    img3d_over2_reshaped = img3d_over2.copy()
                    if view == 'S':                # Rotate image for sagittal or coronal views
                        img3d_over2_tmp = np.transpose(img3d_over2,(1,2,0))
                        img3d_over2_reshaped = img3d_over2_tmp[:,::-1,:]
                        
                    if view == 'C':
                        img3d_over2_tmp = np.transpose(img3d_over2,(0,2,1))
                        img3d_over2_reshaped = img3d_over2_tmp[:,::-1,:]
                
                
                ### Select slices to show
                img3d_sel_slice = img3d_under_masked_reshaped
                if params.NumOverlay == 1:
                    img3d_sel_slice = img3d_over_reshaped
                    
                    print('eee')
                    print(np.sum(img3d_sel_slice))
                    
                if params.NumOverlay == 2:
                    img3d_sel_slice = img3d_over_reshaped + img3d_over2_reshaped
                
                sl_nonzero = np.where(np.sum(np.sum(img3d_sel_slice,0),0) > params.MinVox)[0]
                num_nonzero = sl_nonzero.size
                
                if num_nonzero == 0:        # No slice to show; continue        
                    msg="No slice to show in view: " + params.View + ", for subject " + str(sub_index) + ":" + sub_id
                    writeLog(errFile, msg)
                    continue
                    
                if params.NumSlice < 0:        # Select using absolute value of params.NumSlice as step size
                    if -1*params.NumSlice > num_nonzero:        # Not enough slices, just show the one in middle
                        sl_sel = sl_nonzero[num_nonzero/2]
                    else:
                        sl_sel = sl_nonzero[np.arange(params.NumSlice/2.0, num_nonzero, -1*params.NumSlice)[1:].round().astype(int)]
                        
                else:                        # Select using params.NumSlice as total number of slices
                    if num_nonzero <= params.NumSlice:
                        sl_sel = sl_nonzero
                    else:
                        hstep = float(num_nonzero) / params.NumSlice
                        sl_sel = sl_nonzero[np.arange(-1.0 * hstep / 2.0, num_nonzero, hstep)[1:].round().astype(int)]

                if isinstance(sl_sel,np.int64):
                    sl_sel = np.array([sl_sel])

                for index, sl_sel_tmp in enumerate(sl_sel):

                    # Get underlay slice
                    img2d_under = np.array(img3d_under_reshaped[:,:,sl_sel_tmp].transpose()).astype(float)
                    
                    # Scale underlay image between 0 and 1
                    img2d_under = (img2d_under - img2d_under.min()) / (img2d_under.max() - img2d_under.min())

                    # Resize underlay slice
                    img2d_under = zoom(img2d_under, (scX,scY), order=1)
                    
                    # Create final images and save
                    if params.NumOverlay == 0:
                        pil_under = imolay.singleImage(img2d_under)
                        outimg_noolay_suffix = sub_id + '_orient_' + view + '_slice_' + str(index)
                        pil_under.convert('RGB').save(dir_snapshots_full + os.sep + outimg_noolay_suffix + '.png')
                        outimg_witholay_suffix=''
                        
                    elif params.NumOverlay == 1:
                        img2d_over = np.array(img3d_over_reshaped[:,:,sl_sel_tmp].transpose()).astype(float)
                        #img2d_over[np.where(img2d_over>0)]=1

                        # Resize overlay slice
                        img2d_over = zoom(img2d_over, (scX,scY), order=0)

                        pil_under , pil_fused = imolay.overlayImage(img2d_under, img2d_over, params.Transp, params.IsEdge)

                        outimg_noolay_suffix = sub_id + '_orient_' + view + '_slice_' + str(index)
                        pil_under.convert('RGB').save(dir_snapshots_full + os.sep + outimg_noolay_suffix + '.png')

                        outimg_witholay_suffix = sub_id + '_orient_' + view + '_slice_' + str(index) + '_witholay'
                        pil_fused.convert('RGB').save(dir_snapshots_full + os.sep + outimg_witholay_suffix + '.png')

                    elif params.NumOverlay == 2:
                        img2d_over = np.array(img3d_over_reshaped[:,:,sl_sel_tmp].transpose()).astype(float)
                        img2d_over2 = np.array(img3d_over2_reshaped[:,:,sl_sel_tmp].transpose()).astype(float)

                        # Resize overlay slice
                        img2d_over = zoom(img2d_over, (scX,scY), order=0)
                        img2d_over2 = zoom(img2d_over2, (scX,scY), order=0)

                        pil_under , pil_fused = imolay.overlayImageDouble(img2d_under, img2d_over, img2d_over2, params.Transp, params.IsEdge)

                        outimg_noolay_suffix = sub_id + '_orient_' + view + '_slice_' + str(index)
                        pil_under.convert('RGB').save(dir_snapshots_full + os.sep + outimg_noolay_suffix + '.png')

                        outimg_witholay_suffix = sub_id + '_orient_' + view + '_slice_' + str(index) + '_witholay'
                        pil_fused.convert('RGB').save(dir_snapshots_full + os.sep + outimg_witholay_suffix + '.png')

                    # Keep image information for later creation of html files
                    outimg_caption = 'Slice: ' + view + '_' + str(sl_sel[index] + 1)
                    outimg_caption_all.append(outimg_caption)
                    
                    outimg_noolay_suffix_all.append(outimg_noolay_suffix)
                    outimg_witholay_suffix_all.append(outimg_witholay_suffix)

                sl_sel_all.append(sl_sel)
                    
            ### Keep image information for later creation of html files
            mriqc_img_info = {'sub_index' : sub_index, 'sub_id' : sub_id, 'fname_under' : fname_under, 'fname_over' : fname_over, 'fname_over2' : fname_over2, 'views' : params.View, 'sl_sel_all' : sl_sel_all, 'outimg_noolay_suffix_all' : outimg_noolay_suffix_all, 'outimg_witholay_suffix_all' : outimg_witholay_suffix_all, 'outimg_caption_all' : outimg_caption_all}
                
            mriqc_img_info_all.append(mriqc_img_info)

        pickle.dump( mriqc_img_info_all, open( dir_snapshots_full + os.sep + 'mriqc_img_info_all.pickle', "wb" ) )
        
    ### If slices were already extracted, use saved meta data
    else:
        mriqc_img_info_all = pickle.load( open( dir_snapshots_full + os.sep + 'mriqc_img_info_all.pickle', "rb" ) )

        writeLog(logFile,'''------------------------Reading pre-calculated snapshots ... (''' + time.strftime("%Y-%m-%d_%H-%M-%S") + ')/n') 



def create_html_report():
    ###################################################################
    ### CREATE HTML REPORTS ###########################################

    ### Create html files for: 
    ###        - single snapshots with/without overlay for each selected slice, 
    ###        - all snapshots from a subject together
    ###        - all subjects together

    writeLog(logFile,'''------------------------ 
    Creating mriqc reports ... (''' + time.strftime("%Y-%m-%d_%H-%M-%S") + ')')

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

        writeLog(logFile,"        Sub " + str(itemindex) + ": " + sub_id) 
        
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
        
    writeLog(logFile,'''------------------------ 
    Mriqc reports created ... (''' + time.strftime("%Y-%m-%d_%H-%M-%S") + ')\n') 


def create_report(param_filelist, param_config, param_outdir):

    ### Read input file list
    try:
        dfFiles = pd.read_csv(param_filelist)
    except:
        sys.exit("Could not read list file (" +  param_filelist + "). Aborting operations !!!");
    file_types = dfFiles.columns.values
    labelID = file_types[0]

    ### Read params from config file
    try:
        dfConf = pd.read_csv(param_config, comment='#', index_col='PARAM_NAME').fillna('')
    except:
        sys.exit("Could not read config file (" +  param_config + "). Aborting operations !!!");
    params = dfConf[dfConf.columns[0]].T
    df_tmp = pd.to_numeric(params, errors='coerce')
    params = df_tmp.combine_first(params)

    params = check_params(params, file_types)

    ### Check output file
    outReportName = param_outdir + os.sep + 'mriqc_report.html'
    if os.path.exists(outReportName):
        sys.exit("Output report exists, delete it first and rerun: " + outReportName);

    ### Create out dirs
    create_out_dirs(param_outdir)
    
    ### Create snapshots
    create_snapshots()
    
    ### Create report
    create_html_report()
    
    
    
