.. _ref_config:

*************
Configuration
*************

All user parameters for QC report creation are collected in a single csv file, **config.csv**. Initial configuration file with default values can be created using the  *mrisnapshot_prep_data* helper function, or can be copied from a configuration template.

.. warning::
    **config.csv** consists of two columns, *ParamName* and *ParamValue*. Values in *ParamName* column should not be edited. Values in *ParamValue* column can be edited to customize the QC report.

Configuration file includes groups of variables that modify different aspects of QC report creation:

**1. Image list parameters**

.. csv-table::
   :file: ../resources/config_file_descr_images.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1

These parameters control the selection of input images. The report will display 
snapshots from either a single underlay image (required), or with one or two 
overlay masks (optional).

Users can input a mask image (optional). A mask is primarily used for
slice selection, but it can also be used for cropping the selected 
snapshots (optional).

Image column names can be modified (e.g. *T1* instead of *ulay_col*). However, the suggested use is to keep the 
default values.
   
**2. Slice selection parameters**

.. csv-table::
   :file: ../resources/config_file_descr_sliceselection.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1

An important task in QC report creation is the slice selection and extraction 
of final snapshots. We use a set of parameters to control this task for various 
use cases, such as:

* QC of raw scans

* QC of segmentations focusing on the whole image

* QC of segmentations focusing only on the segmented structures

* QC of selected labels in segmentations

Slice selection is done independently in each view plane. Users can select any combinations of 3 view planes, axial, sagittal and coronal (*view_plane*).

For slice selection, we first create a *"foreground mask"*. This mask is primarily derived from the mask image; from the overlay image(s) if no mask exists; and from the non-zero voxels of the underlay image if no mask and overlay images exist

Users can select a subset of labels from overlay images (*sel_vals_olay*, *sel_vals_olay2*). By default, all labels from an overlay image are used for slice selection and display. If a user selects specific labels, all other values are reset to 0, and the updated overlay image is used for slice selection and display.

.. note::
    A typical use case scenario is the QC of specific regions of interest (ROI). For example, users may want to create QC snapshots that only show the hippocampus 
    
    With the example dataset, this can be tested using the _ROIMASK.nii.gz images as overlay, and setting "sel_vals_olay = 47+48", the indices for the left hippocampus (47) and right hippocampus (48)

If a slice has less than *min_vox* foreground pixels, it will not be included in the foreground mask for the specific view plane.

.. note::
    This parameter was introduced to allow discarding small structures during the QC of segmentations. For example, to discard slices with small lesions in the QC of WM lesion segmentation

From the *foreground mask*, a subset of slices are selected using either the *num_slice* parameter (fixed number of slices with a regular step size) or the *step_size_slice* (fixed step size with a variable number of slices)

.. note::
    When the size of the overlay or mask image is predictable and not variable (e.g. brain mask), it may be suitable to set the *num_slice* parameter, in order to display a fixed number of equally spaced slices.
    
    For segmentation labels with variable size (e.g WM lesions) users may prefer to set the *step_size_slice* parameter to display more or less slices depending on the segmented volumes
    
    For this second case, the typical usage is not to input a mask image, so that the slice selection is limited to segmented structures on the overlay image (e.g. show only slices with lesions). However, the users may also provide a mask to select slices independent of the overlay images (e.g. showing slices that are in the brain mask with or without lesions)    

Selected slices can be displayed in full size, cropped to the mask (bounding box of the foreground volume in the mask), or cropped to overlay images (bounding box of the foreground volume in the overlay images).

Users can also extend the cropping area by setting the *padding_ratio*. For example *padding_ratio = 0.1* will increase the padding area by %10 of its initial size


**3. Visualization parameters**

.. csv-table::
   :file: ../resources/config_file_descr_view.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1

These parameters control few visual aspects of the display. Users can binarize the overlay images (*bin_olay*); extract the edges of the overlays (*is_edge*); change the transparency of the overlays (*alpha_olay*); and change the contrast of the underlay image by mapping the minimum and maximum image intensities to selected intensity percentile values calculated from the complete image (*perc_low*, *perc_high*)

.. note:
    If an overlay image has multiple labels, extraction of edges is done independently for each label

**4. Report parameters**

.. csv-table::
   :file: ../resources/config_file_descr_report.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1

These parameters control few aspects of the final QC report. Users can create a single *.html* report for the complete set, instead of an independent *.html* file for each subject (*is_out_single*); prefer not to include the QC form in the report (*is_out_noqc*); and set the size of the snapshots in the report (*img_width*)

.. note::
    Default QC report is a collection of independent *.html* files. This format was selected as default, because it provides a better separation of subjects and easier annotation using the quick navigation shortcuts (left / right arrows)
    
.. note::
    The style of the QC report is managed by a stylesheet that is saved in *QCReport/subjects/scripts/pagestyle.css*. Users can edit this file to modify the style of the QC report without rerunning the report creation (e.g. change the snapshot size by editing the *.column:width* value)
    
The QC form includes 2 check boxes and 1 edit box. The names of these objects are used as the columns of the final QC annotation *.csv* file. Users can rename these objects in a way that will better reflect their task (*label_checkbox1*, *label_checkbox2*, *label_editbox*). For example, for the QC of raw T1 scans users may prefer to rename them as "ImgHasMotion", "ImgCorrupt" and "Details"


