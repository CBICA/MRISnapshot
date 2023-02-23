.. _ref_config:

*************
Configuration
*************

All user parameters for QC report creation are collected in a single csv file, **config.csv**. Initial configuration file with default values can be created using the  *mrisnapshot_prep_data* helper function, or can be copied from a configuration template.

.. warning::
    **config.csv** consists of two columns, *ParamName* and *ParamValue*. Values in *ParamName* column should not be edited. Values in *ParamValue* column can be edited to customize the QC report.

Configuration file includes groups of variables that modify different aspects of QC report creation:

**1. Image list parameters**

These parameters control the selection of input images. The report will display 
snapshots from either a single underlay image (required), or with one or two 
overlay masks (optional).

Users can input a mask image (optional). A mask is primarily used for
slice selection, but it can also be used for cropping the selected 
snapshots (optional).

Image column names can be modified (e.g. *T1* instead of *ulay_col*). However, the suggested use is to keep the 
default values.

.. csv-table::
   :file: ../resources/config_file_descr_images.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1
   
**2. Slice selection parameters**

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
    
    With the example dataset, this can be tested using the _ROIMASK.nii.gz images as overlay, and setting "sel_vals_olay = 47+48", the indices for the left hippocampus (47) and right hippocampus (48))

If a slice has less than *min_vox* foreground pixels, it will not be included in the foreground mask for the specific view plane.

.. note::
    This parameter was introduced to allow discarding small structures during the QC of segmentations. For example, to discard slices with small lesions in the QC of WM lesion segmentation

From the *foreground mask*, a subset of slices are selected using either the *num_slice* parameter (fixed number of slices with a regular step size) or the *step_size_slice* (fixed step size with a variable number of slices)

.. note::
    When the size of the overlay or mask image is predictable and less variable (e.g. brain mask) it may be suitable to set the *num_slice* parameter to display a fixed number of equally spaced slices.
    
    For segmentation labels with variable size (e.g WM lesions) users may prefer to set the *step_size_slice* parameter to display more or less slices depending on the segmented volumes
    

Selected slices can be displayed in full, cropped to mask 
(bounding box of the mask foreground), or cropped to overlay (bounding box of 
the overlay foreground). Users can extend the cropping area by setting the 
*padding_ratio*. For example *padding_ratio = 0.1* will add a padding of 
%10 times size of the bounding box.

.. csv-table::
   :file: ../resources/config_file_descr_sliceselection.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1

**3. Visual parameters**

These parameters control various visual aspects of the snapshots.

.. csv-table::
   :file: ../resources/config_file_descr_view.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1
   
   
**4. Report parameters**

These parameters control various aspects of the QC report.

.. csv-table::
   :file: ../resources/config_file_descr_report.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1
