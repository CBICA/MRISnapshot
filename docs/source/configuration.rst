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

Users can input a mask image (optional). A mask is primarily used in
slice selection, but it can also be used for cropping the selected 
snapshots (optional).

Image column names can be modified. However, the suggested use is to keep the 
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

Slice selection is done independently in any combination of 3 *view_planes* 
(axial, sagittal and coronal)

Slices are selected based on "foreground voxels", determined primarily from the 
mask image, otherwise (if no mask image) from the overlay image, and otherwise 
(if no mask and overlay images) from the non-zero voxels of the underlay image

A slice can not be selected if it has less than *min_vox* foreground pixels 
(This parameter was introduced to allow discarding small structures 
during the QC of segmentations, e.g. to eliminate small lesions)  

From the slices with enough foreground pixels, a subset is selected using 
either the *num_slice* parameter (fixed number of slices with a regular step 
size) or the *step_size_slice* (fixed step size with a variable number of 
slices)

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
