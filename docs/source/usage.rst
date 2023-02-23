*****
Usage
*****

Usage of MRISnapshot for a specific dataset and QC task is very similar to the simple example in :ref:`ref_quickstart`. This process is done in 3 steps:

**1. Data preparation**

Data preparation involves creation of two files in an empty directory that will be used as the output folder for the final report.

* **list_images.csv**: a list of images

* **config.csv**: a list of configuration parameters and their values

.. warning::
    The names of the two data files are constant and should not be modified.

The helper script: 

.. code-block:: console

    mrisnap_prep_data [-h] -i INDIR -s S_ULAY -d OUTDIR [--mask S_MASK] [--olay S_OLAY] [--olay2 S_OLAY2]

is provided for data preparation.

    **1.A. Structured data (input images with consistent and unified path and naming format)**

    The target images for QC are all in the same parent folder (with any level of nesting), and they follow a consistent naming format, specifically *{ScanID}{ScanSuffix}*, with unique and consistent suffixes for different types of input images
    
    Data preparation can be directly done using the helper script without additional work. The example in :ref:`ref_quickstart` illustrates this case.

    **1.B. Unstructured data (input images with non-unified paths and/or naming formats)**

    Target images are located in different directories without a common parent folder, and/or image names are not consistent. A typical example is initial raw scans extracted from dicoms.
    
    Users can still run the helper script for initialization. However, they will need to edit the image list, or create a new  image list with the full path and name of all input images.

    .. warning::
        Image list file should always include a first column named as "ID". This field is used as an identifier for the participant or the scan session.
        
        Other columns in the image list are "UnderlayImg", "MaskImg", "OverlayImg" and "OverlayImg2". We suggest to keep these column names. However, this is not a requirement. If a column name is modified, users should edit the related field in the configuration file.
    
**2. Selection of QC report configuration parameters**

User parameters for configuring the QC report are listed in the **config.csv** file. The helper script creates a configuration file with default values. Users can edit the values in this file to customize the QC report creation for the specific task. See :ref:`ref_config` for a detailed list of all parameters.

**3. Creation of the QC report**

Run the report creation script with the output folder created in previous steps as the single input argument.

.. code-block:: console

    mrisnap_create_report [-h] -d OUTDIR
