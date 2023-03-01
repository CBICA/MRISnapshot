*****
Usage
*****

Usage of MRISnapshot for a specific dataset and QC task is very similar to the simple example in :ref:`ref_quickstart`. This process is done in 3 steps:

**1. Data preparation**

    Data preparation involves creation of two files in an empty directory that will be used as the output folder for the final report.

    * **list_images.csv**: a list of images

    * **config.csv**: a list of configuration parameters and their values

    .. warning::
        The names of these two files are constant and should not be modified.

    The helper script: 

    .. code-block:: console

        mrisnap_prep_data [-h] -i INDIR -s S_ULAY -d OUTDIR [--mask S_MASK] [--olay S_OLAY] [--olay2 S_OLAY2]

    is provided for data preparation.

**1.A. Structured data (input images with consistent and unified path and naming format)**

    The target images for QC are all in the same parent folder (with any level of nesting), and they follow a consistent naming format, specifically *{ScanID}{ScanSuffix}*, with unique and consistent suffixes for different types of input images
    
    Data preparation can be directly done using the helper script without additional work. The example in :ref:`ref_quickstart` illustrates this case.
    
    .. warning::
        The script uses the underlay image names detected in the input folder to create the initial list of scan id's. The file name after removing the suffix will be considered as the scan id for an image, and images with the same scan id will be removed from the list (e.g. two files with the same name in two different sub-folders).
        
        Overlay and mask images will be added to the list **only if they have the same exact scan id as the underlay image**. Users should enter the complete suffix for all image types to make sure that the matching by scan id will work as intended. 
        

**1.B. Unstructured data (input images with non-unified paths and/or naming formats)**

    Target images are located in different directories without a common parent folder, and/or image names are not consistent. A typical example is the QC of initial raw scans extracted from dicoms.
    
    Users can still run the helper script for initialization. However, they will need to edit the image list, or create a new  image list with the full path and name of all input images.

    .. warning::
        Image list file should include a first column named as *"ScanID"*. This field is used as an identifier for the participant or the scan session.
        
    .. note::
        Other columns in the image list are *"UnderlayImg"*, *"MaskImg"*, *"OverlayImg"* and *"OverlayImg2"*. Users can use names other than these default ones (e.g. *T1* instead of *ulay_col*). If they do so, they should edit the configuration file to make sure that the column name matches the related field in the configuration file.
    
**2. Selection of QC report configuration parameters**

    User parameters for configuring the QC report are listed in the **config.csv** file. The helper script creates a configuration file with default values. Users can edit the values in this file to customize the QC report creation for the specific task. See :ref:`ref_config` for explanations, and detailed lists of all parameters.

**3. Creation of the QC report**

    Run the report creation script, providing as input the name of the folder with data files created in previous steps. The QC report will be created in the same folder.

    .. code-block:: console

        mrisnap_create_report [-h] -d OUTDIR

    The QC report can be viewed using the command:

    .. code-block:: console

        google-chrome {output_directory}/QCReport/qcreport.html

**4. Logs for missing or problematic images**

    QC report creation process aims to be robust. So, problematic cases will not interrupt report creation. If a case fails, it will be skipped, and this will be reported in a log file with a brief description of the problem, to allow users additional data 
    verifications. Possible reasons for failure are:

    * Image missing (underlay, mask or overlay)

    * Image could not be read

    * Overlay or mask image that is used for slice selection has no foreground voxels

    A log file that lists the QC status of all subjects will be saved in the output folder in :

        .. code-block:: console

            QCReport/log_qc_images_all.csv

    A log file that lists the QC status of failed subjects will be saved in the output folder in :

        .. code-block:: console

            QCReport/log_qc_images_fail.csv



    
