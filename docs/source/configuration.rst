.. _ref_config:

*************
Configuration
*************

All user parameters for QC report creation are collected in a single csv file, **config.csv**. Initial configuration file with default values can be created using the  *mrisnapshot_prep_data* helper function, or can be copied from a configuration template.

.. warning::
    **config.csv** consists of two columns, *ParamName* and *ParamValue*. Values in *ParamName* column should not be edited. Values in *ParamValue* column can be edited to customize the QC report.

Configuration file includes groups of variables that modify different aspects of QC report creation:

**1. Image selection parameters**

These parameters control the images that will be displayed in each snapshot.

.. csv-table::
   :file: ../resources/config_file_descr_images.csv
   :widths: 50, 50, 30, 60, 80
   :header-rows: 1
   
**2. Slice selection parameters**

These parameters control the slices selected for creating the snapshots.

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
