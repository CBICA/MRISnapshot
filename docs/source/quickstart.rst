**********
Quickstart
**********

The typical use case for ``MRISnaphot`` is to assist QC of a dataset with raw and processed MRI images. Manual QC, i.e. by viewing 3D images and/or overlaid derived masks, may be time consuming for large datasets. ``MRISnaphot`` allows users to create snapshots of selected underlay and overlay images in a configurable way.

QC report creation involves 2 steps:

* Creating an input list of images and configuration file with parameters

* Applying the report creation script

.. code-block:: console

    (.venv) $ mrisnapshot_prep_data
    (.venv) $ mrisnapshot_create_report
