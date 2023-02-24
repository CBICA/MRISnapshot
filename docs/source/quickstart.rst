.. _ref_quickstart:

**********
Quickstart
**********

After installation, users can quickly run a small example:

1. Go to examples script directory.

.. code-block:: console

    cd examples/Scripts

2. Run the command for data preparation.

.. code-block:: console

    mrisnapshot_prep_data -i ../Input/Scans -s _T1_DS.nii.gz -d ../Output/Test --mask _T1_ICVMASK.nii.gz --olay _T1_ICVMASK.nii.gz --olay2 _T1_BRAINMASK.nii.gz

3. Run the command for report creation.

.. code-block:: console

    mrisnapshot_create_report -d ../Output/Test

.. note::
    Alternatively, instead of steps 2 and 3, users can run the example script provided.
    
    .. code-block:: console

        ./run_example.sh

4. View the QC report (using your favorite browser).

.. code-block:: console

    google-chrome ../Output/Test/QCReport/qcreport.html

5. Navigate subjects using **PREV** and **NEXT** buttons (or left and right arrows), and annotate them by editing the QC form fields.

.. |br| raw:: html

    <br />
    
6. Save the QC annotations using the **SAVE REPORT** button.

.. note::
    Depending on the browser type and configuration, the path and naming of the output *.csv* file may show differences. For example, the output file may be directly saved to a default download folder.



