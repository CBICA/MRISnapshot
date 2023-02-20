Description:
------------
Configuration template to create a QC report for the brain mask

Details:
--------
- Overlays the edges of the brain mask on the T1 scan
- Shows snapshots for 5 slices in Axial and Sagittal views
- Slice selection is done using the overlay image as a mask
- Underlay image contrast is scaled between 2 and 98 percentiles
- Snapshot size on the QC report is 300x300 pixels

Commands to run:
----------------
## Create img list
mrisnapshot_prep_data -i ../../Scans -s _T1_DS.nii.gz -d . --olay _T1_BRAINMASK.nii.gz
## Create report
mrisnapshot_create_report -d .
## View report
google-chrome QCReport/qcreport.html 

Usage of the report:
--------------------
- Left / write arrows: Navigate between subjects (or use the PREV / NEXT buttons)
- Click on a snapshot: Open it on a larger image (middle click to open on a new tab)
- Click on the larger image: Hide / show the overlay
- SAVE REPORT button: Save the report on file (e.g. qcreport.csv)

