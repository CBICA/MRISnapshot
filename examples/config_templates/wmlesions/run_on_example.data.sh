## Commands to use this template with example scans

## Create img list
mrisnapshot_prep_data -i ../../scans -s _T1.nii.gz -d . --olay _T1_WMLMASK.nii.gz

## Create report
mrisnapshot_create_report -d .

## View report
google-chrome ./QCReport/qcreport.html
