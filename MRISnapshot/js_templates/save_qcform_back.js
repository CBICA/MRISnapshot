// Write qc results to report
function saveQcAndReport() {
	saveQc();
	saveReportToFile();
}

function saveSubject() {
    var qcpassCheck = document.getElementsByName('QC_PASS');
    var qcfailCheck = document.getElementsByName('QC_FAIL');
    var subjidText = document.getElementsByName('SUBJ_ID');
    var notesText = document.getElementsByName('QC_NOTES');
    
    textToWrite = subjidText[0].textContent;
    if (qcpassCheck[i].checked)  {
        textToWrite = textToWrite + ',1';
    }
    else {
        textToWrite = textToWrite + ',0';
    }
    if (qcfailCheck[0].checked)  {
        textToWrite = textToWrite + ',1';
    }
    else {
        textToWrite = textToWrite + ',0';
    }
    textToWrite = textToWrite + ',' + notesText[0].value;

    var textFileAsBlob = new Blob([textToWrite], {type:'text/plain'});
    var fileNameToSaveAs = 'mriqcreport_' + subjidText[0].textContent + '.csv';

    var downloadLink = document.createElement("a");
    downloadLink.download = fileNameToSaveAs;
    downloadLink.innerHTML = "Download File";
    if (window.webkitURL != null)
    {
            // Chrome allows the link to be clicked
            // without actually adding it to the DOM.
            downloadLink.href = window.webkitURL.createObjectURL(textFileAsBlob);
    }
    else
    {
            // Firefox requires the link to be added to the DOM
            // before it can be clicked.
            downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
            downloadLink.onclick = destroyClickedElement;
            downloadLink.style.display = "none";
            document.body.appendChild(downloadLink);
    }

    downloadLink.click();
    
    
    
}

// Write qc results to session var
function saveQc() {
	var qcpassCheck = document.getElementsByName('QC_PASS');
	var qcfailCheck = document.getElementsByName('QC_FAIL');
	var subjidText = document.getElementsByName('SUBJ_ID');
	var notesText = document.getElementsByName('QC_NOTES');
	
	for (i = 0; i < subjidText.length; i++) {
		textToWrite = subjidText[i].textContent;
		if (qcpassCheck[i].checked)  {
				textToWrite = textToWrite + ',1';
		}
		else {
				textToWrite = textToWrite + ',0';
		}
		if (qcfailCheck[i].checked)  {
				textToWrite = textToWrite + ',1';
		}
		else {
				textToWrite = textToWrite + ',0';
		}
		textToWrite = textToWrite + ',' + notesText[i].value;
		
		sessionStorage.setItem("id_" + subjidText[i].textContent,textToWrite);
	}
}

// Load qc results from session var
function loadQc() {
	var qcpassCheck = document.getElementsByName('QC_PASS');
	var qcfailCheck = document.getElementsByName('QC_FAIL');
	var subjidText = document.getElementsByName('SUBJ_ID');
	var notesText = document.getElementsByName('QC_NOTES');
	var textTmp;
	var tmpArray;
	
	for (i = 0; i < subjidText.length; i++) {
		textTmp = sessionStorage.getItem("id_" + subjidText[i].textContent) || '';
		tmpArray = textTmp.split(",");
		if (tmpArray.length > 1) {
			qcpassCheck[i].checked = parseInt(tmpArray[1]);
			qcfailCheck[i].checked = parseInt(tmpArray[2]);
			notesText[i].value = tmpArray[3];
		}
	}
}

// Write saved results to file
function saveReportToFile() {
	var textToWrite = 'ID,QC_PASS,QC_FAIL,NOTES';
	for(var i in sessionStorage) {
		if (i.substr(0,3)=="id_") {
    			textToWrite = textToWrite + sessionStorage[i];
    		}
	}
	
	var textFileAsBlob = new Blob([textToWrite], {type:'text/plain'});
	var fileNameToSaveAs = document.getElementById("inputFileNameToSaveAs").value;

	var downloadLink = document.createElement("a");
	downloadLink.download = fileNameToSaveAs;
	downloadLink.innerHTML = "Download File";
	if (window.webkitURL != null)
	{
		// Chrome allows the link to be clicked
		// without actually adding it to the DOM.
		downloadLink.href = window.webkitURL.createObjectURL(textFileAsBlob);
	}
	else
	{
		// Firefox requires the link to be added to the DOM
		// before it can be clicked.
		downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
		downloadLink.onclick = destroyClickedElement;
		downloadLink.style.display = "none";
		document.body.appendChild(downloadLink);
	}

	downloadLink.click();
}

function destroyClickedElement(event)
{
	document.body.removeChild(event.target);
}
