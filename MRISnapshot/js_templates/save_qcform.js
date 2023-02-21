// Write qc results to report
function saveQcAndReport() {
    saveQc();
    saveReportToFile();
}

function saveSubject() {
    var Check1 = document.getElementsByName('CHECK1');
    var Check2 = document.getElementsByName('CHECK2');
    var subjidText = document.getElementsByName('SUBJ_ID');
    var notesText = document.getElementsByName('TXT1');
    
    textToWrite = subjidText[0].textContent;
    if (Check1[i].checked)  {
        textToWrite = textToWrite + ',1';
    }
    else {
        textToWrite = textToWrite + ',0';
    }
    if (Check2[0].checked)  {
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
    var Check1 = document.getElementsByName('CHECK1');
    var Check2 = document.getElementsByName('CHECK2');
    var subjidText = document.getElementsByName('SUBJ_ID');
    var notesText = document.getElementsByName('TXT1');
    
    for (i = 0; i < subjidText.length; i++) {
        textToWrite = subjidText[i].textContent;
        if (Check1[i].checked)  {
                textToWrite = textToWrite + ',1';
        }
        else {
                textToWrite = textToWrite + ',0';
        }
        if (Check2[i].checked)  {
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
    var Check1 = document.getElementsByName('CHECK1');
    var Check2 = document.getElementsByName('CHECK2');
    var subjidText = document.getElementsByName('SUBJ_ID');
    var notesText = document.getElementsByName('TXT1');
    var textTmp;
    var tmpArray;
    
    for (i = 0; i < subjidText.length; i++) {
        textTmp = sessionStorage.getItem("id_" + subjidText[i].textContent) || '';
        tmpArray = textTmp.split(",");
        if (tmpArray.length > 1) {
            Check1[i].checked = parseInt(tmpArray[1]);
            Check2[i].checked = parseInt(tmpArray[2]);
            notesText[i].value = tmpArray[3];
        }
    }
}

// Write saved results to file
function saveReportToFile() {
    var textToWrite = window.textHdr;
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
