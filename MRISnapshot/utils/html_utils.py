#!/usr/bin/env python

######################################################################
### prepare html code for single snapshot
def htmlMainPage(link_html_name):

    TEXT_HTML_STR = '''<!DOCTYPE html> <html> <body onLoad="triggerJS();"> 
<script> 
 function triggerJS(){ location.replace("''' +  link_html_name + '''"); }
</script>
</body> </html>'''

    return TEXT_HTML_STR

######################################################################
### prepare html code for single snapshot
def htmlSnapshot(img_name, fig_caption, onclick_html_name, imgw):

    TEXT_HTML_STR = '''<a href="''' + onclick_html_name + '''">
<img src="''' + img_name + '''" style="width:''' + str(imgw) + '''">
<figcaption>''' + fig_caption + '''</figcaption></a>
<center><button id="backBut">BACK</button></center>
<script src="../scripts/load_back.js" type="text/javascript" ></script>
'''

    return TEXT_HTML_STR


######################################################################
### prepare html code for subject, page header
def htmlSubjectPrefix(sub_id, sub_no, sub_count, fname_under, fname_over, fname_over2):
    
    TEXT_HTML_STR = '''<html><head><script src="scripts/save_qcform.js" type="text/javascript" ></script>
<script src="scripts/misc_func.js" type="text/javascript" ></script>
<script src="scripts/shortcut.js" type="text/javascript" ></script>
<link rel="stylesheet" href="scripts/pagestyle.css">
</head><body>
<hr><br>
<center><h3>SUBJECT: ''' + str(sub_no) + ' / ' + str(sub_count) + ''', ID: <label id="SUBJ_ID" name="SUBJ_ID">
''' + sub_id + '''</label></h3></center><br>
<h3>View cmd:</h3> 
fsleyes ''' + fname_under + ''' ''' + fname_over +  ''' ''' + fname_over2 + '''<br>
'''

    return TEXT_HTML_STR

######################################################################
### prepare html code for subject, add snapshot
def htmlSubjectAddImage(img_name, fig_caption, onclick_html_name):

    TEXT_HTML_STR = '''<div class="column">
  <div class="img"><a href="''' + onclick_html_name + '''"><img src="''' + img_name + '''"/></a></div>
  <p>''' + fig_caption + '''</p>
  </div> 
'''
    return TEXT_HTML_STR

######################################################################
### prepare html code that contains qc form
def htmlQCForm(txtCHECK1, txtCHECK2, txtTXT1):

    TEXT_HTML_STR = '''<center>
QC: <input type="checkbox" id="CHECK1" name="CHECK1" value="CHECK1"> ''' + txtCHECK1 + ''' &nbsp
<input type="checkbox" id="CHECK2" name="CHECK2" value="CHECK2"> ''' + txtCHECK2 + ''' &nbsp &nbsp &nbsp &nbsp &nbsp
''' + txtTXT1 + ''' <input type="text" id="TXT1" name="TXT1" value="" class="field left"> &nbsp &nbsp &nbsp &nbsp &nbsp
SAVE AS: <input id="inputFileNameToSaveAs"></input> &nbsp
<button onclick="saveQcAndReport()">SAVE REPORT</button>
<br>
</center>'''

    return TEXT_HTML_STR

######################################################################
def html_Navig(prevLink, nextLink):

    TEXT_HTML_STR = '''<center>
<a href="''' + prevLink + '''"><button title="Shortcut:Left Arrow">PREV</button></a>
<a href="''' + nextLink + '''"><button title="Shortcut:Right Arrow">NEXT</button></a>
<br>
</center>

<script>
shortcut.add("Left",function() {window.location.href="''' + prevLink + ''' "}); 
shortcut.add("Right",function(){window.location.href="''' + nextLink + ''' "}); 
</script>

'''    

    return TEXT_HTML_STR

#######################################################################
def html_stylesheet(imgWidth):

    TEXT_HTML_STR = '''.container{
 margin:10px;
 padding:10px;
}

.column{
 width:''' + str(imgWidth) + '''px;
 height:''' + str(imgWidth) + '''px;
 min-width:100px;
 padding:5px;    
 margin:5px;
 
 background:#fff;
 font-size:14px;
 
 display:inline-block;
 position:relative;
 float:center;
}
 
.column img{
 width:100%;
 height:100%;
 margin:auto;
 background:#ccc;
 display:block;
}

'''

    return TEXT_HTML_STR


