
// Save page URL, scroll position and qc results before leaving
window.onbeforeunload = function() {
	var subjidText = document.getElementById('SUBJ_ID').textContent;
	sessionStorage.setItem("global_scrollPosition", window.pageYOffset);
	sessionStorage.setItem("currentPage",document.URL);
	saveQc();
};

// window.onbeforeunload = function() {
// 	var subjidText = document.getElementById('SUBJ_ID');
// 	sessionStorage.setItem(subjidText.textContent + "_scrollPosition", window.pageYOffset);
// 	sessionStorage.setItem("currentPage",document.URL);
// 	saveQc();
// };


// Set page scroll position and qc results
window.onload = function () {
	window.scrollBy(0,sessionStorage.getItem("global_scrollPosition"));
	loadQc();
}

// window.onload = function () { 
// 	var subjidText = document.getElementById('SUBJ_ID');
// 	window.scrollBy(0,sessionStorage.getItem(subjidText.textContent + "_scrollPosition"));
// 	loadQc();
// }


