// Set back button and page scroll position
window.onload = function () {
   window.scrollBy(0,sessionStorage.getItem("page_scrollPosition"));
   document.getElementById("backBut").onclick = function () {
   document.location=sessionStorage.getItem("currentPage");};
}
