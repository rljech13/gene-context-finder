/* A static report has no dangling portable-app download and starts no tools. */
const fs=require('fs'),assert=require('node:assert/strict'),{webcrypto}=require('node:crypto'),{JSDOM}=require('jsdom');
const code=fs.readFileSync('genecontext/review_realign.txt','utf8');
async function check(portable){
 const seed={stage:3,review_id:'static-control',portable_available:portable,selected:['a','b'],visible:['a','b'],sequences:[{id:'a',label:'A',group:'g',protein:'MA',cds:'ATGGCTTAA'},{id:'b',label:'B',group:'g',protein:'MT',cds:'ATGACTTAA'}]};
 const dom=new JSDOM('<section id="interactiveReview"></section><script type="application/json" id="interactiveReviewData">'+JSON.stringify(seed)+'</script>',{url:'file:///software-control/report.html',runScripts:'outside-only',beforeParse(w){Object.defineProperty(w,'crypto',{value:webcrypto});w.fetch=()=>{throw Error('Static viewing must not request an API or download')}}});
 dom.window.eval(code);await new Promise(resolve=>setTimeout(resolve,10));
 const document=dom.window.document;assert.equal(document.querySelector('a[href="review_app.zip"]').hidden,!portable);
 assert.match(document.getElementById('realignStatus').textContent,portable?/Start_review/:/Saved comparisons are available offline/);
 assert.match(document.body.textContent,/Export new choices and results before closing/);
 dom.window.close();
}
(async()=>{await check(false);await check(true);console.log(JSON.stringify({status:'passed',static_download_hidden:true,portable_download_retained:true,no_calculation_on_open:true}))})().catch(error=>{console.error(error);process.exit(1)});
