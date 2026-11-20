/* Shared integration coverage for the current inline model comparison. */
const check=require('./stage2_inline_models_dom.cjs');module.exports=check;if(require.main===module)console.log(JSON.stringify(check(process.argv[2])));
