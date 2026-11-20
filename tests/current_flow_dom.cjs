/* Generate a synthetic native checkpoint and inspect the actual generic page. */
const fs=require('node:fs'),os=require('node:os'),path=require('node:path'),{execFileSync}=require('node:child_process');
const root=fs.mkdtempSync(path.join(os.tmpdir(),'gcf-current-flow-'));
try{
const code='from pathlib import Path\nimport sys\nfrom test_project_review import workspace\np,source,choices,review=workspace(Path(sys.argv[1]))\nprint(review["report"])';
const report=execFileSync('python3',['-c',code,root],{env:{...process.env,PYTHONPATH:'tests'+path.delimiter+'.'},encoding:'utf8'}).trim();
process.stdout.write(execFileSync(process.execPath,['tests/sequential_page_dom.cjs',report],{encoding:'utf8'}));
}finally{fs.rmSync(root,{recursive:true,force:true})}
