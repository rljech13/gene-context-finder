"""Smaller transport preserves the application, data, DOM markup and exports."""
import base64
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest

from genecontext.compact_report import compact_html, run

NODE=shutil.which('node') or os.environ.get('GENE_CONTEXT_TEST_NODE')


class CompactReport(unittest.TestCase):
    def fixture(self):
        data=dict(review_id='same',payload=['NANOG & <evidence>']*1000)
        application='const DATA='+json.dumps(data).replace('<','\\u003c')+'; globalThis.answer=DATA;'
        html='<html><body><main>Original markup</main><script>'+application+'</script></body></html>'
        return data,application,html

    def test_application_roundtrip_and_surrounding_markup_are_exact(self):
        _,application,html=self.fixture();result=compact_html(html)
        encoded=re.search(r"atob\('([^']+)'\)",result)[1]
        self.assertEqual(gzip.decompress(base64.b64decode(encoded)).decode(),application)
        self.assertTrue(result.startswith('<html><body><main>Original markup</main>'))
        self.assertTrue(result.endswith('</script></body></html>'))
        self.assertLess(len(result),len(html)//3)

    def test_cli_retains_original_and_refuses_replacement_or_broken_relative_assets(self):
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/'report.html';target=Path(d)/'report.compact.html'
            source.write_text(self.fixture()[2]);before=source.read_bytes()
            self.assertEqual(run(SimpleNamespace(report=source,out=target)),0)
            self.assertEqual(source.read_bytes(),before)
            with self.assertRaises(FileExistsError):run(SimpleNamespace(report=source,out=target))
            with self.assertRaises(ValueError):run(SimpleNamespace(report=source,out=source))
            with self.assertRaises(ValueError):run(SimpleNamespace(report=source,out=source.parent/'elsewhere/report.html'))

    @unittest.skipUnless(NODE,'requires Node for actual native decompression/bootstrap')
    def test_actual_javascript_bootstrap_preserves_global_data_and_capability_fallback(self):
        data,_,html=self.fixture();result=compact_html(html)
        script=re.search(r'<script>(.*?)</script>',result,re.S)[1]
        harness=r'''
const vm=require('node:vm');
const input=JSON.parse(require('node:fs').readFileSync(0,'utf8'));
(async()=>{
 for(const supported of [true,false]){
  let removed=false,notice=false;
  const message={remove(){removed=true},prepend(){notice=true}};
  let context;
  const document={getElementById(){return message},createElement(){return {}},createTextNode(s){return s},
   body:{appendChild(script){vm.runInContext(script.textContent,context)}}};
  context=vm.createContext({DecompressionStream:supported?DecompressionStream:undefined,Blob,Response,Uint8Array,atob,document});
  await vm.runInContext(input.script,context);
  if(supported){
   if(JSON.stringify(context.answer)!==JSON.stringify(input.data)||!removed)throw Error('Changed global data');
  }else if(!notice||context.answer)throw Error('Missing capability fallback');
 }
})().catch(e=>{console.error(e);process.exitCode=1});
'''
        result=subprocess.run([NODE,'-e',harness],input=json.dumps(dict(script=script,data=data)),capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
