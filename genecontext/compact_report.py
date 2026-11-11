"""Optional standalone gzip transport; scientific evidence and renderer stay intact."""
import base64
import gzip
from html import escape
from pathlib import Path
import re
from urllib.parse import quote


def compact_html(text, fallback='report.html'):
    scripts = list(re.finditer(r'<script>\s*(.*?)</script>', text, flags=re.DOTALL))
    if len(scripts) != 1:
        raise ValueError('Expected one inline application script in a Gene Context report')
    script = scripts[0]
    encoded = base64.b64encode(gzip.compress(script.group(1).encode(), mtime=0)).decode()
    loader = '''<div id="gcf-compact-loading" role="status">Loading report… <a href="FALLBACK">Open the original report</a></div>
<script>
(async function(){
  const message=document.getElementById('gcf-compact-loading');
  try {
    if(typeof DecompressionStream==='undefined') throw new Error('This browser cannot decompress this report.');
    const bytes=Uint8Array.from(atob('PAYLOAD'), c=>c.charCodeAt(0));
    const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
    const source=await new Response(stream).text();
    const application=document.createElement('script');
    application.textContent=source;
    document.body.appendChild(application);
    message.remove();
  } catch(error) {
    message.prepend(document.createTextNode(String(error)+' '));
  }
})();
</script>'''.replace('FALLBACK', escape(fallback, quote=True)).replace('PAYLOAD', encoded)
    return text[:script.start()] + loader + text[script.end():]


def run(args):
    source = Path(args.report).resolve(); target = Path(args.out).resolve()
    if target == source or target.parent != source.parent:
        raise ValueError('Choose a new sibling HTML file, keeping relative sequence/export links valid')
    text = compact_html(source.read_text(), quote(source.name))
    with target.open('x') as handle:
        handle.write(text)
    print(f'Compact report: {target} ({target.stat().st_size:,} bytes; original {source.stat().st_size:,})')
    return 0
