"""Hydrate saved checkpoints before executing classic report scripts.

Keep separate script tags: native review scripts share global lexical bindings.
"""
from pathlib import Path
import re

SCRIPT = re.compile(r'(<script\b[^>]*>)(.*?)(</script\s*>)', re.I | re.S)


def finalize(html):
    if 'id="interactiveReviewData"' not in html:
        return html
    html = re.sub(r'<script id="gcfReviewStorageLoader">.*?</script>', '', html, flags=re.S)
    if '// gcf-review-storage-v2' not in html:
        runtime = Path(__file__).with_name('review_storage.txt').read_text()
        html = html.replace('</head>', '<script>'+runtime+'</script></head>', 1)
    head, sep, body = html.partition('</head>')
    if not sep or '</body>' not in body:
        raise ValueError('Expected a complete interactive report document')

    def defer(match):
        start, source, end = match.groups()
        if re.search(r'\b(?:type|src)\s*=', start, re.I):
            return match[0]
        return start[:-1]+' type="application/x-gcf-review-script">'+source+end

    body = SCRIPT.sub(defer, body)
    loader = """<script id="gcfReviewStorageLoader">
(async()=>{
 await window.GCFStorageReady;
 await window.GCFPortableReady;
 for(const old of document.querySelectorAll('script[type="application/x-gcf-review-script"]')){
  const script=document.createElement('script');
  for(const a of old.attributes)if(a.name!=='type')script.setAttribute(a.name,a.value);
  script.textContent=old.textContent;old.replaceWith(script);
 }
 document.documentElement.setAttribute('data-gcf-review-ready','true');
 document.dispatchEvent(new Event('gcf-review-ready'));
})().catch(error=>{
 const warning=document.createElement('p');warning.setAttribute('role','alert');
 warning.textContent='Could not restore the saved review: '+error.message;
 document.querySelector('main')?.prepend(warning);
});
</script>"""
    return head+sep+body.replace('</body>', loader+'</body>', 1)
