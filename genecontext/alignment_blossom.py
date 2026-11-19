"""Presentation-only GeCoS Blossom shading for the protein MSA viewers.

Exact 20-residue palette from Jalview's Blossom documentation:
https://www.jalview.org/help/html/colourSchemes/gecos-blossom.html
Kunzmann et al. (2020), doi:10.1186/s12859-020-3526-6.
Unknown residues, stops and gaps have no invented Blossom colour.
"""
BLOSSOM = {
    'A': '#8bc4b4', 'C': '#0893fe', 'D': '#5fa505', 'E': '#dbb501',
    'F': '#f74fa8', 'G': '#00d382', 'H': '#ff5701', 'I': '#9abaf3',
    'K': '#fea527', 'L': '#cda5dc', 'M': '#f5a1b8', 'N': '#b5c206',
    'P': '#10d631', 'Q': '#bf8527', 'R': '#fc9502', 'S': '#7e9d59',
    'T': '#00a29c', 'V': '#87c0e4', 'W': '#fe08fb', 'Y': '#ff4e7a',
}

MARKER = '<!-- gcf-gecos-blossom-v2 -->'


def colour_alignment_html(html):
    """Add offline shading without touching payloads, decisions or alignment text.

    The observer reapplies shading after switching sequences, pairs or groups.
    CDS keeps the viewer's existing nucleotide/consensus shading.
    """
    if MARKER in html:
        return html
    if 'id="alignment"' not in html or 'id="alignmentKind"' not in html:
        return html
    css = '#alignment.gecos-protein .res{background:#fff;color:#111}'
    css += ''.join('#alignment.gecos-protein .res[data-aa="'+aa+'"]{background:'+colour+'}'
                   for aa, colour in BLOSSOM.items())
    css += '#alignment.gecos-protein .res.gap{background:transparent;color:#bdc8c4}'
    script = r"""
(()=>{
const box=document.getElementById('alignment'),kind=document.getElementById('alignmentKind');
if(!box||!kind)return;
const paint=()=>{
 const protein=kind.value==='protein';
 box.classList.toggle('gecos-protein',protein);
 box.dataset.colourScheme=protein?'gecos-blossom':'original-cds';
 for(const residue of box.querySelectorAll('.res')){
  if(protein)residue.dataset.aa=residue.textContent.toUpperCase();
  else delete residue.dataset.aa;
 }
 const legend=document.getElementById('alignmentLegend');
 if(legend&&protein){
  const text=legend.textContent,at=text.indexOf('Dashes:');
  legend.textContent='Protein: GeCoS Blossom (Jalview); colour identifies the amino acid, not identity or confidence. '+(at<0?'Gaps and unknown residues are uncoloured.':text.slice(at));
 }
 const technical=box.closest('section')?.querySelector('.legend');
 if(technical)technical.textContent=protein?'Protein: GeCoS Blossom (Jalview); colour identifies the amino acid, not identity or confidence. Gaps and unknown residues are uncoloured. Position numbers count ungapped residues.':'CDS: nucleotide colours identify A/C/G/T. Dashes are alignment gaps. Position numbers count ungapped bases.';
};
new MutationObserver(paint).observe(box,{childList:true,subtree:true});
kind.addEventListener('change',paint);paint();
})();
"""
    fragment = MARKER+'\n<style>'+css+'</style>\n<script>'+script+'</script>\n'
    if '</body>' not in html:
        raise ValueError('Alignment report has no closing body')
    return html.replace('</body>', fragment+'</body>', 1)
