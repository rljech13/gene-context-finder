"""Traceable locus names, kept separate from protein similarity and review groups.

This layer annotates existing evidence. It neither assigns cross-species
orthology nor changes saved choices, sequences, physical groups or alignments.
"""
from copy import deepcopy
import re


# Exact species + source gene, not a sequence-similarity or coordinate rule.
# The historical Nanogl1 / modern LOC correspondence is a curated local label;
# it does not assert expression or a conserved pluripotency function.
CURATED_COPIES = {
    ('Gallus gallus', '101748883'): dict(
        label='NANOG-like', family_hint='NANOG',
        aliases=['LOC101748883', 'Nanogl1', 'Nanog-like 1'],
        sources=['https://www.ncbi.nlm.nih.gov/gene/101748883',
                 'https://doi.org/10.1002/dvdy.20927'],
        note='Curated chicken NANOG-like locus; its role in pluripotency is not established here.'),
}


def _species(record):
    genome_species = str(record.get('genome_id', '')).split('_GC')[0].replace('_', ' ')
    value=record.get('species', '')
    # Stage-1 genomic loci use a list of query-source species in this field;
    # it is not the target organism and cannot establish a curated copy name.
    return genome_species if '_GC' in str(record.get('genome_id', '')) else value if isinstance(value,str) else ''


def _gene_id(record, annotation):
    for value in [record.get('ncbi_gene_id'), record.get('annotation_gene_id'),
                  record.get('original_gene'), record.get('gene_id'), annotation]:
        match = re.fullmatch(r'(?:gene-)?(?:LOC|GeneID:)?(\d+)', str(value or ''))
        if match:
            return match[1]
    return ''


def _homology_names(record):
    info = record.get('identification', {})
    # proposed_name can be a protein accession when the database has no symbol;
    # an accession alone is not a gene-family name.
    names = [h.get('gene_symbol', '') for h in info.get('near_best', [])]
    # Preserve -like, numbered paralogs and taxon-specific spelling. Only case
    # duplicates are collapsed, as in the underlying protein identification.
    return list({str(n).casefold(): str(n) for n in names if n}.values())


def identity(record):
    """Return a conservative identity description for one annotated locus/model."""
    annotation = str(record.get('annotation_gene_name') or record.get('annotation_name')
                     or record.get('name') or record.get('gene_symbol') or '')
    gene_id = _gene_id(record, annotation)
    homology = _homology_names(record)
    annotated = record.get('source') != 'tiberius'
    curated = CURATED_COPIES.get((_species(record), gene_id)) if annotated else None
    result = dict(label='', display_name='', status='unresolved', basis='unresolved',
                  original_name=annotation, original_gene=record.get('annotation_gene_id')
                  or record.get('original_gene') or record.get('gene_id') or '',
                  ncbi_gene_id=gene_id if annotated else '', family_hint='',
                  homology_names=homology, sources=[], notes=[],
                  orthology_confirmed=False, function_status='not_assessed',
                  database_name_differs=False)
    if curated:
        result.update(label=curated['label'], display_name=curated['label']+' ('+(annotation or 'GeneID '+gene_id)+')',
                      status='curated_locus_label', basis='species_and_source_gene_id',
                      family_hint=curated['family_hint'], sources=list(curated['sources']),
                      aliases=list(curated['aliases']), notes=[curated['note']])
    elif annotated and annotation and not re.fullmatch(r'LOC\d+', annotation):
        result.update(label=annotation, display_name=annotation, status='annotation_label',
                      basis='provided_annotation')
    elif homology:
        hints = ' / '.join(homology)
        label = hints+' family candidate'
        result.update(label=label, display_name=label+(' ('+annotation+')' if annotated and annotation else ''),
                      status='homology_only', basis='protein_similarity_only',
                      notes=['Protein similarity does not identify a particular paralog.'])
    else:
        label = annotation if annotated and annotation else 'Unassigned candidate'
        result.update(label=label, display_name=label,
                      notes=['Specific copy identity is unresolved.'])
    if result['status'] in {'curated_locus_label', 'annotation_label'} and homology:
        result['database_name_differs'] = all(n.casefold() != result['label'].casefold() for n in homology)
        if result['database_name_differs']:
            result['notes'].append('Protein-database names differ from this locus label; they do not replace it.')
    return result


def model_identities(data, metadata=None):
    """Use an annotated locus for predictions only with an exact CDS structure.

    Shared best hits, an overlapping interval or a broad physical-group merge
    cannot transfer a copy name. Native contig, strand, full phased CDS/stop
    coordinates and an available matching CDS fingerprint are required.
    """
    models = data.get('models', [])
    results = {m['model_id']: identity(m) for m in models}
    def signature(m):
        return (m.get('genome_id'), m.get('contig'), m.get('strand'),
                tuple(tuple(b) for b in m.get('cds_blocks', [])),
                tuple(tuple(b) for b in m.get('stop_blocks', [])), m.get('cds_sha256'))
    annotations = {}
    for m in models:
        result = results[m['model_id']]
        if m.get('source') == 'annotation' and m.get('cds_sha256') and result['status'] in {'curated_locus_label', 'annotation_label'}:
            annotations.setdefault(signature(m), []).append(result)
    for m in models:
        if m.get('source') != 'tiberius' or not m.get('cds_sha256'):
            continue
        matches = annotations.get(signature(m), [])
        identities = {(r['label'], r['original_gene']) for r in matches}
        if len(identities) != 1:
            continue
        result = deepcopy(matches[0])
        result.update(status='annotation_structure_match', basis='exact_native_cds_and_sequence',
                      homology_names=_homology_names(m))
        result['database_name_differs'] = bool(result['homology_names']) and all(
            n.casefold() != result['label'].casefold() for n in result['homology_names'])
        result['notes'].append('Prediction has the same phased native CDS, stop blocks and CDS fingerprint as the annotated model.')
        results[m['model_id']] = result
    return results


def export_identities(root, records, identities, id_field):
    """A separate derived table; retain the exact source/decision identifiers."""
    from .sequence import write_tsv
    rows = []
    for record in records:
        value = identities[record[id_field]]
        rows.append(dict(record_id=record[id_field], genome_id=record.get('genome_id', ''),
            region_id=record.get('region_id', ''), copy_identity=value['label'],
            copy_identity_status=value['status'], copy_identity_basis=value['basis'],
            original_annotation_name=value['original_name'], original_gene=value['original_gene'],
            ncbi_gene_id=value['ncbi_gene_id'], family_hint=value['family_hint'],
            protein_database_names=';'.join(value['homology_names']),
            database_name_differs=value['database_name_differs'],
            orthology_confirmed=False, function_status=value['function_status'],
            sources=';'.join(value['sources']), notes=' '.join(value['notes'])))
    fields = list(rows[0]) if rows else ['record_id','genome_id','copy_identity','copy_identity_status']
    write_tsv(root/'gene_identities.tsv', rows, fields)
