"""Retain BLAST's own gapped HSP strings, without reconstructing an alignment."""


def alignment_fields(fields, query, subject_step=1):
    if len(fields) == 12:
        return {}
    if len(fields) != 14:
        raise ValueError('Expected 12 legacy or 14 alignment BLAST fields')
    qseq, sseq = fields[12:14]
    qa, qb, sa, sb = map(int, fields[4:8])
    if (len(qseq) != int(fields[3]) or len(sseq) != len(qseq)
            or any(a == b == '-' for a, b in zip(qseq, sseq))
            or len(qseq.replace('-', '')) != abs(qb-qa)+1
            or len(sseq.replace('-', ''))*subject_step != abs(sb-sa)+1):
        raise ValueError('BLAST alignment strings do not match HSP coordinates/columns')
    # All stage-1 queries are searched on their forward orientation.
    if qa > qb or qseq.replace('-', '').upper() != query[qa-1:qb].upper():
        raise ValueError('BLAST alignment differs from the exact supplied query')
    identity = sum(a.upper() == b.upper() and a != '-' for a, b in zip(qseq, sseq)) / len(qseq)
    if abs(identity*100-float(fields[2])) > .0011:
        raise ValueError('BLAST alignment identity differs from reported identity')
    return dict(alignment=dict(qseq=qseq, sseq=sseq,
        query_start=qa, query_end=qb, subject_start=sa, subject_end=sb,
        subject_step=subject_step, provenance='saved_blast_output'))
