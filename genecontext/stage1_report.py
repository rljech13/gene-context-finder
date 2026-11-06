"""Offline, locus-first review of existing BLAST evidence and early synteny."""
import json
from pathlib import Path

TEMPLATE = Path(__file__).with_name("stage1_review.txt").read_text()

def write_report(root, data):
    payload = json.dumps(data).replace('<', chr(92)+'u003c').replace('&', chr(92)+'u0026')
    (root/'report.html').write_text(TEMPLATE.replace('__DATA__', payload))
