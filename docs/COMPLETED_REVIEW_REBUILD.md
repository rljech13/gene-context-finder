# Rebuild a completed report in the current review form

`python -m genecontext.legacy_review --results ORIGINAL --out NEW --title GENE`
creates a new presentation. It never modifies the original result directory.
The original prepared model report must match its seal, selected model IDs must
match the consumed decision snapshot, and every imported protein/CDS must match
its native SHA-256. `--confirmed-decisions FILE` accepts a later final-choice
file only when both original review ID and sequence-set ID match.

The new form restores model cards, the exact saved protein alignment, the
separate STOP3 synteny view, final choices, and post-synteny reference approval.
Codon display is back-translated from exact native CDS over the saved protein
alignment columns; existing scientific alignment files remain unchanged.
No annotation, MAFFT, tree or folding job runs during migration. Historical
STOP4 trees can be included as separately labelled evidence. An unavailable
tree record never supplies an invented topology.

Original review IDs, set ID, comments, choices and sequence hashes are recorded
in `presentation_rebuild.json`. Derived review IDs isolate new drafts from
original browser storage. Matching old browser drafts and final choices can be
imported without overwriting their original keys. Imported proposals retain
their original review status; reference approval is always a separate action.

The Python `package()` helper adds verified platform launchers, cached analysis
records, the matching STOP4 backend wheel and `Build_STOP4.py`, retaining source
scientific files and original UI/seals. Both the standalone local kit and full
archive are verified member by member. The original archives remain available.
Hosted completed pages display the imported results; explicit local launcher
sessions enable subsequent recalculation. Existing active studies are untouched.

The pinned templates in `vendor/legacy_review/SOURCE.json` identify the current
review UI snapshot and the adapted files used for this compatibility renderer.
