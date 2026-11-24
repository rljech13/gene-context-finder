# Final reference approval after synteny

Current order: Stage 1 loci → preparation including Tiberius → Stage 2 model choice,
alignment, tree and confirmed-set handoff → **Stage 3 synteny/context review** →
Stage 4 dossier/structures. The confirmation inside Stage 2 hands a model set to
synteny; it does not authorize AlphaFold.

At the end of the separate synteny page, each confirmed target/candidate model
has **«Approve model as reference»**. This explicit action selects the final
reference after context review; at most one alternative model of a physical
locus can be approved at a time. Distinct physical copies remain distinct.
Revoke the old reference before choosing another variant of that same locus.
Queries and neighbor-only anchors have no approval action.

**«Prepare AlphaFold Server submission»** is disabled until approval. It
verifies the full native protein hash and unchanged confirmed-set handoff,
then downloads one monomer Server JSON; upload and launch happen manually on
the official Server. Opening the page or approving a reference never submits
a job. Server result handling still preserves all five returned samples and
shows one best structure.

Neighbor-pair review in Stage 3 stays separate from the Stage 2 model-set handoff.
Changing these context decisions requires fresh reference approval, without
sending the user back through model selection or rerunning the tree.

Approval is stored under `gcf-synteny-reference-v1-<source-review-id>`, separately
from model decisions, alignment/tree selections, outgroups and undo. Old choices
are never rewritten or promoted automatically. A changed model set, comments,
source/protein identity, revoked approval or mismatched Stage 2 handoff blocks
preparation. The downloaded `stage3.synteny.references.json` carries the confirmed
handoff, separately recorded Stage 3 neighbor-pair choices and exact approved subset for `af3-prepare` / `af3-campaign`.

The additive `genecontext.synteny_reference.attach(html)` publisher hook appends
one independent script after the separate synteny renderer. It is idempotent,
does not rewrite embedded scientific data and stays inactive in Stage 2. Existing
pages should be patched from their latest contents with a compare-before-write
check; never rebuild a live page from an older full HTML snapshot. Concurrent
renderers and selection/undo code must be preserved byte-for-byte.

Completed legacy studies retain their native exports and explicit historical
approval route to Stage 4. Their alignment/comparison sheet no longer presents
itself as the new synteny stop. User-confirmed historical decisions are not
withdrawn by moving the UI button.
