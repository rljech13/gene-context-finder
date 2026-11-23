# Report presentation and portable review

The generic package generates the same current presentation for arbitrary
gene/family labels. The pluripotency deployment is an application, not a
required server or input catalogue.

The review order is model choice, saved protein/CDS alignment, tree and outgroup,
confirmation, then STOP3 synteny and final outcomes. STOP3 is a distinct step
on the same page. Existing `?view=synteny` links open its section. Full model
accessions and stable short identifiers are shown consistently in the model
list, alignment, outgroup controls and final dossier. Exon/splice marks use the
provided native CDS coordinates; protein-only references remain unmarked.

Alignment and tree comments belong to exact result identifiers and are carried
into the STOP4 handoff and dossier. Changing a model or rebuilding an alignment
does not silently reuse an old final decision. No calculation starts on opening
an HTML file or launching the viewer.

## Viewing and calculation

| Operation | Requirements |
| --- | --- |
| Read a generated report | A current web browser and its linked result files |
| Interactive local alignment/tree review | Extracted portable package and matching native launcher; no Python, Conda or research server |
| New genome search and model preparation | Python package plus the external scientific tools and suitable predictor models described in the setup guide |

Build a launcher with `python scripts/build_review.py --out RUNTIME`. Set
`GCF_REVIEW_BIN_DIR` to that directory and run `gene-context review-package`
on a saved model-review workspace. The default packages the current platform;
build with `--platform all` and package with `--platform all` to include the
six configured OS/CPU launchers. Move the entire extracted folder, not only
its HTML. On a different OS/CPU, include that destination's launcher.

The viewer uses a suitable MAFFT on PATH, or downloads the pinned portable
distribution only when alignment is requested. Other explicitly requested
calculations similarly require their scientific engines. Architectures without
a compatible portable engine require a native installation. This is separate
from the ability to read the report or restore saved decisions. Full analysis
on Windows uses WSL; the portable viewer runs natively on Windows.

## Saving and recovery

Small choices retain their original browser keys. Large alignment, tree and
one-step undo checkpoints use IndexedDB, with migration of matching older
localStorage records. Conflicting older copies are retained. If persistent
browser storage is unavailable, file viewing remains available with a notice
to export changes before closing.

The native viewer saves this combined state to `.review/session-state.json`
in the extracted package folder. It restores that file before report scripts
run, including on a new localhost port or in a new browser context. The file
is bound to the exact dataset checksum. Another gene's state is excluded.
The save indicator reports pending writes and errors; Finish and clean up
waits for saving. A session file is limited to 128 MiB.

## Verification

`tests/installed_report_smoke.py` must run against an installed wheel from
outside the checkout. It checks the current rendered assets, exact handoff,
STOP4 comments, package contents and preservation of source checkpoints.

`tests/report_portability_browser.cjs PROOF_DIRECTORY` launches the compiled
viewer in Chromium, runs MAFFT on synthetic models, checks migration, saves
an 18 MiB capacity-control record, reloads, shuts down and restarts on a new
origin/browser context. It verifies recovered alignment/comments, preservation
of another gene's choices and absence of duplicate calculations. Desktop and
narrow screenshots are saved with the verification record.

The `installed-reports` CI matrix runs this check on Linux, macOS and Windows.
Install browser-test dependencies with `npm ci` and `npx playwright install chromium`.
These controls exercise software behavior; they do not measure biological
specificity or prediction accuracy on new genes or taxa.
