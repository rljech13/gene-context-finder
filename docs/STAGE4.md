# STOP 4 prototype — candidate dossier

STOP 4 reads a completed, sealed STOP 3 result. It does not search genomes,
predict genes, resume queues, change STOP 2 selections or approve orthology.
The main view is one compact primary-model table. Extra copies are in a
collapsed separate-cases table. Every row separates target identity from
coding-model quality, with the exact source, locus, choice reason and limits.

```sh
gene-context stage4 --results STOP3 --decisions stage3.final.decisions.json --out STOP4
```

Without `--decisions`, the command uses STOP 3's saved finalization state.
Unreviewed primary choices are rejected unless `--preview` is explicitly set.
Preview reports have an unmistakable banner and never become manual approval.
The saved reference panel and source annotation can support a candidate, but
this prototype does not automatically resolve target/paralog clades or mRNA
existence. A complete predicted CDS is not an experimentally confirmed transcript.

Output: offline `report.html`, `final_table.tsv`, `special_cases.tsv`, exact
primary/case protein and CDS FASTAs, per-model sequences, per-group trees and
alignments, optional structures/measurements, compact evidence, manifest and seal.
Genomic windows and thousands of intermediate alignment files are not copied.

## Gene trees

Each comparison group is independent. Saved selected models, separate cases and
supplied protein queries are used; discarded STOP 2 models are not reintroduced.
Identical full sequences share a leaf with all source aliases retained.
The existing MSA is reused unless additional references require MAFFT.
After deduplication, default column occupancy is 0.5. Empty taxa after masking
are flagged and omitted. The mask and aliases are exported. Prototype minimums
are 4 unique taxa, 30 retained columns, 5 variable columns; these are software
admission rules, not guarantees of phylogenetic information.

IQ-TREE chooses a substitution model with ModelFinder (BIC), estimates ML
topology with UFBoot (1000 replicates, BNNI, seed 42), and exports Newick/SVG.
`--threads 2`, `--bootstrap 1000`, `--occupancy .5`, `--iqtree PATH` can be set.
`--tree none` explicitly skips trees. Missing tools are reported as missing;
tool failures stop the build instead of masquerading as completed trees.

Optional `--reference-panel panel.json`:

```json
{"schema":"stage4_reference_panel_v1","entries":[
  {"id":"related_copy","group":"GENE","gene":"PARALOG","species":"Species name",
   "sequence":"MA...","role":"paralog","source":"Verified accession/version"},
  {"id":"outer_copy","group":"GENE","gene":"GENE","species":"Outside ingroup",
   "sequence":"MA...","role":"outgroup","source":"Verified accession/version",
   "root_basis":"Why this homolog is outside this ingroup"}
]}
```

Replace the illustrative sequences with actual canonical amino-acid sequences.
Roles are `target`, `paralog`, `outgroup`. IDs must be distinct. Outgroups require
a biological basis and cannot collapse with an ingroup sequence. Rooting follows
this explicit per-family profile; otherwise the tree is unrooted. Display order
does not define a root. Nonreversible root tests and duplication/loss
reconciliation are deliberately not part of this lightweight prototype.

## Interactive 3D structures

With the pinned local 3Dmol.js bundle, imported structures have a **View 3D**
link. The included viewer supports rotation, zoom, cartoon/ball-and-stick/CA
representations, local structural-similarity or separate pLDDT coloring, and PDB download. Coordinates are
embedded locally; there are no runtime CDN or structure-service requests.
Identical proteins share a viewer option; selected models and query references
are labeled separately. The default colors show Cα deviation after the actual
US-align superposition against the named reference. Fixed distance bins are
<1, 1–2, 2–4, 4–8 and ≥8 Å; residues not compared are gray. The reference can
be changed or optionally shown as a translucent overlay. Prediction confidence
is a separate mode and is never relabeled as similarity.

Only residues retained at pLDDT ≥70 in each input participate in this prototype's
similarity map. Both full-protein coverages and the two TM-scores are shown on
the picture. TM-score normalization uses the retained lengths, not full proteins.
Residue mappings are validated against both input sequences; every derived Cα
distance must agree with the US-align `-do` output within 0.0011 Å. The matrix,
raw output, original coordinates, engine hash and mappings are retained. These
colors are deviations under one rigid fit, not a domain-independent local score.

New builds retain both the CA coordinates used for structural comparisons and
the original chain's atomic coordinates for viewing. A CA-only input is explicitly
labeled as a backbone trace rather than displayed as an atomic model.

To attach the viewer to a sealed prototype and reproduce its structural pairs
with explicit transforms (trees, sequences and manual choices stay unchanged):

```sh
gene-context stage4-viewer --results STOP4 --out STOP4_WITH_3D --structures structures.json --usalign /path/to/USalign
```

The optional structure catalogue can restore full atomic coordinates. Every CA
residue, position and confidence must match the structure already measured in
STOP 4. The source result and its manual choices stay unchanged. The viewer bundle
and license are shipped under `genecontext/vendor/`, pinned by hash in
`3Dmol-source.json`. Browser WebGL support is required for interactive rendering.

## Exact structures

### Preferred new predictions: AlphaFold Server

Folding export is allowed **only after STOP3 and explicit approval of a final
primary reference model** using **«Утвердить модель как эталон»** in the stage-3
primary-target table. A saved reviewed choice alone, saved candidate choices, a
completed alignment or a proposed primary are insufficient. The approval is bound to the exact protein
SHA-256 and reviewed final gene/genome decision. Alternatives, special cases,
queries and unreviewed suggestions are not sent automatically.

Prepare every explicitly approved final reference by default; identical full
sequences share one job. Show **one best structure per protein**. AlphaFold Server calculates
five samples per job; its public upload format does not expose a setting for
three samples. The original ZIP preserves all raw samples as evidence.

For a whole multi-gene/multi-species collection:

```bash
gene-context af3-approval --results STOP3 --decisions final.json --out final.references.json
gene-context af3-campaign --fasta all.selected.proteins.faa \
  --approval final.references.json \
  --metadata models.by_sequence_id.json --out SERVER_CAMPAIGN \
  --completed prior/verified-result --first-batch-size 28
```

The optional metadata object is keyed by FASTA ID and may contain `gene`,
`species`, `original_model` and source locators. Completed imports are checked
against their seal, archive and exact coordinate sequence before reuse. Pending
proteins are distributed across genes and split into batches of 30 (at most 100
draft jobs per upload); the first batch can match the current remaining quota.
`jobs.tsv` maps every selected model to its unique protein job. `campaign.json`
records completed, pending and unsupported proteins. Preparing batches does not
submit them or create a daily schedule.

Use the exact selected proteins from the saved review, not a replacement sequence
looked up by gene name. The Server's documented JSON upload is the initial
transport; the prototype does not rely on internal browser RPC endpoints.

```bash
gene-context stage4 --results STOP3 --decisions final.json --out STOP4_SERVER \
  --structure-provider alphafold-server --structure-limit 2
# Or prepare independently from an exact selected FASTA:
gene-context af3-prepare --fasta final.proteins.faa --approval final.references.json --out SERVER_REQUESTS --limit 2
```

In AlphaFold Server, upload `upload.server.json`, open each saved draft and submit
it. One exact monomer is one job. Inputs retain seed 42; the Server dialect uses
string seeds and accepts 4–5000 canonical residues. Unsupported proteins are
reported without substituting unknown residues. Identical proteins share a job.

After completion, download the original job ZIP and import it into a new directory:

```bash
pip install 'gene-context-finder[structures]'
gene-context af3-import --requests SERVER_REQUESTS/requests.json \
  --archive downloaded-job.zip --out SERVER_RESULT \
  --source-url https://alphafoldserver.com/fold/JOB_ID
```

The small optional Gemmi parser is needed only for this mmCIF import; ordinary
pipeline installation stays dependency-free. The importer checks the returned
job name, seed and full protein against sealed input files, selects the highest
returned ranking score, validates exact coordinate sequence and PAE dimensions,
and emits the ordinary `stage4_structures_v1` catalogue. Rounded score ties use
the lower Server confidence-rank index (`model_0` first), so the one displayed
structure follows the Server's ranking. All returned samples,
MSAs, templates and output terms remain in the original `source.zip`; selected
raw mmCIF, confidence files and conversion details are also retained. Chain-A PDB
conversion changes decimal precision; the raw structure is preserved. The
catalogue can be supplied with `stage4 --structures SERVER_RESULT/catalogue.json`.
Use a new STOP4 directory and preserve the existing report and decisions.

### Extension boundary

The core consumes `stage4_structures_v1`; provider adapters are separate modules
(`alphafold_db.py`, `alphafold_server.py`). New providers or Server formats can
add an adapter without changing selected model IDs or review decisions. Request
and result manifests carry provider/model, adapter/schema versions, exact protein
SHA-256, source model IDs and an `extensions` object. Request files, original
archives, selected coordinates and outputs are hashed; old attempts are retained.
Unknown returned Server request versions are recorded rather than silently
rewritten; supported monomer/coordinate/metric checks still apply. A failed import
leaves an explicit failed attempt, never a ready receipt. This boundary supports
a later official submission API without making it a requirement for a first trial.

Two real Server jobs for user-selected Gallus gallus (annotation) and Pavo
cristatus (Tiberius) POU5F3 proteins completed on 2026-11-19. Each input was 389 aa
with seed 42; five samples were retained per job. This validates a small real
round trip, not biological functional equivalence or all future service formats.

### Experimental API preparation and retrieval

STOP4 can explicitly look up existing exact-sequence monomers through AlphaFold
DB, or prepare requests for genuinely new AlphaFold 3 predictions. They are
separate operations and have separate statuses:

```sh
gene-context stage4 --results STOP3 --decisions final.json --out STOP4_WITH_AFDB \
  --structure-provider afdb --structure-limit 3
gene-context stage4 --results STOP3 --decisions final.json --out STOP4_AF3_REQUESTS \
  --structure-provider af3-prepare --structure-limit 3
```

`afdb` uses the documented sequence-MD5 endpoint and verifies the complete amino
acid sequence in both the API metadata and PDB coordinates. It downloads a
versioned existing model, not a new folding calculation. Coordinates, metadata,
source accession/organism and optional validated PAE are retained. The pilot
checks at most three unique proteins by default (up to 50 explicitly), with
primary models first. Identical sequences reuse the same cache artifact. The
initial retrieval limit is 1500 amino acids and PDB chain A; unsupported lengths,
missing exact models, API failures and budget omissions are shown explicitly.
No similar protein is substituted. A structure from another accession or species
can be reused only for the same full sequence and keeps its true provenance.

`af3-prepare` writes per-protein native AF3 JSON (for local inference or a
configured Google Cloud endpoint) and AlphaFold Server upload JSON. Status is
`awaiting_prediction`, never “predicted”. It starts no cloud resources, sends
no prediction requests and consumes no paid credits. Endpoint/account setup is
a separate choice; see [service access comparison](STRUCTURE_SERVICES.md).

These steps can also be run separately against an exact protein FASTA:

```sh
gene-context structure-fetch --fasta STOP3/selected.proteins.faa --out AFDB_CATALOGUE
gene-context af3-prepare --fasta final.proteins.faa --approval final.references.json --out AF3_REQUESTS --limit 2
```

The AFDB catalogue is accepted by the ordinary offline `--structures` importer.
Use `--structure-cache PATH` for a shared verified cache, `--structure-offline`
to forbid requests, and `--structure-refresh` to explicitly retain a fresh
response alongside the pinned cache. Network errors are not cached as missing
models. Verified not-found responses remain pinned until explicit refresh.
PAE matrices are dimension/value checked and downloadable; domain interpretation
and a PAE plot are outside this prototype. pLDDT remains separate from measured
structural similarity. The default `stage4` command is offline. Both structure
providers work through `project.json` stage4 options without changing saved
STOP1–3 decisions. Downloading or predicting a structure never approves a gene.

### Import an existing catalogue

With the default provider, no folding jobs, model-weight downloads, public MSA
submissions or network requests occur. Existing exact structures are imported with
`--structures structures.json`; a later folding/cache stage can produce this
same catalogue. The core imports one full PDB CA chain per protein; the Server
adapter converts and verifies mmCIF first. Partial experimental structures and
arbitrary complex/mmCIF imports need a later extension.

```json
{"schema":"stage4_structures_v1","entries":[
 {"sequence_id":"MODEL_ID_OR_reference_q0001","file":"model.pdb","chain":"A",
  "file_sha256":"SHA256_OF_INPUT_PDB","protein_sha256":"SHA256_OF_FULL_PROTEIN",
  "source":"AlphaFold DB","accession":"UNIPROT_ID","source_url":"https://...",
  "confidence":"plddt_bfactor"}
]}
```

The actual CA-chain sequence must equal the entire exact saved protein, not just
the supplied checksum. No homolog structure is relabeled as the selected model.
Omit `confidence` when B factors are not pLDDT. Optional PAE is verified and retained;
it is not interpreted as a domain annotation.
US-align (`--usalign PATH`) measures all available same-group pairs for full
chains and independently filtered pLDDT≥70 coordinates. Both normalized TM-scores,
RMSD, aligned residues and coverage of both original full proteins are exported.
The main table links the best measured confident-reference comparison, with
same-sequence controls excluded from independent support. Filtering is not
domain annotation; disordered regions are not called incorrect genes.

“Consistent candidate” uses saved target-reference comparisons with at least
50% coverage on both sides, compatible specific source naming and no known
alternative-copy conflict. It is a transparent prototype label, not a
calibrated probability or orthology assertion. Missing competitors, reviewed
neighbor anchors, transcripts and structures remain explicit missing evidence.
