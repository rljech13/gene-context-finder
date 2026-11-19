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

No folding jobs, model-weight downloads, public MSA submissions or network
requests occur in this command. Existing exact structures are imported with
`--structures structures.json`; a later folding/cache stage can produce this
same catalogue. Prototype imports one full PDB CA chain per protein; partial
experimental structures and mmCIF need a later extension.

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
Omit `confidence` when B factors are not pLDDT. PAE is not imported or assessed.
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
