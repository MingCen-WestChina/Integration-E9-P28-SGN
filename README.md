# Integration E9–P28 SGN

Analysis and figure code for the developmental atlas of mouse auditory neurons.
The atlas contains 7,355 embryonic cells and 985 postnatal cells, covering E9.5
through P28. Short stage labels in filenames follow the original analysis files.

## Before running

This is a collection of the study's analysis scripts and notebooks. Configure
each input and output location before running the relevant section. A setting
such as `INPUT = ROOT` is a placeholder: replace that occurrence with the full
location of the required file. Different placeholders usually refer to different
files or directories. Comments identify expected intermediate files and their
producing scripts. Set executable and helper-script locations in the same way.

The code retains the original analysis parameters and validation checks. Fixed
cell membership and shared-gene metadata are provided as readable TSV/JSON files
under `fig1/metadata`; the original checksums are still checked when loading them.
The TF reference used for the nine-segment analysis is under `fig3/metadata`.
Per-run filenames required by STAR, FastQC, HDF5, archives, and temporary-file
protocols remain part of those formats.

Candidate-analysis `Inputs.h5` and its formal count H5AD are provided in the
release asset linked in the Fig4 section below. Other expression matrices,
downloaded sequencing reads, fitted models, processed H5AD objects, and cached
result stores must be supplied separately. The panel notebooks read these
intermediate results. Run the corresponding upstream analyses first when
rebuilding them. This repository does not contain the original Illustrator
layout edits.

E10 count export also takes the approved 73-cell selection H5AD and an export
contract text file (`SCHEMA_FILE`, used for provenance hashing). The public
accession manifests do not replace these inputs.

## Code map

| Folder | Entry points and purpose |
| --- | --- |
| `fig1` | Stage clustering, embryonic CSS integration, postnatal scVI integration, whole-atlas CSS integration, UMAP and marker dot plots. |
| `fig2` | Monocle3 and Slingshot, stage comparisons and GO analysis, marker trends, and the associated supplementary plots. |
| `fig3` | Accepted global tree, local CVG tree, conditional graft, path and segment gene analyses, and individual panels A–J. |
| `fig4` | SCENIC, regulon dynamics, RNA-associated modules, human genetics, VGN analyses, and individual Fig4/FS4/FS5 panels. |

### Fig1

1. `fig1_E10_download.py`, `fig1_E10_pair_validation.py`, and
   `fig1_E10_ena_sequence_sha256.py` establish the approved E10.5 FASTQ cohort and
   its sequence baseline. `fig1_E10_recover_quality.py` recovers full-quality
   Normalized SRA reads. Public accession manifests are in `metadata/e10`.
2. `fig1_E10_build_index.sh` builds the matching GRCm38/GENCODE M25 STAR index.
   `fig1_E10_gtf_bed12.py` supplies the transcript reference for RSeQC.
   `fig1_E10_trim_pilot.py` compares trimming settings;
   `fig1_E10_full_counts.py` reuses the nine verified pilot results and counts the
   complete 73-cell cohort. `fig1_E10_Counts.py` checks the count matrix and
   exports the count-based stage object.
3. `fig1_E9_clustering.ipynb` through `fig1_E18_clustering.ipynb` retain the
   approved embryonic selections. These files use the final reviewed source
   versions; the review metadata are external readable files.
4. `fig1_P1_clustering.ipynb`, `fig1_P7_clustering.ipynb`, and
   `fig1_P28_clustering.ipynb` prepare the postnatal selections.
   `fig1_Postnatal_preparation.py` harmonizes these inputs using
   GRCm39/Ensembl release 115 and writes the prepared postnatal object.
5. Run the embryonic, postnatal, and whole integration notebooks with their
   required stage objects. `fig1_IntegrationPlotModify.ipynb` produces the final
   integration and marker displays.

The reviewed embryonic cell counts are 12, 73, 331, 818, 1,674, 161, 3,028, 71,
and 1,187 for E9, E10, E11, E13, E14, E15, E16, E17, and E18, respectively.
P1, P7, and P28 contribute 314, 387, and 284 cells.

### Fig2 and associated supplements

`fig2_moncle3-1.ipynb` and `fig2_slingshot-2.ipynb` infer the two trajectories.
`fig2_embry_heatmap.ipynb` and `fig2_Stage-GO.R` implement stage comparisons and
enrichment. `fig2_slingshot_genealogy.ipynb` keeps the endpoint-constrained model
and its marker-trend plots together in labelled sections.
`FS2_embry_violin.ipynb` and `FS3_hierarchical_clustering.ipynb` provide the
associated supplementary displays. This folder contains seven code files.

### Fig3

1. `fig3_tree_preparation.py` reads the accepted embryonic CSS integration
   (`CSS_3000_Embryonic_Inte_adata.h5ad`) and writes `prepared.pkl`, using the
   original PCA30 and global RMS normalization.
2. `fig3_global_tree_analysis.py` reads `prepared.pkl` and produces
   `selected.pkl` and the original `joint-tree.png` rendering. It reproduces the
   accepted 7,355-cell cold-start model: 160 initial nodes, lambda 3, sigma 0.095,
   1,000 iterations and seed 1. The fitted tree contracts to 153 nodes.
3. The first section of `fig3_path_gene_analysis.ipynb` packages these outputs
   as `tree_results.pkl`. Its later sections retain the complete path GAM
   engine, path-specific selection, source checks and Ic endpoint edgeR tests.
   `REBUILD_PATH_ANALYSIS` explicitly runs those original functions and writes
   their working statistical tables. The default figure sections read the
   separately supplied formal `fig3CG_results.h5` and its validation JSON;
   rebuilding the working tables does not silently replace that accepted store.
4. `fig3_local_CVG_analysis.py` uses the original GSE178931 QC counts and
   `tree_results.pkl`. Its 1,978 curated reference records are readable in the
   source, so it does not execute or read an unpublished notebook. Set
   each of the seven `QC_FILES` locations independently and set
   `LOCAL_FIT_SCRIPT` to the included `fig3_local_CVG_fit.py`. These scripts
   prepare and fit the 478-cell CVG component in `local_results.pkl`; the
   accepted local key is HVG2000, 20 nodes, lambda 5 and seed 42, with sigma 0.1
   and 1,000 iterations.
5. `fig3_tree_graft_analysis.ipynb` contains the original `build_stitched`
   registration and conditional graft. `rebuild_from_upstream_components()`
   reads `tree_results.pkl` and `local_results.pkl` to reconstruct its numerical
   component. The saved-result and final display sections use the supplied
   `tree_results-1.pkl` archive and preserve its original component checks.
6. `fig3_segment_gene_analysis.ipynb` reads that final tree, the original CSS
   count object, native source denominators and the unchanged TF TSV in
   `metadata`. Set `REBUILD_ANALYSIS=True` only to recompute the complete
   nine-segment association and five-path smoothing. It writes
   `Trajectory dot plot_results-2.h5` and `Trajectory dot plot_genes-2.tsv`.

The path-GAM, Ic-edgeR and nine-segment R engines retain their original source
bytes. In particular, the nine-segment engine still has the exact checksum
required by the unchanged H/I panel readers. Files starting with `fig3A_`
through `fig3J_` draw the individual panels from the corresponding accepted
result stores. These files provide the original analysis steps and cached-result
workflows, rather than a single command that rebuilds every publication archive.

### Fig4, FS4, and FS5

`fig4_SCENIC_analysis.ipynb` and its released worker implement network inference,
motif pruning, and AUCell scoring. `fig4A_regulon_dynamics_analysis.py`, its
global and early-fork engines, and `fig4_regulon_supplement_analysis.ipynb`
implement the accepted regulon-dynamics analysis and additional predefined
signatures. The primary model uses the accepted `tree_results-1.pkl`, with
7,355 measured cells and 7,348 modeled cells; the seven connector cells do not
define a biological spline state. Fixed AUCell/signature inputs and the retained
late-fork results are required separately. The supplied native-scoring and
late-fork support engines retain their original methods.

`fig4_candidate_regulator_analysis.ipynb` and
`fig4_SGN_program_analysis.ipynb` implement candidate selection, modules,
seven-library holdouts, all 19 program models, SCENIC activity, and the final
Neurod1 reference/project-cohort analyses. Run candidate analysis before the SGN
program analysis: `Candidates.h5` is produced from the audited `Inputs.h5` and
formal common-count input. The SGN analysis then produces `Activity.h5`,
`Fig4B-D.h5`, and `SGN-lines.h5`. Its plain Python/R workers and embedded exact
source are both supplied.

`Inputs.h5` and its matching formal count H5AD are provided in
[fig4_candidate_inputs.zip](https://github.com/MingCen-WestChina/Integration-E9-P28-SGN/releases/download/v1.0.0/fig4_candidate_inputs.zip)
from GitHub Release `v1.0.0`. The ZIP preserves these project-relative locations:

- `Processed Data/07_Candidate_regulors/Inputs.h5`
- `Processed Data/04_Embryonic_cluster&Inte/CSS_3000_Embryonic_Inte_adata.h5ad`

Extract the ZIP into your chosen project root, keeping both `Processed Data`
subfolders intact. This preserves the project-relative `common_counts_file`
reference stored in `Inputs.h5`. Configure `PROJECT_ROOT` and every source
`ROOT` placeholder for your installation, including input/output locations,
the released notebook's own location, executable paths, and a work directory.
The archive includes `DATA_README.txt` and `SHA256SUMS.txt` for these two files.

`Inputs.h5` is an audited intermediate package containing native TF RNA/count
and measurement fields, AUCell values, metadata, motif sets, reviewed program
annotations, and tree/embedding correspondence. The complete gene-count matrix
is read from the accompanying H5AD. The standalone assembly script for
`Inputs.h5` was not located in the retained project source.
This archive supplies only the candidate-analysis input pair for 7,355 embryonic
cells (E9.5-E18.5) and 11,467 genes; it is not a complete package of paper results
and contains no P7 data. The complete downstream candidate algorithms are
included; this is not a raw-data-to-panel automatic workflow.

`fig4_human_genetics_analysis.ipynb` contains the full acquisition, MAGMA/scDRS,
FinnGen/GTEx colocalization, prior/tissue sensitivities, and ortholog workflow.
Its identity helpers are `fig4_Integration_genes.py` and
`fig4_Integration_axis.py`; the shared trajectory source is the released
`fig3_segment_gene_analysis.ipynb`.

The separate VGN sequence is `fig4_VGN_input_analysis.py` →
`fig4_VGN_trajectory_analysis.py` → `fig4_VGN_network_analysis.py` →
`fig4_VGN_modules_analysis.py` → `fig4_VGN_program_analysis.ipynb`.
These stages produce `Early-data.h5`, the accepted 1,010-cell `Trajectory.h5`,
`Network.h5`, `Modules.h5`, and `VGN-results.h5`, respectively. The superseded
diffusion-pseudotime fitting workflow is omitted from the input producer.

The individual Fig4/FS4/FS5 panel notebooks consume these result stores.
`FS4B_Neurod1.ipynb` defaults to the explicit `original_archive` source-identity
profile for the supplied frozen cache. Select `released_rebuild` only for
results rebuilt from the supplied upstream source after configuring its path
placeholders. Each profile checks one exact set of six worker fingerprints;
the profiles do not accept alternative hashes within a check.
If editing `ROOT` changes a `PAPER_B_SOURCES` string, derive the
`released_rebuild` profile's exact fingerprints from the reviewed configured
strings while retaining its strict assertions, because the template fingerprints
no longer apply.

## Inputs and software

Public expression datasets are GSE178931, GSE231513, GSE165502, GSE195500, and
GSE202920. The P28 inputs are GSM6138444 and GSM6138445. P7 is the study's newly
generated dataset and must be supplied separately. The human-genetics analysis
uses FinnGen R13 and GTEx v8. SCENIC requires the matching mouse TF catalog,
motif rankings, and motif annotations specified in the source.

The main Python dependencies include NumPy, pandas, SciPy, AnnData, Scanpy,
scikit-learn, scvi-tools/PyTorch, SimplePPT, NetworkX, pySCENIC/arboreto/ctxcore,
h5py, and Matplotlib. Individual sections use additional packages named in their
imports. R sections use Monocle3, Slingshot, mgcv, edgeR, coloc, and the annotation
and enrichment packages loaded by those sections. Some acquisition/counting
scripts require Linux command-line tools and cannot run as Windows-only scripts.

Study versions recorded in the manuscript include Cutadapt 5.1, STAR 2.7.11b,
Monocle3 1.4.27, Slingshot 2.18.0, edgeR 4.8.2, mgcv 1.9.4,
pySCENIC 0.12.1, and arboreto 0.1.6. The workflows were developed in separate
environments; this list is not a single environment lockfile.

The local release was checked for syntax, retained parameters, source provenance,
empty notebook outputs, and readable reference-metadata integrity. Model fitting
and full figure regeneration were not rerun during code preparation.
