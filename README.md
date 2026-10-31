# Integration E9–P28 SGN

Analysis source for the developmental atlas of mouse auditory neurons, covering
7,355 embryonic cells and 985 postnatal cells from E9.5 through P28. Short stage
labels in filenames follow the original project files.

## Organisation

Figure code is grouped by manuscript panel. Each panel has one Python notebook.
An additional R file is retained when the original analysis used a separate R
program. Shared original source is copied into each relevant panel notebook;
this duplication keeps the source for a panel together.

The original analysis and plotting blocks, functions, parameter values, worker
source strings, and statistical checks are retained. Changes concern file
organisation, explanatory headings, English comments, and local path placeholders.
Some original functions produce several panels together. Those functions remain
intact, and the notebook heading identifies which panel the file corresponds to.

Set every file or directory placeholder, such as `INPUT = ROOT`, to its own local
location before using the corresponding section. Comments identify the expected
intermediate inputs. Different placeholders generally refer to different paths.
Original per-run tool filenames, archive keys, and HDF5 keys are retained.

## Fig1

The `fig1` folder retains stage clustering, E10.5 raw-read recovery and counting,
postnatal preparation, embryonic and whole-atlas CSS integration, postnatal scVI
integration, and the final integration/marker plots. Its files are unchanged by
the panel-code merge.

The nine embryonic stage notebooks use the final reviewed selections. Readable
cell, gene, and review metadata are under `fig1/metadata`, with the original
payload checksums retained. The stage cell counts are 12, 73, 331, 818, 1,674,
161, 3,028, 71, and 1,187 for E9, E10, E11, E13, E14, E15, E16, E17, and E18.
P1, P7, and P28 contribute 314, 387, and 284 cells.

For E10.5, the released download, paired-read validation, sequence-digest,
Normalized SRA recovery, STAR-index, BED12, trimming-pilot, full-count, and
count-object scripts retain the original processing settings. The public
accession manifests are in `fig1/metadata/e10`. Count export also requires the
approved 73-cell selection H5AD and its export-contract text input
(`SCHEMA_FILE`, used for provenance hashing).

`fig1_Postnatal_preparation.py` retains the GRCm39/Ensembl release 115 preparation
used for the accepted postnatal object. Follow the producer comments in the
integration notebooks when configuring their intermediate inputs.

## Fig2 and supplementary figures

| Panel | Files |
| --- | --- |
| Fig2A | `fig2A_Monocle3.ipynb` |
| Fig2B | `fig2B_Slingshot.ipynb` |
| Fig2C | `fig2C_Stage_heatmap.ipynb`, `fig2C_Stage_GO.R` |
| Fig2D | `fig2D_Genealogy.ipynb` |
| FS2 | `FS2_Stage_violin.ipynb` |
| FS3A | `FS3A_Stage_clustering.ipynb` |
| FS3B | `FS3B_State_clustering.ipynb` |

Fig2D keeps the original endpoint-constrained model and expression-trend code in
one notebook. Fig2C retains the original separate R enrichment program.
FS3A and FS3B each retain the original joint plotting source, which draws both
panels; their headings identify the intended panel and the final all-stage,
average-linkage/Pearson section. No single-panel drawing function was added.

## Fig3

| Panel | Notebook |
| --- | --- |
| Fig3A | `fig3A_tree.ipynb` |
| Fig3B | `fig3B_tree.ipynb` |
| Fig3C | `fig3C_tree_heatmap.ipynb` |
| Fig3D | `fig3D_tree_heatmap.ipynb` |
| Fig3E | `fig3E_tree_heatmap.ipynb` |
| Fig3F | `fig3F_tree_heatmap.ipynb` |
| Fig3G | `fig3G_tree_heatmap.ipynb` |
| Fig3H | `fig3H_Trajectory_dot_plot.ipynb` |
| Fig3I | `fig3I_Trajectory_dot_plot.ipynb` |
| Fig3J | `fig3J_Gene_expression_plots.ipynb` |

The corresponding original tree, local-CVG, conditional-graft, path-GAM,
Ic-endpoint, segment-association, and plotting source blocks are grouped within
these notebooks. Original embedded R engines keep their source bytes and strict
checksum checks. The readable TF reference remains in `fig3/metadata`.

## Fig4 and supplementary figures

The `fig4` folder contains the seven main-panel notebooks
`fig4A_regulon_dotplot.ipynb`, `fig4B_Nhlh1_activity.ipynb`,
`fig4C_Nhlh1_modules.ipynb`, `fig4D_Nhlh1_network.ipynb`,
`fig4E_Gene_expression_plots.ipynb`, `fig4F_candidate_dotplot.ipynb`, and
`fig4G_drug_evidence.ipynb`.

Supplementary notebooks are `FS4A_SGN_stage_expression.ipynb`,
`FS4B_Neurod1.ipynb`, `FS4C_Neurod2.ipynb`, `FS4D_Foxo6.ipynb`,
`FS5A_VGN_stage_expression.ipynb`, `FS5B_Tbx3_RNA.ipynb`,
`FS5C_Tbx3_module.ipynb`, and `FS5D_Tbx3_network.ipynb`.

Each notebook groups the relevant existing analysis source with its original
panel drawing code. This includes SCENIC/regulon dynamics, candidate and SGN
program analyses, human genetics and evidence curation, or the separate VGN
workflow as appropriate. Original shared programs can produce additional panels;
they have not been rewritten into new panel-specific algorithms.

## Original execution contexts

Read section headings before running cells. The notebooks combine original
programs that were used in separate execution contexts. Original standalone CLI
source sections retain their argument parsers, entry points, and subprocess
interfaces; use those sections in their original script context. No new runner
or dispatcher was introduced to convert them into a sequential Run All workflow.

The original rendering sections read the accepted intermediate result stores.
Rebuilding upstream working tables does not automatically replace every formal
publication archive. Expression matrices, models, accepted tree/result stores,
and other section-specific inputs must be supplied where the original code
requires them. Illustrator layout edits are outside this source collection.

FS4B retains explicit source-identity profiles. Use `original_archive` for the
original frozen result store, and `released_rebuild` for results from the copied
worker source. Editing a worker string, including its local paths, changes its
source hash. Derive the exact expected fingerprints from the reviewed configured
source while retaining the strict assertions.

## Candidate-analysis data package

The candidate-analysis data archive is available in the `v1.0.0` Release:
[fig4_candidate_inputs.zip](https://github.com/MingCen-WestChina/Integration-E9-P28-SGN/releases/download/v1.0.0/fig4_candidate_inputs.zip).
The published asset retains the verified original ZIP bytes.

It contains the original, unchanged files:

- `Processed Data/07_Candidate_regulors/Inputs.h5`
- `Processed Data/04_Embryonic_cluster&Inte/CSS_3000_Embryonic_Inte_adata.h5ad`

Keep both `Processed Data` subfolders when extracting, so the relative count-file
reference inside `Inputs.h5` resolves correctly. Configure the selected
notebook's path placeholders separately. The archive includes English usage
notes and SHA256 checksums. Cell and gene ordering match, and the count-file
checksum matches the provenance stored inside `Inputs.h5`.

This pair supplies the candidate-analysis intermediate inputs used by the
Fig4B–D and related SGN supplementary sections. It contains no P7 data and does
not supply every publication result store or model. `Inputs.h5` contains
previously aligned RNA/count and measurement fields, AUCell values, metadata,
motif-supported gene sets, reviewed annotations, and trajectory correspondence.
Its standalone assembly script was not located in the retained project source.
The downstream candidate algorithms are included.

## Data and software

Public datasets are GSE178931, GSE231513, GSE165502, GSE195500, and GSE202920;
the P28 libraries are GSM6138444 and GSM6138445. P7 is the newly generated study
dataset and must be supplied separately. Human-genetics sections use FinnGen R13
and GTEx v8. SCENIC requires the TF catalog, motif rankings, and annotations
specified in the original source.

Python dependencies are listed in each source block's imports. They include
NumPy, pandas, SciPy, AnnData, Scanpy, scikit-learn, scvi-tools/PyTorch,
SimplePPT, NetworkX, pySCENIC/arboreto/ctxcore, h5py, and Matplotlib. R sections
use packages including Monocle3, Slingshot, mgcv, edgeR, and coloc. Figure QA
source additionally uses Pillow, pypdf, and Poppler. Some original counting and
acquisition programs require Linux command-line tools.

Study versions recorded in the manuscript include Cutadapt 5.1, STAR 2.7.11b,
Monocle3 1.4.27, Slingshot 2.18.0, edgeR 4.8.2, mgcv 1.9.4,
pySCENIC 0.12.1, and arboreto 0.1.6. Workflows used separate environments;
this is not a single environment lockfile.

The organisation was checked against the existing source blocks and notebook
format. No model fitting, statistical analysis, or figure generation was run
while merging these files. The merged source has been reviewed before publication.
