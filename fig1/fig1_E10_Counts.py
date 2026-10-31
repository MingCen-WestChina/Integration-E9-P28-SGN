#!/usr/bin/env python3
# Source: Counts-E10.5.py. Analysis parameters and validation gates are retained.
# Configure every ROOT placeholder independently before running.
# Per-run tool filenames and atomic temporary-file handling are unchanged.
"""Read and diagnose the selected E10.5 raw read counts without filtering or export."""

# Shared paths and category mappings
import sys
sys.dont_write_bytecode = True
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

# Set each storage placeholder below to its own local file or directory.
ROOT = Path()
# Set the file or directory for work.
WORK = ROOT
PLATE_ORDER = [f"GSE231513_plate_{i}" for i in range(1, 7)]
PLATE_COLORS = dict(zip(PLATE_ORDER, ["#6687A5", "#C2896B", "#8C9E75", "#A58DA6", "#B8A263", "#659F9B"]))
ANNOTATION_ORDER = ["O-2", "O-3"]
# Set a writable Matplotlib cache directory.
os.environ["MPLCONFIGDIR"] = str(ROOT)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def md5(path):
    value = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


# Load verified counts, completed run records and the unchanged source object
# Uses counts_unstranded.tsv.gz from fig1_E10_full_counts.py.
COUNTS_FILE = ROOT
# Uses full_counts_results.json from fig1_E10_full_counts.py.
SUMMARY_FILE = ROOT
# Set the GENCODE M25 primary-assembly GTF gzip file.
GTF_FILE = ROOT
# Set the approved 73-cell E10.5 selection H5AD input.
SOURCE_FILE = ROOT
INPUT_MD5 = {COUNTS_FILE: "ca63f58d71b2da398a95628c8d28a468", SUMMARY_FILE: "a539130b067b327aaf1bab0b4da7d9c1",
             GTF_FILE: "c5125258a0a2c5250ddb4c192abbf4e8"}
EXPECTED_SHAPE = (73, 55487)  # Fixed selected cohort and the completed STAR reference
EXPECTED_SUM, EXPECTED_NNZ = 74875349, 545499


def load_inputs():
    """Return counts, var, source_adata, qc, run_summary in unchanged source obs order."""
    for path, expected in INPUT_MD5.items():
        require(path.is_file() and md5(path) == expected, f"Input MD5 differs: {path}")
    run_summary = json.loads(SUMMARY_FILE.read_text())
    require(run_summary["completed_cells"] == run_summary["expected_cells"] == EXPECTED_SHAPE[0]
            and run_summary["reused_pilot_cells"] == 9 and run_summary["completed_new_cells"] == 64,
            "The full-count run is incomplete.")
    require(run_summary["configuration"]["selected_branch"] == "adapter_q20"
            and run_summary["configuration"]["STAR_count_column"] == 2, "Unexpected trimming branch/count column.")
    require(run_summary["matrix"]["md5"] == INPUT_MD5[COUNTS_FILE], "Summary/matrix checksum mismatch.")
    cell_ids, rows = [], []
    with gzip.open(COUNTS_FILE, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        require(header[0] == "cell_id_global", "Unexpected count-matrix header.")
        gene_ids = header[1:]
        for line in handle:
            cell, values = line.rstrip("\n").split("\t", 1)
            require(re.fullmatch(r"[0-9]+(?:\t[0-9]+)*", values), f"Noninteger/negative count: {cell}")
            values = np.fromstring(values, sep="\t", dtype=np.int64)
            require(values.size == len(gene_ids), f"Wrong gene count: {cell}")
            cell_ids.append(cell)
            rows.append(values)
    counts = sparse.csr_matrix(np.stack(rows))
    require(counts.shape == EXPECTED_SHAPE and counts.nnz == EXPECTED_NNZ and int(counts.sum()) == EXPECTED_SUM,
            "Count dimensions, sum or nonzero count differ from the completed run.")
    require(len(set(cell_ids)) == len(cell_ids) and len(set(gene_ids)) == len(gene_ids), "Duplicate matrix identifiers.")
    source_adata = ad.read_h5ad(SOURCE_FILE)
    require(source_adata.obs_names.is_unique and set(source_adata.obs_names) == set(cell_ids), "Original and reconstructed cell sets differ.")
    counts = counts[pd.Index(cell_ids).get_indexer(source_adata.obs_names)].tocsr()
    records = {entry["cell_id_global"]: entry["record"] for entry in run_summary["cells"]}
    require(len(run_summary["cells"]) == len(records) == EXPECTED_SHAPE[0] and set(records) == set(cell_ids),
            "Completed-run cell identities differ.")
    require(source_adata.obs["stagewise_annotation_code"].astype(str).value_counts().to_dict() == {"O-2": 44, "O-3": 29},
            "Original selected annotation counts changed.")
    for cell in source_adata.obs_names:
        record, obs = records[cell], source_adata.obs.loc[cell]
        require(record["cell"] == obs["original_cell_id"] and record["plate"] == obs["plate_id"]
                and record["annotation"] == obs["stagewise_annotation_code"], f"Original metadata mapping differs: {cell}")
    var = read_gene_annotations(gene_ids)
    qc = calculate_qc(counts, var, source_adata, records)
    return counts, var, source_adata, qc, run_summary


# Annotate genes and calculate diagnostic metrics without drop thresholds
MT_CHROMOSOME = "chrM"          # Mitochondrial genes use the verified GENCODE genomic location
HB_SYMBOL_PATTERN = r"^(Hba|Hbb)-"  # Explicit current-project hemoglobin symbol rule
RIBO_SYMBOL_PATTERN = r"^Rp[sl]"     # Requested symbol-prefix diagnostic, not a curated ribosomal gene set


def read_gene_annotations(gene_ids):
    annotations = []
    with gzip.open(GTF_FILE, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            attributes = dict(re.findall(r'(\S+) "([^\"]*)";', fields[8]))
            gene_id = attributes["gene_id"]
            annotations.append({"gene_id_versioned": gene_id, "gene_id_base": re.sub(r"\.\d+$", "", gene_id),
                                "gene_symbol": attributes["gene_name"], "chromosome": fields[0], "gene_type": attributes["gene_type"]})
    var = pd.DataFrame(annotations).set_index("gene_id_versioned", drop=False)
    require(var.index.is_unique and set(var.index) == set(gene_ids), "GTF gene IDs do not exactly match the STAR matrix.")
    var = var.loc[gene_ids].copy()
    var.index.name = "gene_id_versioned_index"
    var["is_mt"] = var["chromosome"].eq(MT_CHROMOSOME)
    var["is_hb"] = var["gene_symbol"].str.match(HB_SYMBOL_PATTERN)
    var["is_ribo_symbol"] = var["gene_symbol"].str.match(RIBO_SYMBOL_PATTERN)
    return var


def calculate_qc(counts, var, source_adata, records):
    qc = pd.DataFrame(index=source_adata.obs_names.copy())
    for key in ("original_cell_id", "plate_id", "stagewise_annotation_code"):
        qc[key] = source_adata.obs[key].copy()
    qc["run_accession"] = [records[cell]["run"] for cell in qc.index]
    qc["original_read_pairs"] = [records[cell]["original_pairs"] for cell in qc.index]
    qc["n_reads"] = 2 * qc["original_read_pairs"]  # Individual raw reads, two reads per original pair
    for field in ("retained_pairs", "unique_mapped_pairs", "multimapped_pairs"):
        qc[field] = [records[cell][field] for cell in qc.index]
    qc["total_assigned_counts"] = np.asarray(counts.sum(axis=1)).ravel()
    qc["detected_genes"] = np.asarray((counts > 0).sum(axis=1)).ravel()
    require((qc["original_read_pairs"] > 0).all() and (qc["total_assigned_counts"] > 0).all(), "Zero input or zero assigned counts.")
    for cell in qc.index:
        recorded = records[cell]["gene_counts_by_column"]["unstranded"]
        require(qc.at[cell, "total_assigned_counts"] == recorded["assigned_pairs"]
                and qc.at[cell, "detected_genes"] == recorded["detected_genes"], f"Matrix/STAR gene-yield mismatch: {cell}")
    for name, numerator in (("retained", "retained_pairs"), ("unique_mapped", "unique_mapped_pairs"), ("assigned", "total_assigned_counts")):
        qc[f"{name}_fraction_original_pairs"] = qc[numerator] / qc["original_read_pairs"]
    qc["assigned_fraction_retained_pairs"] = qc["total_assigned_counts"] / qc["retained_pairs"]
    qc["unique_mapped_fraction_retained_pairs"] = qc["unique_mapped_pairs"] / qc["retained_pairs"]
    for name, mask in (("mt", "is_mt"), ("hb", "is_hb"), ("ribo_symbol", "is_ribo_symbol")):
        require(var[mask].any(), f"Empty reference feature set: {mask}")
        qc[f"pct_{name}"] = 100 * np.asarray(counts[:, var[mask].to_numpy()].sum(axis=1)).ravel() / qc["total_assigned_counts"]
    maximum = counts.max(axis=1).toarray().ravel()
    qc["pct_top_gene"] = 100 * maximum / qc["total_assigned_counts"]
    qc["top_gene_id"] = var.index[np.asarray(counts.argmax(axis=1)).ravel()]
    qc["top_gene_symbol"] = var.loc[qc["top_gene_id"], "gene_symbol"].to_numpy()
    qc["pct_ercc"] = np.nan
    qc["pct_reporter"] = np.nan
    qc.attrs["rules"] = {"matrix_unit": "unstranded read-pair assignments; not UMI counts",
                         "feature_fraction_denominator": "all GENCODE M25 assigned gene counts per cell",
                         "mt": f"GENCODE M25 chromosome == {MT_CHROMOSOME}", "hb": HB_SYMBOL_PATTERN,
                         "ribo_symbol": RIBO_SYMBOL_PATTERN + "; requested diagnostic prefix rule, may include pseudogenes",
                         "ERCC_reporter": "Not quantified: controls were absent from the counting reference; NaN is not zero",
                         "cell_filtering": "No thresholds applied; all 73 selected cells retained"}
    require(set(qc["plate_id"].astype(str)) == set(PLATE_ORDER), "Expected all six source plates.")
    return qc


# Render four complementary diagnostic panels
QC_FIGSIZE = (9.4, 6.5)          # Figure size in inches; affects this QC plot only
QC_FIGURE_DPI = 1201            # Required raster export resolution
POINT_SIZE = 20
POINT_ALPHA = 0.80
FONT_SIZE = 8
JITTER_WIDTH = 0.17
JITTER_SEED = 42             # Only the horizontal display jitter is randomized
# Set the file or directory for qc figure file.
QC_FIGURE_FILE = ROOT


def plot_qc(qc, output=QC_FIGURE_FILE):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    rng = np.random.default_rng(JITTER_SEED)
    with plt.rc_context({"font.size": FONT_SIZE, "axes.spines.top": False, "axes.spines.right": False}):
        fig, axes = plt.subplots(2, 2, figsize=QC_FIGSIZE)
        for plate in PLATE_ORDER:
            frame = qc.loc[qc["plate_id"].astype(str).eq(plate)]
            axes[0, 0].scatter(frame["total_assigned_counts"], frame["detected_genes"], s=POINT_SIZE,
                               color=PLATE_COLORS[plate], alpha=POINT_ALPHA, edgecolors="white", linewidths=0.35)
        axes[0, 0].set(xscale="log", xlabel="Assigned read-pair counts (log scale)", ylabel="Detected genes", title="A  Library complexity")
        groups = [qc.loc[qc["plate_id"].astype(str).eq(plate)] for plate in PLATE_ORDER]
        for ax, metric, scale, title, label in ((axes[0, 1], "pct_mt", 1, "B  Mitochondrial proportion", "Mitochondrial gene counts (%)"),
                                               (axes[1, 0], "assigned_fraction_original_pairs", 100, "C  Gene assignment yield", "Assigned / original read pairs (%)")):
            ax.boxplot([frame[metric].to_numpy() * scale for frame in groups], positions=np.arange(1, 7), widths=0.50,
                       showfliers=False, medianprops={"color": "0.2"}, boxprops={"color": "0.55"},
                       whiskerprops={"color": "0.55"}, capprops={"color": "0.55"})
            for i, (plate, frame) in enumerate(zip(PLATE_ORDER, groups), 1):
                ax.scatter(i + rng.uniform(-JITTER_WIDTH, JITTER_WIDTH, len(frame)), frame[metric] * scale,
                           s=POINT_SIZE, color=PLATE_COLORS[plate], alpha=POINT_ALPHA, edgecolors="white", linewidths=0.35)
            ax.set(xticks=np.arange(1, 7), xticklabels=[f"P{i}\n(n={len(frame)})" for i, frame in enumerate(groups, 1)],
                   xlabel="Source plate", ylabel=label, title=title)
        ax = axes[1, 1]
        groups = [qc.loc[qc["stagewise_annotation_code"].astype(str).eq(annotation)] for annotation in ANNOTATION_ORDER]
        ax.boxplot([frame["pct_top_gene"] for frame in groups], positions=[1, 2], widths=0.40, showfliers=False,
                   medianprops={"color": "0.2"}, boxprops={"color": "0.55"}, whiskerprops={"color": "0.55"}, capprops={"color": "0.55"})
        for i, frame in enumerate(groups, 1):
            ax.scatter(i + rng.uniform(-JITTER_WIDTH, JITTER_WIDTH, len(frame)), frame["pct_top_gene"], s=POINT_SIZE,
                       c=[PLATE_COLORS[str(p)] for p in frame["plate_id"]], alpha=POINT_ALPHA, edgecolors="white", linewidths=0.35)
        ax.set(xticks=[1, 2], xticklabels=[f"{name}\n(n={len(frame)})" for name, frame in zip(ANNOTATION_ORDER, groups)],
               ylabel="Counts in the most abundant gene (%)", xlabel="Existing stagewise annotation", title="D  Dominance of the top gene")
        handles = [Line2D([], [], linestyle="", marker="o", color=PLATE_COLORS[plate], markersize=5, label=f"Plate {i}")
                   for i, plate in enumerate(PLATE_ORDER, 1)]
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.946), ncol=6, frameon=False)
        fig.suptitle("E10.5 reconstructed read-count QC | 73 selected cells", y=0.985, fontsize=10)
        fig.subplots_adjust(left=0.09, right=0.985, bottom=0.105, top=0.85, wspace=0.34, hspace=0.46)
        fig.savefig(output, dpi=QC_FIGURE_DPI, facecolor="white")
        plt.close(fig)
    return Path(output)



# Shared inputs and confirmed marker groups
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

MARKER_GROUPS = {
    "Otic": ["Six1", "Eya1", "Foxg1", "Fgf10"],
    "Neuroblast": ["Neurod1", "Insm1", "Mki67", "Cdk1", "Ccnb1"],
    "Neuronal": ["Tubb3", "Dcx", "Myt1", "Ncam1", "Nefl", "Nefm", "Gap43", "Pou4f1", "Ebf1", "St18", "Zeb1", "Zeb2", "Isl1"],
    "Alternative": ["Phox2a", "Phox2b", "Tlx2", "Dbh", "Slc18a3", "Crabp1", "Meis1", "Dcc", "Pou3f2", "Hoxb2"],
    "Blood": ["Gata1", "Klf1", "Alas2", "Ptprc", "Csf1r", "Spi1", "Lst1", "Tyrobp"],
    "Megakaryocyte": ["Pf4", "Ppbp", "Itga2b"],
    "Vascular": ["Kdr", "Cdh5", "Pecam1"],
}
HIGHLIGHT_CELLS = ["E3.6", "G9.4"]


# Normalize raw counts for display only and sort the plotted cells
NORMALIZATION_TOTAL = 10000   # Denominator includes every quantified GENCODE gene


def marker_values():
    counts, var, source_adata, qc, run_summary = load_inputs()
    requested = [gene for genes in MARKER_GROUPS.values() for gene in genes]
    require(len(requested) == len(set(requested)), "A marker occurs in more than one requested group.")
    symbol_counts = var["gene_symbol"].value_counts()
    duplicates = {gene: int(symbol_counts[gene]) for gene in requested if symbol_counts.get(gene, 0) > 1}
    require(not duplicates, f"Ambiguous marker symbols; explicit gene-ID selection is required: {duplicates}")
    missing = [gene for gene in requested if symbol_counts.get(gene, 0) == 0]
    print("Missing markers: " + (", ".join(missing) if missing else "none"), flush=True)
    groups = {group: [gene for gene in genes if gene not in missing] for group, genes in MARKER_GROUPS.items()}
    require(all(groups.values()), "At least one complete marker group is absent from the reference.")
    cells = qc[["original_cell_id", "plate_id", "stagewise_annotation_code"]].copy()
    cells["annotation_order"] = cells["stagewise_annotation_code"].astype(str).map(dict(zip(ANNOTATION_ORDER, range(len(ANNOTATION_ORDER)))))
    cells["plate_order"] = cells["plate_id"].astype(str).map(dict(zip(PLATE_ORDER, range(len(PLATE_ORDER)))))
    cells = cells.sort_values(["annotation_order", "plate_order", "original_cell_id"], kind="stable")
    require(cells[["annotation_order", "plate_order"]].notna().all().all(), "Unknown annotation or source plate.")
    cell_order = source_adata.obs_names.get_indexer(cells.index)
    gene_order = [int(np.flatnonzero(var["gene_symbol"].to_numpy() == gene)[0]) for genes in groups.values() for gene in genes]
    totals = np.asarray(counts.sum(axis=1)).ravel()
    require((totals > 0).all(), "Cannot normalize a zero-count cell.")
    values = counts[cell_order][:, gene_order].toarray().astype(np.float64)
    values = np.log1p(values * (NORMALIZATION_TOTAL / totals[cell_order])[:, None]).T
    require(np.isfinite(values).all() and (values >= 0).all(), "Invalid transformed marker values.")
    require(all((cells["original_cell_id"].astype(str) == cell).sum() == 1 for cell in HIGHLIGHT_CELLS), "A highlighted cell is absent or duplicated.")
    return values, groups, cells, missing


# Plot all cells with one shared expression scale
MARKER_FIGSIZE = (9.2, 8.8)          # Inches; increased height keeps the 46 gene rows readable
MARKER_FIGURE_DPI = 1201
GENE_FONT_SIZE = 6.5
CELL_FONT_SIZE = 5.0
GROUP_FONT_SIZE = 6.0
COLORMAP = "YlOrRd"
ANNOTATION_COLORS = ["#BCCDD8", "#D7B89C"]
HIGHLIGHT_COLOR = "#7B2C35"
# Set the file or directory for marker figure file.
MARKER_FIGURE_FILE = ROOT


def plot_marker_qc(values, groups, cells):
    gene_labels, gene_positions, group_bounds, blocks = [], [], [], []
    source_start, plot_start = 0, 0
    for group, genes in groups.items():
        blocks.append(values[source_start:source_start + len(genes)])
        gene_labels.extend(genes)
        gene_positions.extend(range(plot_start, plot_start + len(genes)))
        group_bounds.append((group, plot_start, plot_start + len(genes) - 1))
        source_start += len(genes)
        plot_start += len(genes)
        if source_start < values.shape[0]:
            blocks.append(np.full((1, values.shape[1]), np.nan))
            plot_start += 1
    heatmap = np.vstack(blocks)
    cmap = plt.get_cmap(COLORMAP).copy()
    cmap.set_bad("white")
    fig = plt.figure(figsize=MARKER_FIGSIZE)
    main = fig.add_axes([0.185, 0.15, 0.755, 0.675])
    annotation_ax = fig.add_axes([0.185, 0.865, 0.755, 0.022])
    plate_ax = fig.add_axes([0.185, 0.832, 0.755, 0.023])
    color_ax = fig.add_axes([0.54, 0.927, 0.30, 0.010])
    image = main.imshow(heatmap, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=float(values.max()))
    main.set_yticks(gene_positions, gene_labels, fontsize=GENE_FONT_SIZE)
    main.set_xticks(range(len(cells)), cells["original_cell_id"].astype(str), rotation=90, fontsize=CELL_FONT_SIZE)
    main.tick_params(axis="both", length=0, pad=2)
    main.set_xlabel("Individual selected cells, ordered by existing annotation, plate and cell ID", fontsize=7, labelpad=8)
    for group, start, end in group_bounds:
        middle = (start + end) / 2
        main.plot([-0.112, -0.103], [start - 0.4, start - 0.4], transform=main.get_yaxis_transform(), color="0.35", lw=0.6, clip_on=False)
        main.plot([-0.112, -0.103], [end + 0.4, end + 0.4], transform=main.get_yaxis_transform(), color="0.35", lw=0.6, clip_on=False)
        main.plot([-0.112, -0.112], [start - 0.4, end + 0.4], transform=main.get_yaxis_transform(), color="0.35", lw=0.6, clip_on=False)
        main.text(-0.18, middle, group, transform=main.get_yaxis_transform(), ha="center", va="center",
                  fontsize=GROUP_FONT_SIZE, color="0.25", clip_on=False)
    annotation_ax.imshow([cells["annotation_order"].to_numpy()], aspect="auto", interpolation="nearest",
                         cmap=ListedColormap(ANNOTATION_COLORS), vmin=0, vmax=1)
    plate_ax.imshow([cells["plate_order"].to_numpy()], aspect="auto", interpolation="nearest",
                   cmap=ListedColormap([PLATE_COLORS[p] for p in PLATE_ORDER]), vmin=0, vmax=5)
    for ax, label in ((annotation_ax, "Annotation"), (plate_ax, "Plate")):
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(-0.012, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=6.5)
        for spine in ax.spines.values():
            spine.set_visible(False)
    annotation_values = cells["stagewise_annotation_code"].astype(str).to_numpy()
    for annotation in ANNOTATION_ORDER:
        positions = np.flatnonzero(annotation_values == annotation)
        annotation_ax.text(float(positions.mean()), 0, f"{annotation} (n={len(positions)})", ha="center", va="center", fontsize=7)
    block_start = 0
    plate_values = cells["plate_order"].to_numpy()
    for boundary in range(1, len(cells) + 1):
        if boundary == len(cells) or plate_values[boundary] != plate_values[block_start] or annotation_values[boundary] != annotation_values[block_start]:
            plate_ax.text((block_start + boundary - 1) / 2, 0, str(int(plate_values[block_start]) + 1), ha="center", va="center", fontsize=5.5)
            if boundary < len(cells):
                width = 1.1 if annotation_values[boundary] != annotation_values[block_start] else 0.30
                main.axvline(boundary - 0.5, color="0.38", linewidth=width)
                plate_ax.axvline(boundary - 0.5, color="white", linewidth=0.7)
            block_start = boundary
    for cell in HIGHLIGHT_CELLS:
        position = int(np.flatnonzero(cells["original_cell_id"].astype(str).to_numpy() == cell)[0])
        main.add_patch(Rectangle((position - 0.5, -0.5), 1, heatmap.shape[0], fill=False, linewidth=0.7, edgecolor=HIGHLIGHT_COLOR))
        main.get_xticklabels()[position].set_color(HIGHLIGHT_COLOR)
        main.get_xticklabels()[position].set_fontweight("bold")
        annotation_ax.annotate(cell, xy=(position, 0.5), xycoords=("data", "axes fraction"), xytext=(position, 1.5),
                               textcoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=6.5, color=HIGHLIGHT_COLOR,
                               arrowprops={"arrowstyle": "-", "color": HIGHLIGHT_COLOR, "lw": 0.7}, annotation_clip=False)
        print(f"Highlighted cell: {cell}; plotted column={position + 1}; annotation={annotation_values[position]}; plate={plate_values[position] + 1}", flush=True)
    colorbar = fig.colorbar(image, cax=color_ax, orientation="horizontal")
    colorbar.ax.tick_params(labelsize=6, length=2)
    colorbar.set_label("log1p(read-pair counts per 10,000)", fontsize=6.5, labelpad=2)
    fig.text(0.185, 0.957, "E10.5 cell-level marker QC", fontsize=11, ha="left")
    fig.text(0.185, 0.934, "Shared scale across genes; all 73 selected cells", fontsize=7, ha="left")
    fig.savefig(MARKER_FIGURE_FILE, dpi=MARKER_FIGURE_DPI, facecolor="white")
    plt.close(fig)
    return MARKER_FIGURE_FILE



# Shared paths and verified cohort constants
import argparse
from datetime import datetime, timezone
import importlib.metadata
import platform
import uuid

# Set the file or directory for output.
OUTPUT = ROOT
# Set the H5AD export-contract text file used for provenance hashing.
SCHEMA_FILE = ROOT
SCHEMA_VERSION = "2.0.0"
CONTRACT_KEYS = ("matrix_contract", "selection_contract", "integration_contract", "source_provenance", "export_contract")
EXPORT_EXPECTED_SHAPE = (73, 55438)
EXPORT_EXPECTED_SUM, EXPORT_EXPECTED_NNZ, EXPORT_EXPECTED_ZERO_GENES = 72854613, 544404, 30930
FEATURE_BLOCKER = "versioned_cross_dataset_gene_harmonization_not_run"


def json_ready(value):
    """Convert historical numpy/pandas metadata into portable strict JSON."""
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray, pd.Index)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if value is None or value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return value


def json_text(value):
    return json.dumps(json_ready(value), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ids_sha256(values):
    return hashlib.sha256(json_text([str(v) for v in values]).encode("utf-8")).hexdigest()


def frame_sha256(frame):
    payload = {"index": frame.index.tolist(), "columns": frame.columns.tolist(),
               "values_by_column": [frame[key].tolist() for key in frame.columns]}
    return hashlib.sha256(json_text(payload).encode("utf-8")).hexdigest()


def matrix_sha256(matrix):
    require(sparse.isspmatrix_csr(matrix) and matrix.dtype == np.dtype("int32"), "Expected int32 CSR for hashing.")
    value = matrix.copy()
    value.sum_duplicates()
    value.sort_indices()
    value.eliminate_zeros()
    digest = hashlib.sha256(b"CSR-int32-v1\0")
    for array, dtype in ((value.shape, "<i8"), (value.indptr, "<i8"),
                         (value.indices, "<i8"), (value.data, "<i4")):
        digest.update(np.asarray(array, dtype=dtype).tobytes(order="C"))
    return digest.hexdigest()


# Preserve source metadata and replace obsolete matrix/QC state
SOURCE_QC_FIELDS = """
n_nonzero_source_features row_sum_source_values total_counts total_counts_source
n_features_source total_counts_endogenous n_genes_endogenous pct_ercc ercc_available
pct_reporter reporter_available pct_mito mito_available pct_ribo ribo_available
log10_counts_endogenous log10_genes_endogenous zero_endogenous zero_source technical_only
low_counts low_genes high_mito high_ercc high_counts_review high_genes_review candidate_3mad
high_complexity_review n_core_rules_evaluated qc_candidate_reason n_positive_eligible_features
mean_source_value_eligible zero_eligible_expression n_eligible_features author_n_genes
author_total_count author_pct_mt author_pct_ribosomal log10_author_total_count log10_author_n_genes
external_spike_counts fail_zero_expression fail_low_complexity fail_high_complexity fail_high_mito
fail_high_ercc fail_missing_qc_metric qc_fail_reason qc_fail_primary_reason
qc_pass_before_spike_correction qc_pass_changed_by_spike_correction joint_low_review
joint_high_review single_metric_review review_candidate source_author_qc_status source_object_key
""".split()
PRIMARY_REVIEW = {"E3.6": "low_unique_mapping_fraction_original_pairs; elevated_mitochondrial_fraction",
                  "G9.4": "low_unique_mapping_fraction_original_pairs"}
SECONDARY_REVIEW = {"F11.4": "lowest_assigned_counts_with_preserved_detected_genes",
                    "A2.3": "lowest_detected_genes_with_preserved_assigned_counts",
                    "E10.6": "highest_detected_genes; elevated_no_feature_fraction_original_pairs"}


def build_obs(source, qc, final_counts):
    obs = source.obs.copy(deep=True)
    migration = {}

    def preserve(key):
        if key not in obs:
            return
        historical = "historical_" + key
        require(historical not in obs, f"Historical column collision: {historical}")
        obs[historical] = obs[key].copy()
        obs.drop(columns=[key], inplace=True)
        migration[key] = historical

    def current(key, value):
        nonlocal obs
        preserve(key)
        obs = obs.copy()
        obs[key] = value

    for key in SOURCE_QC_FIELDS + [k for k in obs if k.startswith("before_spike_correction__")]:
        preserve(key)
    # Identity columns in the QC table must match, rather than silently overwrite.
    for key in qc:
        if key in ("original_cell_id", "plate_id", "stagewise_annotation_code"):
            require(obs[key].astype(str).equals(qc[key].astype(str)), f"QC metadata mismatch: {key}")
        else:
            current(key, qc[key].copy())
    updates = {
        "source_file": str(COUNTS_FILE.relative_to(ROOT)),
        "metadata_source": str(SOURCE_FILE.relative_to(ROOT)),
        "metadata_source_role": "accepted_stagewise_annotations_and_cell_identity",
        "technology": "plate_based_full_length_scRNA_seq",
        "lineage_status": "selected_otic_CVG_neuroblast_or_developing_neuron",
        "cell_scope_status": "accepted_73_stagewise_cells_with_reconstructed_counts",
        "matrix_value_type": "raw_integer_counts", "matrix_class": "raw_read_counts",
        "matrix_unit": "unstranded_read_pair_assignments", "raw_counts_available": True,
        "raw_counts_status": "available_verified", "source_normalization": "not_applied",
        "source_log_transform": "not_applied", "self_qc_mode": "reconstructed_read_count_diagnostics",
        "local_entry_rule": "unchanged_73_cell_whitelist_from_accepted_stagewise_object",
        "current_pipeline_qc_applied": True, "current_pipeline_normalization_applied": False,
        "raw_count_model_permission": "verified_counts_require_cross_dataset_gene_harmonization",
        "droplet_background_model_permission": "not_applicable_plate_based",
        "eligible_for_count_qc": True, "eligible_for_processed_audit": False,
        "integer_like": True, "count_status": True,
        "qc_grouping_status": "six_source_plates_for_read_count_diagnostics",
        "qc_rule_set": "diagnostic_review_without_automatic_cell_thresholds",
        "qc_evidence_scope": "full_55487_GENCODE_M25_gene_read_counts_before_mt_hb_removal",
        "count_model_compatible": True, "count_based_integration_allowed": False,
        "integration_data_role": "stagewise_pre_harmonization_read_counts",
        "gene_harmonization_status": "not_run", "model_ready": False,
        "scanvi_supervision_status": "unlabeled", "doublet_detection_assessed": False,
        "doublet_status": "not_assessed_on_reconstructed_counts", "doublet_score": np.nan,
        "count_qc_retained": True, "count_qc_decision": "retained_after_count_qc",
        "blood_lineage_status": "retained_non_blood; count_marker_review_no_blood_combination_evidence",
        "blood_qc_pass": True,
    }
    for key, value in updates.items():
        current(key, value)
    for key in ("predicted_doublet", "doublet_qc_pass", "doublet_review", "qc_pass"):
        current(key, pd.array([pd.NA] * len(obs), dtype="boolean"))
    current("technical_qc_pass", pd.array([True] * len(obs), dtype="boolean"))
    originals = obs["original_cell_id"].astype(str)
    require(set(PRIMARY_REVIEW | SECONDARY_REVIEW).issubset(set(originals)), "A reviewed cell is missing.")
    obs["technical_review"] = originals.isin(PRIMARY_REVIEW)
    obs["secondary_technical_review"] = originals.isin(SECONDARY_REVIEW)
    obs["technical_review_reason"] = originals.map(PRIMARY_REVIEW).fillna("")
    obs["secondary_technical_review_reason"] = originals.map(SECONDARY_REVIEW).fillna("")
    obs.loc[obs["technical_review"], "technical_qc_pass"] = pd.NA
    obs.loc[obs["technical_review"], "count_qc_decision"] = "retained_with_technical_review"
    current("cell_type", pd.Categorical(obs["stagewise_annotation"].astype(str)))
    current("total_counts_analysis", np.asarray(final_counts.sum(axis=1)).ravel().astype(np.int64))
    current("n_genes_analysis", np.asarray((final_counts > 0).sum(axis=1)).ravel().astype(np.int64))
    # Full-reference QC is kept separately from statistics on the exported gene axis.
    current("total_counts_source", qc["total_assigned_counts"].to_numpy())
    current("n_features_source", qc["detected_genes"].to_numpy())
    current("pct_mito", qc["pct_mt"].to_numpy())
    obs["mito_available"] = True
    obs["ercc_available"] = False
    obs["reporter_available"] = False
    # Prefix-based ribosomal QC is named explicitly; a shared ribosomal set is pending.
    obs["pct_ribo"] = np.nan
    obs["ribo_available"] = False
    obs["qc_limitation"] = "doublets_not_assessed; ERCC_and_reporter_absent_from_reference; no_automatic_cell_thresholds"
    return obs, migration


# Retain the complete biological gene axis, including zero-count genes
EXCLUDED_MT_COUNT, EXCLUDED_HB_COUNT = 37, 12


def build_var(source_var, counts):
    require(int(source_var["is_mt"].sum()) == EXCLUDED_MT_COUNT
            and int(source_var["is_hb"].sum()) == EXCLUDED_HB_COUNT
            and not (source_var["is_mt"] & source_var["is_hb"]).any(), "Unexpected excluded feature identities.")
    excluded = source_var["is_mt"] | source_var["is_hb"]
    var = source_var.loc[~excluded].copy()
    final_counts = counts[:, ~excluded.to_numpy()].tocsr()
    require(final_counts.data.size and np.isfinite(final_counts.data).all()
            and final_counts.data.min() >= 0 and final_counts.data.max() <= np.iinfo(np.int32).max
            and np.issubdtype(final_counts.dtype, np.integer), "Counts cannot be safely represented as int32.")
    final_counts = final_counts.astype(np.int32)
    final_counts.sum_duplicates()
    final_counts.sort_indices()
    final_counts.eliminate_zeros()
    var["gene_id_analysis"] = var.index.astype(str)
    var["gene_id_type"] = "ensembl_gene_id_versioned"
    var["source_feature_id"] = var["gene_id_versioned"].copy()
    var["gene_id_original"] = var["gene_id_versioned"].copy()
    var["source_gene_id_type"] = "ensembl_gene_id_versioned"
    var["feature_row_id"] = np.flatnonzero(~excluded.to_numpy()) + 1
    var["feature_type"] = "Gene Expression"
    var["gene_id_repair"] = "none"
    var["integration_gene_id_candidate"] = var["gene_id_base"].copy()
    var["integration_gene_id_candidate_type"] = "ensembl_gene_id"
    var["canonical_gene_id"] = pd.Categorical([None] * len(var), categories=[])
    var["canonical_mapping_status"] = "not_run"
    var["eligible_for_cross_dataset_integration"] = False
    var["cross_dataset_feature_blocker"] = FEATURE_BLOCKER
    var["eligible_for_within_dataset_analysis"] = True
    var["use_for_model_input"] = True
    for key in ("mt_feature", "hemoglobin_feature", "ERCC_feature", "reporter_or_vector_feature", "technical_feature"):
        var[key] = False
    var["source_ribosomal_feature"] = var["is_ribo_symbol"].copy()
    var["ribosomal_feature"] = pd.array([pd.NA] * len(var), dtype="boolean")
    var["ribosomal_classification_status"] = "pending_shared_gene_set"
    var["n_cells_in_export"] = np.asarray((final_counts > 0).sum(axis=0)).ravel().astype(np.int64)
    var["total_counts_in_export"] = np.asarray(final_counts.sum(axis=0)).ravel().astype(np.int64)
    var["detected_in_exported_cells"] = var["total_counts_in_export"] > 0
    excluded_ids = {"mitochondrial": source_var.index[source_var["is_mt"]].tolist(),
                    "hemoglobin": source_var.index[source_var["is_hb"]].tolist()}
    return final_counts, var, excluded_ids


# Construct the five current contracts and immutable source fingerprints
def build_object(output=OUTPUT):
    counts, source_var, source, qc, run_summary = load_inputs()
    final_counts, var, excluded_ids = build_var(source_var, counts)
    obs, migration = build_obs(source, qc, final_counts)
    result = ad.AnnData(X=final_counts.copy(), obs=obs, var=var)
    result.layers["counts"] = final_counts.copy()
    selection = source.uns["selection_contract"]
    source_notebook = str(selection["source_notebook"])
    source_object = str(selection["final_selection_object"])
    # Set the optional reference_manifest.json provenance file.
    reference_manifest = ROOT
    source_files = {"raw_read_counts": COUNTS_FILE, "completed_count_run": SUMMARY_FILE,
                    "reference_gene_annotation": GTF_FILE, "accepted_cell_selection": SOURCE_FILE,
                    "qc_code": Path(__file__), "export_code": Path(__file__), "schema_rules": SCHEMA_FILE}
    if reference_manifest.is_file():
        source_files["reference_manifest"] = reference_manifest
    source_paths = {role: str(path.resolve()) for role, path in source_files.items()}
    source_hashes = {role: file_sha256(path) for role, path in source_files.items()}
    legacy = {key + "_json": json_text(source.uns.get(key, {})) for key in CONTRACT_KEYS}
    result.uns["matrix_contract"] = {
        "schema_version": SCHEMA_VERSION, "raw_counts_available": True, "raw_counts_status": "available_verified",
        "counts_layer": "counts", "counts_layer_semantics": "raw_counts", "count_type": "read",
        "count_unit_detail": "unstranded paired-end read-pair assignments from STAR GeneCounts column 2; not UMI",
        "counts_verification_source": "Normal SRA FASTQ -> adapter_q20 -> STAR 2.7.11b GRCm38/GENCODE M25; verified run records, MD5, source cell IDs and exact GTF gene IDs",
        "counts_dtype": "int32", "counts_sparse_format": "CSR", "X_semantics": "raw_integer_read_counts",
        "X_dtype": "int32", "X_sparse_format": "CSR", "X_normalization_status": "not_applied",
        "X_log1p_status": "not_applied", "X_scaling_status": "not_applied", "X_counts_relationship": "exact_copy",
        "raw_present": False,
    }
    result.uns["selection_contract"] = {
        "schema_version": SCHEMA_VERSION, "source_notebook": source_notebook, "source_object": source_object,
        "selection_rule": "Preserve accepted stagewise O-2/O-3 cell IDs and order from the source H5AD; no cells removed during count QC/export",
        "n_selected_cells": result.n_obs, "cell_type_counts": obs["stagewise_annotation_code"].astype(str).value_counts().to_dict(),
        "cell_type_label_mapping": obs[["stagewise_annotation_code", "cell_type"]].astype(str).drop_duplicates().set_index("stagewise_annotation_code")["cell_type"].to_dict(),
        "annotation_code_counts": obs["stagewise_annotation_code"].astype(str).value_counts().to_dict(),
        "annotation_status": "stagewise_target_selected_cross_dataset_annotation_pending",
        "annotation_source": str(SOURCE_FILE), "label_policy": "unlabeled_candidates",
        "confirmed_blood_cells_removed_upstream": int(selection["confirmed_blood_cells_removed_upstream"]),
        "cells_removed_during_count_qc": 0, "primary_technical_review_original_ids": list(PRIMARY_REVIEW),
        "secondary_technical_review_original_ids": list(SECONDARY_REVIEW),
        "identity_limitation": "O-3 is a developing CVG neuron candidate and is not SGN-specific",
    }
    result.uns["integration_contract"] = {
        "schema_name": "auditory_neuron_stagewise", "schema_version": SCHEMA_VERSION,
        "input_role": "stagewise_pre_harmonization", "counts_layer": "counts", "batch_key": "integration_batch",
        "labels_key": "scanvi_label", "unlabeled_category": "Unknown", "label_policy": "unlabeled_candidates",
        "gene_harmonization_status": "not_run", "cross_dataset_blocker": FEATURE_BLOCKER,
        "batch_stage_confounded": True, "count_model_compatible": True, "model_ready": False,
        "model_blockers": [FEATURE_BLOCKER, "cross_stage_supervision_labels_not_constructed"],
        "count_model_compatibility_scope": "Verified nonnegative integer read counts for the authorized count-model workflow; full-length read counts are not UMI counts",
        "count_based_integration_allowed_semantics": "False until the joint gene axis and model inputs are ready; does not deny that verified counts are available",
        "batch_identity_scope": "Six source plates; plate identity does not establish animal or biological replicate identity",
    }
    metric_definitions = dict(qc.attrs["rules"])
    metric_definitions.update({
        "total_counts_analysis_n_genes_analysis": "Retained 55438-gene model axis, after 37 chrM and 12 Hba/Hbb genes were excluded",
        "total_counts_source_n_features_source": "Full 55487-gene GENCODE M25 count axis, before feature exclusion",
        "pct_mito": "Alias of pct_mt; full-reference count percentage before feature exclusion",
        "pct_ribo": "Unavailable pending a shared ribosomal gene set; pct_ribo_symbol retains the separate prefix diagnostic",
        "n_reads": "Two individual reads per original paired-end fragment; before trimming",
        "technical_qc_pass": "True for 71 retained cells after diagnostic review; NA for E3.6/G9.4 pending technical review; no automatic cutoffs applied",
        "qc_pass": "NA: no single omnibus pass/fail claim; doublet assessment is unavailable",
        "feature_row_id": "One-based gene position in the source count matrix header, excluding its cell_id_global field",
    })
    result.uns["source_provenance"] = {
        "schema_version": SCHEMA_VERSION, "export_date": datetime.now(timezone.utc).isoformat(),
        "source_files": source_paths, "source_file_sha256": source_hashes,
        "source_notebook": source_notebook, "source_dataset": "GSE231513",
        "source_library": list(PLATE_ORDER), "source_object": source_object,
        "source_reference_assembly": "GRCm38 primary assembly", "source_annotation_release": "GENCODE M25 primary comprehensive",
        "source_reference_evidence": str(reference_manifest) if reference_manifest.is_file() else str(GTF_FILE),
        "count_run_configuration_json": json_text(run_summary["configuration"]),
        "software_versions": {"python": platform.python_version(), **{key: importlib.metadata.version(key)
                              for key in ("anndata", "numpy", "pandas", "scipy", "h5py")}},
        "qc_metric_definitions": metric_definitions,
        "qc_decision": "73 retained, zero excluded; primary and secondary flags are review indicators, not failures",
        "marker_review": "Independent count-marker review found Neurod1 and Insm1 detected in all 73 cells and no combined blood-lineage marker evidence; original O-2/O-3 labels retained",
        "limitations": ["No doublet scoring/classification on reconstructed counts",
                        "ERCC and reporter were absent from the counting reference and their fractions are unavailable",
                        "O-3 is not SGN-specific", "Cross-dataset gene harmonization and supervision remain pending"],
        "migration_history": {"obs_column_mapping_json": json_text(migration), "legacy_contracts": legacy,
                              "original_obs_dtypes_json": json_text({k: str(v) for k, v in source.obs.dtypes.items()}),
                              "original_label_policy": "All 73 scanvi_label values Unknown; unchanged",
                              "source_matrix": "Original continuous-expression H5AD preserved unchanged; no old feature axis or embeddings copied"},
    }
    result.uns["export_contract"] = {
        "schema_version": SCHEMA_VERSION, "output_filename": Path(output).name,
        "n_obs": result.n_obs, "n_vars": result.n_vars,
        "saved_slots": ["X", "obs", "var", "layers/counts", "uns/five_contracts"],
        "excluded_feature_classes": ["mitochondrial", "hemoglobin"], "excluded_feature_ids": excluded_ids,
        "ordered_cell_id_sha256": ids_sha256(result.obs_names), "ordered_feature_id_sha256": ids_sha256(result.var_names),
        "X_sha256": matrix_sha256(result.X), "counts_sha256": matrix_sha256(result.layers["counts"]),
        "obs_sha256": frame_sha256(result.obs), "var_sha256": frame_sha256(result.var),
        "X_nnz": result.X.nnz, "counts_nnz": result.layers["counts"].nnz,
        "total_counts": int(result.layers["counts"].sum()), "zero_count_genes_retained": int((var["total_counts_in_export"] == 0).sum()),
        "full_source_gene_count": source_var.shape[0], "full_source_total_counts": int(counts.sum()),
        "fingerprint_encoding": "SHA256; IDs: UTF-8 strict compact JSON ordered string array; matrices: CSR-int32-v1 + NUL, canonical CSR with duplicates summed, indices sorted and zeros removed, shape/indptr/indices signed int64 little-endian followed by data signed int32 little-endian; frames: compact sorted-key JSON index/columns/column values with missing values null",
        "write_read_validation_policy": "Publish only after temporary H5AD readback passes all fingerprints, contracts, matrix/axis and selected metadata checks",
    }
    validate_object(result, source.obs_names)
    return result


# Validate exact count, metadata and contract preservation
def validate_object(adata, expected_cell_ids=None):
    require(adata.shape == EXPORT_EXPECTED_SHAPE, f"Unexpected export shape: {adata.shape}")
    require(adata.obs_names.is_unique and adata.var_names.is_unique, "Duplicate axis IDs.")
    require(adata.obs_names.astype(str).tolist() == adata.obs["cell_id_global"].astype(str).tolist(), "Global cell IDs differ from the index.")
    if expected_cell_ids is not None:
        require(adata.obs_names.tolist() == list(expected_cell_ids), "Original cell order changed.")
    require(adata.obs["stagewise_annotation_code"].astype(str).value_counts().to_dict() == {"O-2": 44, "O-3": 29}, "Annotation counts changed.")
    require(adata.obs["scanvi_label"].astype(str).eq("Unknown").all(), "Supervision labels changed.")
    require(adata.obs["selected_for_integration"].all() and adata.obs["count_qc_retained"].all(), "A selected cell was removed.")
    require(adata.obs["integration_batch"].astype(str).equals(adata.obs["batch_id_global"].astype(str)), "Batch IDs differ.")
    require(set(adata.obs.loc[adata.obs["technical_review"], "original_cell_id"].astype(str)) == set(PRIMARY_REVIEW), "Technical review IDs changed.")
    for field in ("predicted_doublet", "doublet_qc_pass", "doublet_review"):
        require(str(adata.obs[field].dtype) == "boolean" and adata.obs[field].isna().all(), f"Unassessed doublet field must be nullable NA: {field}")
    require(adata.raw is None and set(adata.layers) == {"counts"}, "Unexpected raw or layers.")
    require(not any(len(getattr(adata, key)) for key in ("obsm", "obsp", "varm", "varp")), "Stale embeddings/graphs present.")
    require(set(adata.uns) == set(CONTRACT_KEYS), "Unexpected uns contracts/state.")
    for matrix in (adata.X, adata.layers["counts"]):
        require(sparse.isspmatrix_csr(matrix) and matrix.dtype == np.dtype("int32"), "Counts/X must be int32 CSR.")
        require(np.isfinite(matrix.data).all() and (matrix.data >= 0).all(), "Invalid counts.")
        require(matrix.nnz == EXPORT_EXPECTED_NNZ and int(matrix.sum()) == EXPORT_EXPECTED_SUM, "Count sum/nnz differs from independent verification.")
    require((adata.X != adata.layers["counts"]).nnz == 0, "X and counts differ.")
    require(not adata.var["mt_feature"].any() and not adata.var["hemoglobin_feature"].any(), "Excluded features remain.")
    require(adata.var["canonical_gene_id"].isna().all() and not adata.var["eligible_for_cross_dataset_integration"].any(), "Unverified gene mappings claim readiness.")
    require(adata.var["gene_id_analysis"].astype(str).tolist() == adata.var_names.tolist(), "Gene analysis IDs differ.")
    require(int((np.asarray(adata.X.sum(axis=0)).ravel() == 0).sum()) == EXPORT_EXPECTED_ZERO_GENES, "Zero genes were lost.")
    require(np.array_equal(adata.obs["total_counts_analysis"].to_numpy(), np.asarray(adata.X.sum(axis=1)).ravel()), "Per-cell exported counts differ.")
    for key in CONTRACT_KEYS:
        require(adata.uns[key]["schema_version"] == SCHEMA_VERSION, f"Wrong schema: {key}")
    require(not bool(adata.uns["integration_contract"]["model_ready"]), "Gene harmonization is pending.")
    contract = adata.uns["export_contract"]
    observed = {"ordered_cell_id_sha256": ids_sha256(adata.obs_names), "ordered_feature_id_sha256": ids_sha256(adata.var_names),
                "X_sha256": matrix_sha256(adata.X), "counts_sha256": matrix_sha256(adata.layers["counts"]),
                "obs_sha256": frame_sha256(adata.obs), "var_sha256": frame_sha256(adata.var)}
    for key, value in observed.items():
        require(contract[key] == value, f"Content fingerprint mismatch: {key}")


# Write once, read back, validate, then publish without overwriting
def export_object(adata, output=OUTPUT):
    output = Path(output).resolve()
    require(output.suffix.lower() == ".h5ad", "Output must have the .h5ad extension.")
    require(output != SOURCE_FILE.resolve(), "The original continuous-expression object must remain unchanged.")
    require(not output.exists(), f"Refusing to overwrite an existing output: {output}")
    require(output.parent.is_dir(), f"Output directory does not exist: {output.parent}")
    require(output.name == adata.uns["export_contract"]["output_filename"], "Output contract/filename mismatch.")
    validate_object(adata)
    expected_contracts = json_text(dict(adata.uns))
    temporary = WORK / ("export_counts_" + uuid.uuid4().hex + ".tmp.h5ad")
    require(not temporary.exists(), "Temporary filename collision.")
    try:
        adata.write_h5ad(temporary, compression="gzip")
        restored = ad.read_h5ad(temporary)
        validate_object(restored, adata.obs_names)
        require(json_text(dict(restored.uns)) == expected_contracts, "Contract values changed during H5AD serialization.")
        # Hard-link publication is atomic and fails if the target already exists.
        os.link(temporary, output)
        require(file_sha256(output) == file_sha256(temporary), "Published file checksum mismatch.")
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f"Wrote and read-validated: {output}")
    print(f"File SHA256: {file_sha256(output)}")



# Execute the authorized current stage

def main():
    parser = argparse.ArgumentParser(description="Validate E10.5 recovered counts, plot current QC and export a stagewise count object.")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--validate-only", action="store_true", help="Build and validate in memory without writing H5AD.")
    action.add_argument("--write", action="store_true", help="Write a new H5AD with readback validation; never overwrite.")
    action.add_argument("--verify-existing", action="store_true", help="Validate the existing H5AD and compare its counts with the verified raw source.")
    parser.add_argument("--plots", action="store_true", help="Render the two current QC diagnostic figures under Tmp/e10_counts.")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.verify_existing:
        result = ad.read_h5ad(args.output)
        counts, var, source, qc, run_summary = load_inputs()
        expected_counts, expected_var, excluded = build_var(var, counts)
        validate_object(result, source.obs_names)
        require(result.var_names.equals(expected_var.index), "Saved/source gene order differs.")
        require((result.layers["counts"] != expected_counts).nnz == 0, "Saved/raw source counts differ.")
        for key in ("stagewise_annotation_code", "stagewise_annotation", "integration_batch", "scanvi_label"):
            require(result.obs[key].astype(str).equals(source.obs[key].astype(str)), f"Original annotation/identity differs: {key}")
        for role, path in result.uns["source_provenance"]["source_files"].items():
            require(file_sha256(path) == result.uns["source_provenance"]["source_file_sha256"][role], f"Provenance input changed: {role}")
        print("Existing H5AD, source files, exact raw counts and original annotations verified.")
    else:
        if args.write:
            require(not args.output.exists(), f"Refusing to replace an existing output: {args.output}")
        result = build_object(args.output)
        if args.write:
            export_object(result, args.output)
    print(f"Verified object: {result.n_obs} cells x {result.n_vars} genes; sum={int(result.X.sum())}; nnz={result.X.nnz}")
    print("O-2=44, O-3=29; all 73 retained; verified read counts available; harmonization and training not run.")
    if args.plots:
        counts, var, source, qc, run_summary = load_inputs()
        print(f"Count QC figure: {plot_qc(qc)}")
        values, groups, cells, missing = marker_values()
        print(f"Marker QC figure: {plot_marker_qc(values, groups, cells)}")


if __name__ == "__main__":
    main()
