#!/usr/bin/env python3
"""CodeSpringLab single-cell RNA-seq workflow (Scanpy engine).

Usage: python scrna_pipeline_scanpy.py <samples.tsv> <out_dir> <params.tsv>

The manifest and parameter files are shared with the Seurat engine.  This
script intentionally preserves raw counts in ``layers['counts']`` and uses a
portable Harmony representation for automatic multi-sample integration.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

# Batch jobs and local tests have no interactive display.  Explicitly use the
# non-interactive backend before Scanpy imports matplotlib.
os.environ.setdefault("MPLBACKEND", "Agg")


def require(module: str, package: str | None = None):
    try:
        return __import__(module)
    except ImportError as exc:
        name = package or module
        raise SystemExit(
            f"The Scanpy engine requires {name}. Configure a Python environment with "
            "scanpy, anndata, numpy, pandas, matplotlib, scipy, and (for integration) scvi-tools or harmonypy."
        ) from exc


pd = require("pandas")
np = require("numpy")
sc = require("scanpy")
ad = require("anndata")
mpl = require("matplotlib")
try:
    # Bundled with CodeSpringLab so Scanpy figures use the same established
    # jpplot color treatment regardless of the user's working directory.
    import jpplot
except ImportError as exc:
    raise SystemExit("The bundled jpplot.py color module is missing from CodeSpringLab.") from exc

JP_COLOR_MAP = jpplot.cmapjp()
# Newer Scanpy/Seaborn releases require a named colormap in the global plotting
# settings; passing a ListedColormap object makes violin plots fail.  Register
# the bundled palette once and use its stable name in every Scanpy plot.
JP_COLOR_MAP_NAME = "codespring_jp"
try:
    mpl.colormaps.register(JP_COLOR_MAP, name=JP_COLOR_MAP_NAME)
except ValueError:
    # The name may already be registered when a worker reuses a Python process.
    pass
sc.settings.set_figure_params(color_map=JP_COLOR_MAP_NAME, dpi=160)


def apply_jpplot_colors(adata, *keys):
    """Apply jpplot's palette consistently to categorical Scanpy plots."""
    from matplotlib.colors import to_hex
    for key in keys:
        if key not in adata.obs.columns:
            continue
        values = adata.obs[key]
        is_categorical = isinstance(values.dtype, pd.CategoricalDtype)
        if not (is_categorical or values.dtype == object or pd.api.types.is_string_dtype(values)):
            continue
        if not is_categorical:
            adata.obs[key] = pd.Categorical(values.astype(str))
            values = adata.obs[key]
        categories = [str(value) for value in values.cat.categories if str(value)]
        if not categories:
            continue
        # Evenly spaced colors make the palette deterministic across all
        # generated Scanpy UMAPs, violins, and sample-level scatters.
        positions = np.linspace(0.10, 0.90, len(categories))
        adata.uns[f"{key}_colors"] = [to_hex(JP_COLOR_MAP(position)) for position in positions]


def read_table(path: Path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def params_from(path: Path):
    tab = read_table(path)
    if not {"key", "value"}.issubset(tab.columns):
        raise SystemExit("Parameter file must contain key and value columns.")
    d = dict(zip(tab["key"], tab["value"]))
    def get(key, default=""):
        return str(d.get(key, default)).strip()
    def get_bool(key, default=False):
        return get(key, "true" if default else "false").lower() in {"1", "true", "t", "yes", "y"}
    return {
        "normalization": get("normalization", "lognormalize").lower(),
        "integration": get("integration", "auto").lower(),
        "batch_column": get("batch_column", "batch"),
        "cluster_resolution": float(get("cluster_resolution", "0.6") or 0.6),
        "n_pcs": int(float(get("n_pcs", "30") or 30)),
        "min_features": int(float(get("min_features", "200") or 200)),
        "min_counts": int(float(get("min_counts", "0") or 0)),
        "max_features": int(float(get("max_features", "0") or 0)),
        "max_percent_mt": float(get("max_percent_mt", "20") or 20),
        "min_cells_per_gene": int(float(get("min_cells_per_gene", "3") or 3)),
        "doublet_method": get("doublet_method", "auto").lower(),
        "doublet_rate": float(get("doublet_rate", "0.05") or 0.05),
        "remove_doublets": get_bool("remove_doublets", True),
        "marker_file": get("marker_file", ""),
        "celltype_file": get("celltype_file", ""),
        "scvi_max_epochs": int(float(get("scvi_max_epochs", "400") or 400)),
        "seed": int(float(get("seed", "1234") or 1234)),
    }


def input_kind(path: Path):
    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        return "scanpy_h5ad"
    if path.is_dir():
        return "filtered_10x_matrix"
    raise SystemExit(f"Unsupported input type: {path}. Scanpy accepts .h5ad and filtered 10x matrix directories.")


def read_10x_matrix(path: Path):
    try:
        # Never place a Scanpy cache next to a user-owned source matrix.
        return sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    except FileNotFoundError:
        matrix = next((x for x in (path / "matrix.mtx", path / "matrix.mtx.gz") if x.exists()), None)
        features = next((x for x in (path / "features.tsv", path / "features.tsv.gz", path / "genes.tsv", path / "genes.tsv.gz") if x.exists()), None)
        barcodes = next((x for x in (path / "barcodes.tsv", path / "barcodes.tsv.gz") if x.exists()), None)
        if not matrix or not features or not barcodes:
            raise SystemExit(f"A filtered 10x matrix directory needs matrix.mtx(.gz), features.tsv/genes.tsv(.gz), and barcodes.tsv(.gz): {path}")
        x = sc.read_mtx(matrix).T
        feature_frame = pd.read_csv(features, sep="\t", header=None, dtype=str)
        barcode_frame = pd.read_csv(barcodes, sep="\t", header=None, dtype=str)
        x.var_names = feature_frame.iloc[:, 1 if feature_frame.shape[1] > 1 else 0].astype(str).values
        x.obs_names = barcode_frame.iloc[:, 0].astype(str).values
        if feature_frame.shape[1] > 1:
            x.var["gene_ids"] = feature_frame.iloc[:, 0].astype(str).values
        return x


def counts_like(matrix) -> bool:
    """Conservatively distinguish raw UMI counts from normalized expression."""
    values = getattr(matrix, "data", matrix)
    values = np.asarray(values).ravel()
    if values.size == 0:
        return True
    if values.size > 100_000:
        values = values[np.linspace(0, values.size - 1, 100_000, dtype=int)]
    return bool(np.nanmin(values) >= 0 and np.allclose(values, np.rint(values), rtol=0, atol=1e-8))


def categorical_plot_column(adata):
    preferred = [
        "cell_type", "celltype", "CellType", "annotation", "annotated_cell_type", "predicted.celltype",
        "major_celltype", "consolidated", "subtype", "condition", "Condition", "treatment", "Treatment",
        "group", "Group", "sample", "sample_id", "batch",
    ]
    for column in list(dict.fromkeys(preferred + list(adata.obs.columns))):
        if column not in adata.obs:
            continue
        values = adata.obs[column]
        if pd.api.types.is_numeric_dtype(values):
            continue
        n_values = values.astype(str).replace("", np.nan).dropna().nunique()
        if 1 < n_values <= 40:
            return column
    return None


def input_status_one(adata, sample_id: str, kind: str, raw_count_source: str):
    layer_names = [str(name) for name in adata.layers.keys()]
    embedding_names = [str(name) for name in adata.obsm.keys()]
    annotation_columns = [name for name in [
        "cell_type", "celltype", "CellType", "annotation", "annotated_cell_type", "predicted.celltype",
        "major_celltype", "consolidated", "subtype",
    ] if name in adata.obs.columns]
    return {
        "sample_id": sample_id,
        "input_kind": kind,
        "cells_input": adata.n_obs,
        "features_input": adata.n_vars,
        # Keep this inspection inventory schema aligned with the Seurat
        # runner.  The app intentionally renders one engine-independent
        # yes/no source report, so different names here silently made H5AD
        # results look as if their raw counts were absent.
        "assays_detected": "AnnData",
        "rna_count_layers_detected": raw_count_source,
        "rna_normalized_layers_detected": "; ".join(
            name for name in layer_names if name.lower() not in {"counts", "raw"}
        ),
        "sct_detected": False,
        "integrated_detected": any(
            name.lower() in {"x_harmony", "x_scvi", "x_scanvi"} for name in embedding_names
        ),
        "reductions_detected": "; ".join(embedding_names),
        "raw_count_source": raw_count_source,
        "layers_detected": "; ".join(layer_names),
        "raw_slot_detected": adata.raw is not None,
        "embeddings_detected": "; ".join(embedding_names),
        "pca_detected": any(name.lower() in {"x_pca", "pca"} for name in embedding_names),
        "umap_detected": any(name.lower() in {"x_umap", "umap"} for name in embedding_names),
        "clusters_detected": any(name in adata.obs.columns for name in ["leiden", "louvain", "cluster", "seurat_clusters"]),
        "annotation_columns_detected": "; ".join(annotation_columns),
        "workflow_action": "Raw counts retained; CodeSpring reruns selected QC and downstream processing reproducibly.",
    }


def save_input_umap(adata, sample_id: str, figures: Path):
    if not any(name.lower() in {"x_umap", "umap"} for name in adata.obsm.keys()):
        return
    import matplotlib.pyplot as plt
    color = categorical_plot_column(adata)
    apply_jpplot_colors(adata, color)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color=color, ax=ax, show=False, title="UMAP supplied with input object", frameon=False, color_map=JP_COLOR_MAP_NAME)
    fig.tight_layout()
    safe_sample = "".join(char if char.isalnum() or char in "._-" else "_" for char in sample_id)
    fig.savefig(figures / f"00_input_umap_{safe_sample}.png", dpi=160)
    plt.close(fig)


def read_one(row, figures: Path):
    path = Path(row.input_path).expanduser()
    if not path.exists():
        raise SystemExit(f"Input path does not exist: {path}")
    kind = input_kind(path)
    if kind == "scanpy_h5ad":
        x = sc.read_h5ad(path)
        raw_layer = next((name for name in ("counts", "raw") if name in x.layers and counts_like(x.layers[name])), None)
        if raw_layer is not None:
            raw_count_source = f"layers['{raw_layer}']"
            status = input_status_one(x, str(row.sample_id), kind, raw_count_source)
            save_input_umap(x, str(row.sample_id), figures)
            x.X = x.layers[raw_layer].copy()
        elif x.raw is not None:
            raw_count_source = "raw"
            status = input_status_one(x, str(row.sample_id), kind, raw_count_source)
            save_input_umap(x, str(row.sample_id), figures)
            x = x.raw.to_adata()
        elif counts_like(x.X):
            raw_count_source = "X"
            status = input_status_one(x, str(row.sample_id), kind, raw_count_source)
            save_input_umap(x, str(row.sample_id), figures)
        else:
            raw_count_source = "unavailable"
            status = input_status_one(x, str(row.sample_id), kind, raw_count_source)
            save_input_umap(x, str(row.sample_id), figures)
        if not counts_like(x.X):
            raise SystemExit(
                f"{path} does not expose raw integer-like counts in X, layers['counts'/'raw'], or raw. "
                "A full QC/normalization workflow must start from raw counts; export a raw-count AnnData object or use a filtered 10x matrix."
            )
    else:
        x = read_10x_matrix(path)
        raw_count_source = "X"
        status = input_status_one(x, str(row.sample_id), kind, raw_count_source)
    x.var_names_make_unique()
    x.obs_names_make_unique()
    x.obs["source_barcode"] = x.obs_names.astype(str)
    x.obs_names = [f"{row.sample_id}_{barcode}" for barcode in x.obs_names]
    for key, value in row._asdict().items():
        x.obs[key] = str(value)
    x.obs["input_kind"] = kind
    return x, status


def save_umap(adata, color, path, title=None):
    import matplotlib.pyplot as plt
    apply_jpplot_colors(adata, color)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sc.pl.umap(adata, color=color, ax=ax, show=False, title=title or str(color), frameon=False, color_map=JP_COLOR_MAP_NAME)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_qc_plots(adata, figures: Path, prefix: str = "01_qc"):
    """Write clean, readable pre/post-filter QC figures without an interactive display."""
    import matplotlib.pyplot as plt
    sample_series = adata.obs["sample_id"].astype(str)
    samples = list(dict.fromkeys(sample_series.tolist()))
    display_labels = [sample if len(sample) <= 28 else f"{sample[:25]}…" for sample in samples]
    single_sample = len(samples) == 1
    colors = [JP_COLOR_MAP(position) for position in np.linspace(0.18, 0.82, max(len(samples), 2))]
    metrics = [
        ("n_genes_by_counts", "Detected genes per cell"),
        ("total_counts", "UMIs per cell"),
        ("pct_counts_mt", "Mitochondrial reads (%)"),
    ]
    fig, axes = plt.subplots(1, len(metrics), figsize=(14, 4.8), layout="constrained")
    for axis, (column, label) in zip(axes, metrics):
        groups = [pd.to_numeric(adata.obs.loc[sample_series == sample, column], errors="coerce").dropna().to_numpy() for sample in samples]
        groups = [group if len(group) else np.array([0.0]) for group in groups]
        violin = axis.violinplot(groups, showmedians=True, showextrema=False)
        for body, color in zip(violin["bodies"], colors):
            body.set_facecolor(color)
            body.set_edgecolor("none")
            body.set_alpha(0.82)
        violin["cmedians"].set_color("#1f2937")
        violin["cmedians"].set_linewidth(1.25)
        if single_sample:
            axis.set_xticks([])
        else:
            axis.set_xticks(range(1, len(samples) + 1), display_labels, rotation=28, ha="right")
        axis.set_ylabel(label)
        axis.set_title(label, loc="left", fontweight="bold")
        axis.grid(axis="y", color="#d1d5db", linewidth=0.7, alpha=0.8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(axis="x", length=0)
    if single_sample:
        fig.suptitle(f"QC overview — {display_labels[0]}", fontsize=14, fontweight="bold")
    fig.savefig(figures / f"{prefix}_violin.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.6, 5.6), layout="constrained")
    for sample, color, label in zip(samples, colors, display_labels):
        subset = adata.obs.loc[sample_series == sample]
        ax.scatter(
            pd.to_numeric(subset["total_counts"], errors="coerce"),
            pd.to_numeric(subset["pct_counts_mt"], errors="coerce"),
            s=5, alpha=0.32, color=color, edgecolors="none", rasterized=True, label=label,
        )
    ax.set_xlabel("UMIs per cell")
    ax.set_ylabel("Mitochondrial reads (%)")
    ax.set_title("Library size versus mitochondrial content", loc="left", fontweight="bold")
    ax.grid(color="#d1d5db", linewidth=0.7, alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    if len(samples) > 1:
        ax.legend(title="Sample", frameon=False, markerscale=2.2, loc="best")
    fig.savefig(figures / f"{prefix}_counts_vs_mt.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def write_h5ad_checkpoint(adata, path: Path):
    """Write checkpoints across AnnData/Pandas versions without altering values.

    Some source H5AD files use Pandas' nullable ``StringArray`` for an axis
    index or metadata column. Older widely deployed AnnData releases reject
    that representation when writing a new H5AD. Converting it to ordinary
    Python strings is lossless and keeps an imported processed object usable.
    """
    def normalize_frame(frame):
        if isinstance(frame.index.dtype, pd.StringDtype):
            frame.index = pd.Index(frame.index.astype(str).to_numpy(dtype=object), dtype=object)
        for column in frame.columns:
            if isinstance(frame[column].dtype, pd.StringDtype):
                frame[column] = frame[column].astype(object)

    normalize_frame(adata.obs)
    normalize_frame(adata.var)
    if adata.raw is not None:
        normalize_frame(adata.raw.var)
    adata.write_h5ad(path, compression="gzip")


def save_doublet_plot(adata, figures: Path):
    """Show the score distribution used for a transparent doublet call."""
    import matplotlib.pyplot as plt
    if "doublet_score" not in adata.obs:
        return
    score = pd.to_numeric(adata.obs["doublet_score"], errors="coerce").dropna()
    if score.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(score, bins=min(50, max(10, int(np.sqrt(len(score))))), color=JP_COLOR_MAP(0.72), edgecolor="white")
    ax.set(xlabel="Doublet score", ylabel="Cells", title="Doublet-score distribution before removal")
    fig.tight_layout()
    fig.savefig(figures / "02b_doublet_scores.png", dpi=160)
    plt.close(fig)


def run_scrublet(adata, p, tables: Path, figures: Path):
    """Run Scrublet on raw counts, separately by sample when available."""
    method = p["doublet_method"]
    if method not in {"auto", "none", "scrublet"}:
        raise SystemExit("Scanpy doublet_method must be auto, none, or scrublet.")
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    report = []
    if method == "none":
        report.append({"sample_id": "all", "method": "none", "expected_doublet_rate": p["doublet_rate"], "cells_before": adata.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": "Doublet removal disabled"})
        pd.DataFrame(report).to_csv(tables / "doublet_summary_by_sample.tsv", sep="\t", index=False)
        calls = adata.obs[["sample_id", "doublet_score", "predicted_doublet"]].copy()
        calls.insert(0, "cell", adata.obs_names)
        calls["removed_as_doublet"] = False
        calls.to_csv(tables / "doublet_calls.tsv", sep="\t", index=False)
        adata.uns["codespring_doublets_removed"] = 0
        return adata
    sample_col = "sample_id" if "sample_id" in adata.obs else None
    groups = [(str(name), list(idx)) for name, idx in adata.obs.groupby(sample_col, observed=True).groups.items()] if sample_col else [("all", list(adata.obs_names))]
    for sample_id, cells in groups:
        sub = adata[cells].copy()
        expected_rate = p["doublet_rate"]
        if "expected_doublet_rate" in sub.obs.columns:
            rates = pd.to_numeric(sub.obs["expected_doublet_rate"], errors="coerce").dropna().unique()
            if len(rates) > 1:
                raise SystemExit(f"expected_doublet_rate must have one value per sample ({sample_id}).")
            if len(rates) == 1:
                expected_rate = float(rates[0])
                if not 0 < expected_rate < 1:
                    raise SystemExit(f"expected_doublet_rate must be greater than 0 and less than 1 for {sample_id}.")
        # Scrublet is not informative on very small partitions. Do not invent
        # calls; retain these cells and record why in the summary.
        if sub.n_obs < 100 or sub.n_vars < 20:
            report.append({"sample_id": sample_id, "method": "scrublet", "expected_doublet_rate": expected_rate, "cells_before": sub.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": "Skipped: fewer than 100 cells or 20 genes"})
            continue
        try:
            # The internal Scrublet feature filter can leave fewer than 30
            # dimensions, so choose a conservative PC count for small runs.
            scrublet_pcs = min(20, max(2, min(sub.n_obs, sub.n_vars) - 2))
            sc.pp.scrublet(sub, expected_doublet_rate=expected_rate, n_prin_comps=scrublet_pcs, random_state=p["seed"], verbose=False)
            adata.obs.loc[sub.obs_names, "doublet_score"] = pd.to_numeric(sub.obs["doublet_score"], errors="coerce").values
            predicted = sub.obs["predicted_doublet"].astype(bool)
            adata.obs.loc[sub.obs_names, "predicted_doublet"] = predicted.values
            predicted_n = int(predicted.sum())
            report.append({"sample_id": sample_id, "method": "scrublet", "expected_doublet_rate": expected_rate, "cells_before": sub.n_obs, "predicted_doublets": predicted_n, "removed_doublets": predicted_n if p["remove_doublets"] else 0, "note": ""})
        except Exception as exc:
            if method == "scrublet":
                raise SystemExit(f"Scrublet failed for {sample_id}: {exc}") from exc
            report.append({"sample_id": sample_id, "method": "scrublet", "expected_doublet_rate": expected_rate, "cells_before": sub.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": f"Automatic Scrublet skipped after error: {exc}"})
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
    save_doublet_plot(adata, figures)
    pd.DataFrame(report).to_csv(tables / "doublet_summary_by_sample.tsv", sep="\t", index=False)
    calls = adata.obs[["sample_id", "doublet_score", "predicted_doublet"]].copy()
    calls.insert(0, "cell", adata.obs_names)
    calls["removed_as_doublet"] = calls["predicted_doublet"] & p["remove_doublets"]
    calls.to_csv(tables / "doublet_calls.tsv", sep="\t", index=False)
    adata.uns["codespring_doublets_removed"] = int(calls["removed_as_doublet"].sum())
    if p["remove_doublets"]:
        adata = adata[~adata.obs["predicted_doublet"]].copy()
    if adata.n_obs == 0:
        raise SystemExit("Doublet removal removed every cell; use a lower expected doublet rate or review the input data.")
    return adata


def apply_celltype_mapping(adata, path: Path):
    tab = read_table(path)
    cell_col = next((x for x in ["cell", "cell_id", "barcode", "Cell", "CellID"] if x in tab.columns), None)
    type_col = next((x for x in ["cell_type", "celltype", "CellType", "annotation"] if x in tab.columns), None)
    if not cell_col or not type_col:
        raise SystemExit("Cell-type mapping needs a cell/barcode column and a cell_type column.")
    labels = dict(zip(tab[cell_col].astype(str), tab[type_col].astype(str)))
    values = [labels.get(str(cell), labels.get(str(raw), "Unassigned")) for cell, raw in zip(adata.obs_names, adata.obs["source_barcode"])]
    adata.obs["cell_type"] = pd.Categorical(values)
    adata.obs["annotation_source"] = "provided cell metadata"


def apply_marker_annotation(adata, path: Path, tables_dir: Path):
    tab = read_table(path)
    if not {"cell_type", "gene"}.issubset(tab.columns):
        raise SystemExit("Marker list needs cell_type and gene columns.")
    lists = {}
    genes = set(adata.var_names)
    for cell_type, group in tab.groupby("cell_type", sort=False):
        entries = []
        for value in group["gene"].astype(str):
            entries.extend([g for g in value.replace(";", ",").replace(" ", ",").split(",") if g])
        keep = sorted(set(entries).intersection(genes))
        if keep:
            lists[str(cell_type)] = keep
    if not lists:
        raise SystemExit("None of the supplied marker genes were present in the expression object.")
    for name, genes_for_type in lists.items():
        sc.tl.score_genes(adata, genes_for_type, score_name=f"marker_score__{name}", use_raw=False)
    clusters = adata.obs["cluster"].astype(str)
    score_columns = [f"marker_score__{x}" for x in lists]
    means = adata.obs.assign(cluster=clusters).groupby("cluster", observed=True)[score_columns].mean()
    winners = means.idxmax(axis=1).str.replace("marker_score__", "", regex=False)
    adata.obs["cell_type"] = pd.Categorical(clusters.map(winners).fillna("Unassigned"))
    adata.obs["annotation_source"] = "marker-list cluster scoring"
    means.to_csv(tables_dir / "marker_annotation_cluster_scores.tsv", sep="\t")


def apply_existing_annotation(adata):
    candidates = [
        "cell_type", "celltype", "CellType", "annotation", "annotated_cell_type", "predicted.celltype",
        "major_celltype", "consolidated", "subtype",
    ]
    for column in candidates:
        if column not in adata.obs:
            continue
        values = adata.obs[column].astype(str)
        if values.replace("", np.nan).notna().any():
            adata.obs["cell_type"] = pd.Categorical(values.where(values != "", "Unassigned"))
            adata.obs["annotation_source"] = f"existing input metadata: {column}"
            return True
    return False


def main():
    if len(sys.argv) not in {4, 5}:
        raise SystemExit("Usage: scrna_pipeline_scanpy.py <samples.tsv> <out_dir> <params.tsv> [inspect|qc|preprocess|cluster|annotate|all]")
    samples_path, out_dir, params_path = map(Path, sys.argv[1:4])
    stage = sys.argv[4].lower() if len(sys.argv) == 5 else "all"
    stages = {"inspect", "qc", "preprocess", "cluster", "annotate", "all"}
    if stage not in stages:
        raise SystemExit("Unknown scRNA stage: " + stage)
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = out_dir / "tables"; figures = out_dir / "figures"; objects = out_dir / "objects"; checkpoints = out_dir / "checkpoints"
    for d in (tables, figures, objects, checkpoints): d.mkdir(parents=True, exist_ok=True)
    p = params_from(params_path)
    if p["normalization"] not in {"lognormalize", "log1p"}:
        raise SystemExit("The Scanpy engine supports LogNormalize/log1p normalization. Use the Seurat engine for SCTransform.")
    if p["integration"] not in {"auto", "none", "scvi", "harmony"}:
        raise SystemExit("Scanpy integration must be auto, none, scvi, or harmony.")
    np.random.seed(p["seed"])
    input_checkpoint = checkpoints / "01_input_scanpy.h5ad"
    qc_checkpoint = checkpoints / "02_qc_scanpy.h5ad"
    preprocess_checkpoint = checkpoints / "03_preprocessed_scanpy.h5ad"
    cluster_checkpoint = checkpoints / "04_clustered_scanpy.h5ad"

    def mark_complete(name):
        (out_dir / f"_STAGE_{name.upper()}_COMPLETE").write_text("complete\n")

    def require_checkpoint(path: Path, prior: str):
        if not path.exists():
            raise SystemExit(f"The {prior} stage has not completed. Run {prior} before this stage.")
        return sc.read_h5ad(path)

    # Stage 1: inspect each source object and preserve a raw-count checkpoint.
    if stage in {"inspect", "all"}:
        samples = read_table(samples_path)
        if not {"sample_id", "input_path"}.issubset(samples.columns):
            raise SystemExit("Sample manifest requires sample_id and input_path columns.")
        samples = samples[(samples.sample_id.str.strip() != "") & (samples.input_path.str.strip() != "")].copy()
        if samples.empty:
            raise SystemExit("No valid sample rows were supplied.")
        samples["sample_id"] = samples["sample_id"].astype(str).str.strip()
        if samples["sample_id"].duplicated().any():
            duplicates = ", ".join(samples.loc[samples["sample_id"].duplicated(), "sample_id"].unique())
            raise SystemExit(f"sample_id values must be unique in the input manifest: {duplicates}")
        samples.to_csv(tables / "sample_manifest_used.tsv", sep="\t", index=False)
        inputs = [read_one(row, figures) for row in samples.itertuples(index=False)]
        adata = ad.concat([item[0] for item in inputs], join="outer", merge="same", index_unique=None, fill_value=0)
        adata.var_names_make_unique()
        pd.DataFrame([item[1] for item in inputs]).to_csv(tables / "input_processing_detected.tsv", sep="\t", index=False)
        # Inspect produces an unfiltered QC overview first. Users can then
        # choose thresholds from the actual per-sample distributions instead
        # of applying generic defaults blindly.
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        pre_qc = adata.obs[["sample_id", "n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
        pre_qc.insert(0, "cell", adata.obs_names)
        pre_qc.to_csv(tables / "qc_pre_filter_cell_metrics.tsv", sep="\t", index=False)
        save_qc_plots(adata, figures, prefix="00_qc_pre_filter")
        write_h5ad_checkpoint(adata, input_checkpoint)
        mark_complete("inspect")
        if stage == "inspect":
            return
    else:
        adata = require_checkpoint(input_checkpoint, "input inspection")
        samples = read_table(samples_path)

    # Stage 2: filter cells/genes and handle doublets while the data is raw.
    if stage in {"qc", "all"}:
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        cells_before_qc = adata.obs.groupby("sample_id", observed=True).size().rename("cells_input")
        keep = (adata.obs["n_genes_by_counts"] >= p["min_features"]) & (adata.obs["total_counts"] >= p["min_counts"]) & (adata.obs["pct_counts_mt"] <= p["max_percent_mt"])
        if p["max_features"] > 0:
            keep &= adata.obs["n_genes_by_counts"] <= p["max_features"]
        adata = adata[keep].copy()
        if adata.n_obs == 0:
            raise SystemExit("QC thresholds removed every cell.")
        adata = run_scrublet(adata, p, tables, figures)
        gene_keep = np.zeros(adata.n_vars, dtype=bool)
        feature_reports = []
        for sample_id, cells in adata.obs.groupby("sample_id", observed=True).groups.items():
            sub = adata[list(cells)]
            detected_cells = np.asarray((sub.X > 0).sum(axis=0)).ravel()
            keep_here = detected_cells >= p["min_cells_per_gene"]
            gene_keep |= keep_here
            feature_reports.append({"sample_id": str(sample_id), "genes_before_filtering": int(adata.n_vars), "min_cells_per_gene": p["min_cells_per_gene"], "genes_retained": int(keep_here.sum())})
        if not gene_keep.any():
            raise SystemExit("Minimum cells per retained gene removed every gene.")
        adata = adata[:, gene_keep].copy()
        pd.DataFrame(feature_reports).to_csv(tables / "feature_filtering_by_sample.tsv", sep="\t", index=False)
        adata.layers["counts"] = adata.X.copy()
        qc = adata.obs[["sample_id", "n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_score", "predicted_doublet"]].copy()
        qc.insert(0, "cell", adata.obs_names)
        qc.to_csv(tables / "qc_cell_metrics.tsv", sep="\t", index=False)
        qc_summary = qc.groupby("sample_id", observed=True).agg(cells_after_qc_and_doublets=("total_counts", "size"), median_umis=("total_counts", "median"), median_genes=("n_genes_by_counts", "median"), median_percent_mt=("pct_counts_mt", "median")).reset_index()
        qc_summary.insert(1, "cells_input", qc_summary["sample_id"].map(cells_before_qc).astype(int))
        qc_summary.to_csv(tables / "qc_summary_by_sample.tsv", sep="\t", index=False)
        save_qc_plots(adata, figures)
        write_h5ad_checkpoint(adata, qc_checkpoint)
        mark_complete("qc")
        if stage == "qc":
            return
    elif stage not in {"inspect"}:
        adata = require_checkpoint(qc_checkpoint, "QC and doublet handling")

    # Stage 3: normalization, highly variable features, scaling, and PCA.
    if stage in {"preprocess", "all"}:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        hvg_batch_key = p["batch_column"] if p["batch_column"] in adata.obs.columns and adata.obs[p["batch_column"]].astype(str).nunique() > 1 else None
        hvg_kwargs = {"n_top_genes": min(3000, adata.n_vars), "layer": "counts"}
        fallback_hvg_kwargs = {"n_top_genes": min(3000, adata.n_vars)}
        if hvg_batch_key is not None:
            hvg_kwargs["batch_key"] = hvg_batch_key
            fallback_hvg_kwargs["batch_key"] = hvg_batch_key
        try:
            sc.pp.highly_variable_genes(adata, flavor="seurat_v3", **hvg_kwargs)
        except Exception:
            sc.pp.highly_variable_genes(adata, flavor="cell_ranger", **fallback_hvg_kwargs)
        hvg_table = adata.var.copy(); hvg_table.insert(0, "gene", hvg_table.index.astype(str))
        hvg_columns = [c for c in ["gene", "highly_variable", "means", "variances", "variances_norm"] if c in hvg_table.columns]
        hvg_table.loc[:, hvg_columns].to_csv(tables / "highly_variable_genes.tsv", sep="\t", index=False)
        sc.pp.scale(adata, zero_center=False, max_value=10)
        hvg_count = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else 0
        use_hvg = hvg_count >= 2
        pca_features = hvg_count if use_hvg else adata.n_vars
        max_pcs = max(1, min(adata.n_obs - 1, pca_features - 1))
        n_pcs = min(p["n_pcs"], max_pcs)
        sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", use_highly_variable=use_hvg, zero_center=False)
        pca_ratio = np.asarray(adata.uns["pca"]["variance_ratio"]).ravel()
        pd.DataFrame({"PC": np.arange(1, len(pca_ratio) + 1), "variance_explained": pca_ratio, "percent_variance_explained": 100 * pca_ratio}).to_csv(tables / "pca_variance_explained.tsv", sep="\t", index=False)
        write_h5ad_checkpoint(adata, preprocess_checkpoint)
        mark_complete("preprocess")
        if stage == "preprocess":
            return
    elif stage not in {"inspect", "qc"}:
        adata = require_checkpoint(preprocess_checkpoint, "normalization and PCA")

    # Stage 4: optional technical-batch integration followed by neighbours, UMAP, and clusters.
    if stage in {"cluster", "all"}:
        n_pcs = min(p["n_pcs"], adata.obsm["X_pca"].shape[1])
        integration = p["integration"]
        batch_key = p["batch_column"] if p["batch_column"] in adata.obs.columns else "sample_id"
        batch_count = int(adata.obs[batch_key].astype(str).nunique())
        if integration == "auto":
            # Harmony is the automatic choice because it is available in the
            # standard Scanpy runtime. scVI remains an explicit option when a
            # dedicated scvi-tools environment has been installed.
            integration = "harmony" if p["batch_column"] in adata.obs.columns and batch_count > 1 else "none"
        if integration in {"scvi", "harmony"} and batch_count < 2:
            raise SystemExit(f"{integration} integration requires at least two values in the selected batch column ({batch_key}). Choose none or supply the appropriate technical batch column.")
        representation = "X_pca"
        if integration == "scvi" and batch_count > 1:
            try:
                import scvi
            except ImportError as exc:
                raise SystemExit("scVI integration was requested but scvi-tools is unavailable. Install scvi-tools or choose Harmony/none.") from exc
            scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
            model = scvi.model.SCVI(adata, n_latent=min(30, p["n_pcs"]))
            model.train(max_epochs=p["scvi_max_epochs"], early_stopping=True)
            adata.obsm["X_scVI"] = model.get_latent_representation(); representation = "X_scVI"
        elif integration == "harmony" and batch_count > 1:
            try:
                import scanpy.external as sce
                sce.pp.harmony_integrate(adata, key=batch_key, basis="X_pca", adjusted_basis="X_harmony")
            except Exception as exc:
                raise SystemExit("Harmony integration was requested but harmonypy was unavailable or failed: " + str(exc)) from exc
            representation = "X_harmony"
        sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs - 1)), n_pcs=min(n_pcs, adata.obsm[representation].shape[1]), use_rep=representation)
        sc.tl.umap(adata, random_state=p["seed"])
        sc.tl.leiden(adata, resolution=p["cluster_resolution"], key_added="cluster", random_state=p["seed"])
        adata.uns["codespring_integration"] = integration
        write_h5ad_checkpoint(adata, cluster_checkpoint)
        mark_complete("cluster")
        if stage == "cluster":
            return
    elif stage == "annotate":
        adata = require_checkpoint(cluster_checkpoint, "integration and clustering")

    # Stage 5: annotation, markers, and the final portable results set.
    if stage in {"annotate", "all"}:
        integration = str(adata.uns.get("codespring_integration", p["integration"]))
        if p["celltype_file"] and p["celltype_file"].lower() != "none":
            apply_celltype_mapping(adata, Path(p["celltype_file"]).expanduser())
        elif p["marker_file"] and p["marker_file"].lower() != "none":
            apply_marker_annotation(adata, Path(p["marker_file"]).expanduser(), tables)
        elif not apply_existing_annotation(adata):
            adata.obs["cell_type"] = adata.obs["cluster"].astype("category")
            adata.obs["annotation_source"] = "cluster ID (no annotation supplied)"
        for col, name in [("sample_id", "03_umap_sample.png"), ("cluster", "04_umap_clusters.png"), ("cell_type", "05_umap_cell_types.png")]:
            save_umap(adata, col, figures / name)
        if "condition" in adata.obs.columns and adata.obs.condition.nunique() > 1:
            save_umap(adata, "condition", figures / "06_umap_condition.png")
        try:
            sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", use_raw=True)
            markers = sc.get.rank_genes_groups_df(adata, group=None)
            markers.to_csv(tables / "cluster_markers.tsv", sep="\t", index=False)
            markers.groupby("group", observed=True).head(10).to_csv(tables / "top10_markers_per_cluster.tsv", sep="\t", index=False)
        except Exception as exc:
            pd.DataFrame({"warning": [str(exc)]}).to_csv(tables / "cluster_markers.tsv", sep="\t", index=False)
        cell_metadata = adata.obs.copy()
        if "cell" in cell_metadata.columns: cell_metadata = cell_metadata.rename(columns={"cell": "input_cell"})
        cell_metadata.insert(0, "cell", cell_metadata.index.astype(str))
        cell_metadata.to_csv(tables / "cell_metadata.tsv", sep="\t", index=False)
        umap_table = pd.DataFrame(adata.obsm["X_umap"][:, :2], index=adata.obs_names, columns=["UMAP_1", "UMAP_2"])
        umap_metadata = adata.obs.copy()
        if "cell" in umap_metadata.columns: umap_metadata = umap_metadata.rename(columns={"cell": "input_cell"})
        umap_table = pd.concat([umap_table, umap_metadata], axis=1); umap_table.insert(0, "cell", umap_table.index.astype(str))
        umap_table.to_csv(tables / "umap_coordinates.tsv", sep="\t", index=False)
        adata.obs.groupby(["cluster", "cell_type"], observed=True).size().reset_index(name="cells").to_csv(tables / "cluster_cell_type_sizes.tsv", sep="\t", index=False)
        cell_type_by_sample = adata.obs.groupby(["sample_id", "cell_type"], observed=True).size().reset_index(name="cells")
        cell_type_by_sample["proportion_within_sample"] = cell_type_by_sample["cells"] / cell_type_by_sample.groupby("sample_id", observed=True)["cells"].transform("sum")
        cell_type_by_sample.to_csv(tables / "cell_type_by_sample.tsv", sep="\t", index=False)
        write_h5ad_checkpoint(adata, objects / "processed_scanpy.h5ad")
        (out_dir / "run_summary.txt").write_text("\n".join([
            "engine: scanpy", "normalization: lognormalize", f"integration: {integration}", f"doublet_method: {p['doublet_method']}",
            f"doublets_removed: {int(adata.uns.get('codespring_doublets_removed', 0))}", "input_processing_inventory: tables/input_processing_detected.tsv",
            f"input_samples: {samples.shape[0]}", f"cells_after_qc: {adata.n_obs}", f"clusters: {adata.obs['cluster'].nunique()}", f"annotation_source: {adata.obs['annotation_source'].iloc[0]}",
        ]) + "\n")
        mark_complete("annotate")
        (out_dir / "_COMPLETE").write_text("complete\n")


if __name__ == "__main__":
    main()
