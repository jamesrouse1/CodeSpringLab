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
import re
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


QUALITATIVE_COLORS = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#6A3D9A",
    "#56B4E9", "#B15928", "#E31A1C", "#1B9E77", "#7570B3", "#E7298A",
    "#66A61E", "#A6761D", "#1F78B4", "#FF7F00", "#33A02C", "#984EA3",
    "#A65628", "#F781BF", "#17BECF", "#BCBD22", "#8C564B", "#4D4D4D",
]


def apply_categorical_colors(adata, *keys):
    """Use stable, high-contrast colors for unordered Scanpy metadata."""
    import colorsys
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
        colors = list(QUALITATIVE_COLORS)
        if len(categories) > len(colors):
            extra_n = len(categories) - len(colors)
            colors.extend(
                "#{:02X}{:02X}{:02X}".format(*(round(channel * 255) for channel in colorsys.hsv_to_rgb(i / max(extra_n, 1), 0.68, 0.78)))
                for i in range(extra_n)
            )
        adata.uns[f"{key}_colors"] = colors[:len(categories)]


def read_table(path: Path):
    return pd.read_csv(path, sep="," if str(path).lower().endswith(".csv") else "\t", dtype=str).fillna("")


def _pick_ortholog_column(table, candidates):
    normalized = {re.sub(r"[^a-z0-9]+", "", str(col).lower()): col for col in table.columns}
    for candidate in candidates:
        hit = normalized.get(re.sub(r"[^a-z0-9]+", "", candidate.lower()))
        if hit is not None:
            return hit
    return None


def map_cross_species_genes(genes, set_names, source_species, ortholog_path: Path, expression_genes, audit_path: Path):
    pairs = [(str(g).strip(), str(s)) for g, s in zip(genes, set_names) if str(g).strip()]
    genes = [g for g, _ in pairs]
    set_names = [s for _, s in pairs]
    if source_species == "same":
        return genes, set_names
    orth = read_table(ortholog_path)
    mouse_col = _pick_ortholog_column(orth, ["mouse_gene_symbol", "mouse_symbol", "mgi_symbol", "marker_symbol", "external_gene_name", "mouse_gene_name"])
    human_col = _pick_ortholog_column(orth, ["human_gene_symbol", "human_symbol", "hgnc_symbol", "human_gene_name", "hsapiens_homolog_associated_gene_name"])
    if mouse_col is None or human_col is None:
        raise SystemExit("Ortholog table needs recognizable mouse and human gene-symbol columns (for example mouse_gene_symbol and human_gene_symbol).")
    orth = orth[[mouse_col, human_col]].rename(columns={mouse_col: "mouse", human_col: "human"})
    orth = orth[(orth["mouse"].str.strip() != "") & (orth["human"].str.strip() != "")].drop_duplicates()
    expression = set(map(str, expression_genes))
    mouse_overlap = len(expression.intersection(orth["mouse"]))
    human_overlap = len(expression.intersection(orth["human"]))
    target_species = "mouse" if mouse_overlap > human_overlap else "human" if human_overlap > mouse_overlap else "unknown"
    if target_species == "unknown":
        raise SystemExit("Could not infer whether expression features are mouse or human from the ortholog table. Choose 'Same as the expression dataset' if conversion is unnecessary.")
    if source_species == "auto":
        expression_lookup = {str(g).upper(): str(g) for g in expression_genes}
        other_species = "human" if target_species == "mouse" else "mouse"
        source_upper = orth[other_species].astype(str).str.upper()
        counts = source_upper.value_counts()
        usable = orth[source_upper.map(counts).eq(1)].copy()
        ortholog_lookup = dict(zip(usable[other_species].astype(str).str.upper(), usable[target_species].astype(str)))
        direct = [expression_lookup.get(g.upper(), "") for g in genes]
        mapped = [d or ortholog_lookup.get(g.upper(), "") for g, d in zip(genes, direct)]
        status = ["auto_case_match" if d else "auto_ortholog" if m else "unmapped_or_ambiguous" for d, m in zip(direct, mapped)]
        audit = pd.DataFrame({"set": set_names, "original_gene": genes, "mapped_gene": mapped, "status": status, "source_species": "auto", "target_species": target_species})
        audit.to_csv(audit_path, sep="\t", index=False)
        keep = audit["status"].isin(["auto_case_match", "auto_ortholog"])
        return audit.loc[keep, "mapped_gene"].tolist(), audit.loc[keep, "set"].tolist()
    if source_species == target_species:
        audit = pd.DataFrame({"set": set_names, "original_gene": genes, "mapped_gene": genes, "status": "same_species", "source_species": source_species, "target_species": target_species})
    else:
        source_col, target_col = source_species, target_species
        counts = orth[source_col].value_counts()
        usable = orth[orth[source_col].map(counts).eq(1)]
        lookup = dict(zip(usable[source_col], usable[target_col]))
        mapped = [lookup.get(g, "") for g in genes]
        audit = pd.DataFrame({"set": set_names, "original_gene": genes, "mapped_gene": mapped, "status": ["mapped" if x else "unmapped_or_ambiguous" for x in mapped], "source_species": source_species, "target_species": target_species})
    audit.to_csv(audit_path, sep="\t", index=False)
    keep = audit["status"].isin(["mapped", "same_species"])
    return audit.loc[keep, "mapped_gene"].tolist(), audit.loc[keep, "set"].tolist()


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
        "batch_column": get("batch_column", "sample_id"),
        "cluster_resolution": float(get("cluster_resolution", "0.6") or 0.6),
        "n_pcs": int(float(get("n_pcs", "30") or 30)),
        "n_neighbors": int(float(get("n_neighbors", "15") or 15)),
        "umap_min_dist": float(get("umap_min_dist", "0.3") or 0.3),
        "umap_spread": float(get("umap_spread", "1.0") or 1.0),
        "umap_metric": get("umap_metric", "euclidean").lower(),
        "umap_init_pos": get("umap_init_pos", "spectral").lower(),
        "min_features": int(float(get("min_features", "200") or 200)),
        "min_counts": int(float(get("min_counts", "0") or 0)),
        "max_features": int(float(get("max_features", "0") or 0)),
        "max_percent_mt": float(get("max_percent_mt", "20") or 20),
        "min_cells_per_gene": int(float(get("min_cells_per_gene", "3") or 3)),
        "doublet_method": get("doublet_method", "auto").lower(),
        # Zero requests a capture-specific automatic estimate. Positive values
        # are explicit global overrides.
        "doublet_rate": float(get("doublet_rate", "0") or 0),
        "remove_doublets": get_bool("remove_doublets", True),
        "marker_file": get("marker_file", ""),
        "celltype_file": get("celltype_file", ""),
        "marker_species": get("marker_species", "auto").lower(),
        "marker_ortholog_file": get("marker_ortholog_file", ""),
        "annotation_name": get("annotation_name", "cell_type"),
        "signature_file": get("signature_file", ""),
        "signature_species": get("signature_species", "same").lower(),
        "signature_ortholog_file": get("signature_ortholog_file", ""),
        "de_group_column": get("de_group_column", "condition"),
        "de_reference": get("de_reference", ""),
        "de_comparison": get("de_comparison", ""),
        "de_annotation_column": get("de_annotation_column", ""),
        "de_annotation_values": [x for x in get("de_annotation_values", get("de_annotation_value", "all")).split("||") if x],
        "de_method": get("de_method", "both").lower(),
        "de_covariates": [x.strip() for x in get("de_covariates", "").split(",") if x.strip()],
        "scvi_max_epochs": int(float(get("scvi_max_epochs", "400") or 400)),
        "harmony_theta": float(get("harmony_theta", "2") or 2),
        "harmony_lambda": float(get("harmony_lambda", "1") or 1),
        "harmony_max_iter": int(float(get("harmony_max_iter", "20") or 20)),
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
    apply_categorical_colors(adata, color)
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
    apply_categorical_colors(adata, color)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sc.pl.umap(adata, color=color, ax=ax, show=False, title=title or str(color), frameon=False, color_map=JP_COLOR_MAP_NAME)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def recommended_pcs_from_variance(variance):
    """Estimate a practical PCA elbow while avoiding unrealistically tiny choices."""
    values = np.asarray(variance, dtype=float).ravel()
    usable = min(len(values), 50)
    if usable < 3 or not np.isfinite(values[:usable]).all():
        return min(30, len(values)) if len(values) else 0
    x = np.linspace(0, 1, usable)
    y = values[:usable]
    y_range = y.max() - y.min()
    if y_range <= 0:
        return min(30, usable)
    y = (y - y.min()) / y_range
    line = y[0] + (y[-1] - y[0]) * x
    elbow = int(np.argmax(np.abs(y - line)) + 1)
    return max(10, min(50, usable, elbow))


def save_pca_outputs(adata, figures: Path, recommended_pcs=None):
    """Save compact, interpretable PCA outputs immediately after preprocessing."""
    import matplotlib.pyplot as plt
    variance = np.asarray(adata.uns.get("pca", {}).get("variance_ratio", [])).ravel()
    if len(variance):
        fig, ax = plt.subplots(figsize=(7.2, 4.6), layout="constrained")
        pcs = np.arange(1, len(variance) + 1)
        ax.plot(pcs, 100 * variance, marker="o", markersize=3.5, linewidth=1.5, color=JP_COLOR_MAP(0.12))
        if recommended_pcs and recommended_pcs <= len(variance):
            ax.axvline(recommended_pcs, color=JP_COLOR_MAP(0.92), linestyle="--", linewidth=1.3, label=f"Suggested elbow: PC{recommended_pcs}")
            ax.legend(frameon=False, loc="upper right")
        ax.set(xlabel="Principal component", ylabel="Variance explained (%)", title="PCA variance explained")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#d1d5db", linewidth=0.7)
        fig.savefig(figures / "03_pca_variance_explained.png", dpi=160, bbox_inches="tight")
        plt.close(fig)
    if "X_pca" not in adata.obsm or adata.obsm["X_pca"].shape[1] < 2:
        return
    sample_series = adata.obs["sample_id"].astype(str)
    samples = list(dict.fromkeys(sample_series.tolist()))
    colors = [JP_COLOR_MAP(position) for position in np.linspace(0.18, 0.82, max(len(samples), 2))]
    fig, ax = plt.subplots(figsize=(7.2, 5.4), layout="constrained")
    coordinates = adata.obsm["X_pca"]
    for sample, color in zip(samples, colors):
        keep = (sample_series == sample).to_numpy()
        ax.scatter(coordinates[keep, 0], coordinates[keep, 1], s=5, alpha=0.35, color=color, edgecolors="none", rasterized=True, label=sample)
    x_pct = 100 * variance[0] if len(variance) else 0
    y_pct = 100 * variance[1] if len(variance) > 1 else 0
    ax.set(xlabel=f"PC1 ({x_pct:.1f}% variance)", ylabel=f"PC2 ({y_pct:.1f}% variance)", title="PCA by input sample")
    ax.spines[["top", "right"]].set_visible(False)
    if len(samples) > 1:
        ax.legend(title="Sample", frameon=False, markerscale=2.2, loc="best")
    fig.savefig(figures / "03_pca_by_sample.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def save_dashboard_expression(adata, tables: Path, max_genes=1000):
    """Save a compact normalized expression layer for interactive UMAP gene coloring."""
    from scipy import sparse
    from scipy.io import mmwrite
    import gzip
    counts = adata.layers.get("counts")
    if counts is None:
        return
    hvgs = np.flatnonzero(np.asarray(adata.var.get("highly_variable", False), dtype=bool))
    if not len(hvgs):
        hvgs = np.arange(adata.n_vars)
    selected = hvgs[:min(max_genes, len(hvgs))]
    expression = counts[:, selected].copy()
    library_size = np.asarray(counts.sum(axis=1)).ravel()
    library_size[library_size <= 0] = 1
    if sparse.issparse(expression):
        expression = expression.multiply((1e4 / library_size)[:, None]).tocsr()
        expression.data = np.log1p(expression.data)
    else:
        expression = np.log1p(np.asarray(expression) * (1e4 / library_size)[:, None])
    with gzip.open(tables / "dashboard_gene_expression.mtx.gz", "wb") as handle:
        mmwrite(handle, expression)
    pd.DataFrame({"cell": adata.obs_names.astype(str)}).to_csv(tables / "dashboard_gene_expression_cells.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": adata.var_names[selected].astype(str)}).to_csv(tables / "dashboard_gene_expression_genes.tsv", sep="\t", index=False)


def save_umap_coordinates(adata, tables: Path):
    """Write the interactive dashboard coordinates as soon as UMAP exists."""
    umap_table = pd.DataFrame(adata.obsm["X_umap"][:, :2], index=adata.obs_names, columns=["UMAP_1", "UMAP_2"])
    umap_metadata = adata.obs.copy()
    if "cell" in umap_metadata.columns:
        umap_metadata = umap_metadata.rename(columns={"cell": "input_cell"})
    umap_table = pd.concat([umap_table, umap_metadata], axis=1)
    umap_table.insert(0, "cell", umap_table.index.astype(str))
    umap_table.to_csv(tables / "umap_coordinates.tsv", sep="\t", index=False)


def qc_recommendations(adata):
    """Suggest conservative, distribution-aware QC thresholds for user review.

    The workflow still requires the user to review and may override these values.
    Per-sample suggestions use robust distribution tails; the global row is
    deliberately conservative so one low-complexity sample is not discarded
    before sample-specific QC can be reviewed.
    """
    def rounded(value, direction, minimum=0):
        """Round lower filters down and upper filters up for readable QC choices."""
        value = float(value)
        increment = 100 if abs(value) >= 100 else 10
        rounded_value = np.floor(value / increment) * increment if direction == "down" else np.ceil(value / increment) * increment
        return int(max(minimum, rounded_value))

    rows = []
    for sample_id, sample in adata.obs.groupby("sample_id", observed=True):
        genes = pd.to_numeric(sample["n_genes_by_counts"], errors="coerce").dropna()
        counts = pd.to_numeric(sample["total_counts"], errors="coerce").dropna()
        mitochondrial = pd.to_numeric(sample["pct_counts_mt"], errors="coerce").dropna()
        if not len(genes) or not len(counts) or not len(mitochondrial):
            continue
        gene_q05, gene_q75, gene_q99 = np.quantile(genes, [0.05, 0.75, 0.99])
        count_q05 = np.quantile(counts, 0.05)
        gene_iqr = gene_q75 - np.quantile(genes, 0.25)
        rows.append({
            "sample_id": str(sample_id),
            "cells_before_qc": int(len(sample)),
            "min_features": rounded(gene_q05, "down", minimum=200),
            "min_counts": rounded(count_q05, "down"),
            "max_features": rounded(max(gene_q99, gene_q75 + 3 * gene_iqr), "up"),
            "max_percent_mt": int(min(25, max(10, np.ceil(np.quantile(mitochondrial, 0.95))))),
            "recommendation_basis": "5th percentile lower filters rounded down (tens below 100; hundreds otherwise); 99th percentile/IQR high-gene screen rounded up on the same scale; 95th percentile mitochondrial screen",
        })
    recommendations = pd.DataFrame(rows)
    if recommendations.empty:
        return recommendations
    global_row = {
        "sample_id": "Recommended global",
        "cells_before_qc": int(recommendations["cells_before_qc"].sum()),
        "min_features": int(recommendations["min_features"].min()),
        "min_counts": int(recommendations["min_counts"].min()),
        "max_features": int(recommendations["max_features"].max()),
        "max_percent_mt": int(recommendations["max_percent_mt"].max()),
        "recommendation_basis": "Conservative global values across samples; review the per-sample rows before running QC",
    }
    return pd.concat([pd.DataFrame([global_row]), recommendations], ignore_index=True)


def save_qc_plots(adata, figures: Path, prefix: str = "01_qc", cutoffs=None, state_label=None):
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
    def cutoff_lines(axis, column):
        if not cutoffs:
            return
        lines = []
        if column == "n_genes_by_counts":
            lines.append((cutoffs["min_features"], "≥ minimum"))
            if cutoffs["max_features"] > 0:
                lines.append((cutoffs["max_features"], "≤ maximum"))
        elif column == "total_counts":
            lines.append((cutoffs["min_counts"], "≥ minimum"))
        elif column == "pct_counts_mt":
            lines.append((cutoffs["max_percent_mt"], "≤ maximum"))
        for value, label in lines:
            axis.axhline(value, color=JP_COLOR_MAP(0.92), linestyle="--", linewidth=1.25, zorder=4)
            axis.text(0.99, value, f" {label}: {value:g}", color=JP_COLOR_MAP(0.98), fontsize=8,
                      ha="right", va="bottom", transform=axis.get_yaxis_transform(),
                      bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5})
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
        cutoff_lines(axis, column)
        if column == "pct_counts_mt":
            values = np.concatenate(groups)
            # A handful of all-mitochondrial droplets can stretch this axis to
            # 100% and hide the biologically useful part of the distribution.
            # Show the central 99%, capped at 50%. The underlying observations,
            # exported tables, and filtering values remain completely uncapped.
            display_max = max(10.0, float(np.nanquantile(values, 0.99)))
            display_max = min(50.0, np.ceil(display_max / 5.0) * 5.0)
            axis.set_ylim(bottom=0, top=display_max)
            axis.text(0.01, 0.98, f"Display cap: {display_max:g}%",
                      transform=axis.transAxes, ha="left", va="top", fontsize=8,
                      color="#4b5563", bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5})
        axis.grid(axis="y", color="#d1d5db", linewidth=0.7, alpha=0.8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(axis="x", length=0)
    if single_sample:
        fig.suptitle(f"QC overview — {state_label or display_labels[0]}", fontsize=14, fontweight="bold")
    elif state_label:
        fig.suptitle(f"QC overview — {state_label}", fontsize=14, fontweight="bold")
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
    if cutoffs:
        ax.axvline(cutoffs["min_counts"], color=JP_COLOR_MAP(0.92), linestyle="--", linewidth=1.25)
        ax.axhline(cutoffs["max_percent_mt"], color=JP_COLOR_MAP(0.92), linestyle="--", linewidth=1.25)
        ax.text(0.99, 0.98, f"Applied cutoffs: counts ≥ {cutoffs['min_counts']:g}; mitochondrial reads ≤ {cutoffs['max_percent_mt']:g}%",
                transform=ax.transAxes, ha="right", va="top", fontsize=8, color=JP_COLOR_MAP(0.98),
                bbox={"facecolor": "white", "edgecolor": "#fecaca", "alpha": 0.9, "pad": 3})
    mt_values = pd.to_numeric(adata.obs["pct_counts_mt"], errors="coerce").dropna().to_numpy()
    if len(mt_values):
        display_max = max(10.0, float(np.nanquantile(mt_values, 0.99)))
        display_max = min(50.0, np.ceil(display_max / 5.0) * 5.0)
        ax.set_ylim(bottom=0, top=display_max)
        ax.text(0.01, 0.98, f"Mitochondrial display cap: {display_max:g}%",
                transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#4b5563",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5})
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
    """Run Scrublet on raw counts, independently for each droplet capture."""
    method = p["doublet_method"]
    if method not in {"auto", "none", "scrublet"}:
        raise SystemExit("Scanpy doublet_method must be auto, none, or scrublet.")
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    report = []
    if method == "none":
        report.append({"capture_id": "all", "sample_ids": ",".join(sorted(adata.obs["sample_id"].astype(str).unique())), "method": "none", "expected_doublet_rate": np.nan if p["doublet_rate"] <= 0 else p["doublet_rate"], "rate_source": "automatic" if p["doublet_rate"] <= 0 else "global override", "cells_before": adata.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": "Doublet removal disabled"})
        pd.DataFrame(report).to_csv(tables / "doublet_summary_by_sample.tsv", sep="\t", index=False)
        pd.DataFrame(report).to_csv(tables / "doublet_summary_by_capture.tsv", sep="\t", index=False)
        calls = adata.obs[["sample_id", "capture_id", "doublet_score", "predicted_doublet"]].copy()
        calls.insert(0, "cell", adata.obs_names)
        calls["removed_as_doublet"] = False
        calls.to_csv(tables / "doublet_calls.tsv", sep="\t", index=False)
        adata.uns["codespring_doublets_removed"] = 0
        return adata
    capture_col = "capture_id" if "capture_id" in adata.obs else "sample_id" if "sample_id" in adata.obs else None
    groups = [(str(name), list(idx)) for name, idx in adata.obs.groupby(capture_col, observed=True).groups.items()] if capture_col else [("all", list(adata.obs_names))]
    for capture_id, cells in groups:
        sub = adata[cells].copy()
        sample_ids = ",".join(sorted(sub.obs["sample_id"].astype(str).unique()))
        expected_rate = p["doublet_rate"] if p["doublet_rate"] > 0 else min(0.15, max(0.008, 0.008 * sub.n_obs / 1000.0))
        rate_source = "global override" if p["doublet_rate"] > 0 else "automatic 10x estimate"
        if "expected_doublet_rate" in sub.obs.columns:
            raw_rates = sub.obs["expected_doublet_rate"].replace({"": np.nan, "nan": np.nan, "NA": np.nan, "None": np.nan})
            rates = pd.to_numeric(raw_rates, errors="coerce").dropna().unique()
            if len(rates) > 1:
                raise SystemExit(f"expected_doublet_rate must have one value per capture ({capture_id}).")
            if len(rates) == 1:
                expected_rate = float(rates[0])
                rate_source = "manifest"
                if not 0 < expected_rate < 1:
                    raise SystemExit(f"expected_doublet_rate must be greater than 0 and less than 1 for {capture_id}.")
        # Scrublet is not informative on very small partitions. Do not invent
        # calls; retain these cells and record why in the summary.
        if sub.n_obs < 100 or sub.n_vars < 20:
            report.append({"capture_id": capture_id, "sample_ids": sample_ids, "method": "scrublet", "expected_doublet_rate": expected_rate, "rate_source": rate_source, "cells_before": sub.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": "Skipped: fewer than 100 cells or 20 genes"})
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
            report.append({"capture_id": capture_id, "sample_ids": sample_ids, "method": "scrublet", "expected_doublet_rate": expected_rate, "rate_source": rate_source, "cells_before": sub.n_obs, "predicted_doublets": predicted_n, "removed_doublets": predicted_n if p["remove_doublets"] else 0, "note": ""})
        except Exception as exc:
            if method == "scrublet":
                raise SystemExit(f"Scrublet failed for {capture_id}: {exc}") from exc
            report.append({"capture_id": capture_id, "sample_ids": sample_ids, "method": "scrublet", "expected_doublet_rate": expected_rate, "rate_source": rate_source, "cells_before": sub.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": f"Automatic Scrublet skipped after error: {exc}"})
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
    save_doublet_plot(adata, figures)
    pd.DataFrame(report).to_csv(tables / "doublet_summary_by_sample.tsv", sep="\t", index=False)
    pd.DataFrame(report).to_csv(tables / "doublet_summary_by_capture.tsv", sep="\t", index=False)
    calls = adata.obs[["sample_id", "capture_id", "doublet_score", "predicted_doublet"]].copy()
    calls.insert(0, "cell", adata.obs_names)
    calls["removed_as_doublet"] = calls["predicted_doublet"] & p["remove_doublets"]
    calls.to_csv(tables / "doublet_calls.tsv", sep="\t", index=False)
    adata.uns["codespring_doublets_removed"] = int(calls["removed_as_doublet"].sum())
    if p["remove_doublets"]:
        adata = adata[~adata.obs["predicted_doublet"]].copy()
    if adata.n_obs == 0:
        raise SystemExit("Doublet removal removed every cell; use a lower expected doublet rate or review the input data.")
    return adata


def safe_metadata_name(value: str, default: str = "cell_type") -> str:
    value = re.sub(r"[^A-Za-z0-9_]+", "_", str(value).strip()).strip("_")
    if not value:
        value = default
    if value[0].isdigit():
        value = "annotation_" + value
    return value


def apply_celltype_mapping(adata, path: Path, annotation_name: str):
    tab = read_table(path)
    cell_col = next((x for x in ["cell", "cell_id", "barcode", "Cell", "CellID"] if x in tab.columns), None)
    type_col = next((x for x in ["cell_type", "celltype", "CellType", "annotation"] if x in tab.columns), None)
    if not cell_col or not type_col:
        raise SystemExit("Cell-type mapping needs a cell/barcode column and a cell_type column.")
    labels = dict(zip(tab[cell_col].astype(str), tab[type_col].astype(str)))
    values = [labels.get(str(cell), labels.get(str(raw), "Unassigned")) for cell, raw in zip(adata.obs_names, adata.obs["source_barcode"])]
    adata.obs[annotation_name] = pd.Categorical(values)
    adata.obs[f"annotation_source__{annotation_name}"] = "provided cell metadata"
    adata.obs["annotation_source"] = f"provided cell metadata ({annotation_name})"


def apply_marker_annotation(adata, path: Path, tables_dir: Path, annotation_name: str, marker_species="same", ortholog_path=None):
    tab = read_table(path)
    if not {"cell_type", "gene"}.issubset(tab.columns):
        raise SystemExit("Marker list needs cell_type and gene columns.")
    expanded_genes, expanded_sets = [], []
    for cell_type, group in tab.groupby("cell_type", sort=False):
        for value in group["gene"].astype(str):
            parsed = [g for g in re.split(r"[,;\s]+", value) if g]
            expanded_genes.extend(parsed); expanded_sets.extend([str(cell_type)] * len(parsed))
    if marker_species != "same":
        expression_genes = adata.raw.var_names if adata.raw is not None else adata.var_names
        expanded_genes, expanded_sets = map_cross_species_genes(expanded_genes, expanded_sets, marker_species, Path(ortholog_path).expanduser(), expression_genes, tables_dir / "marker_ortholog_mapping.tsv")
    expanded = pd.DataFrame({"cell_type": expanded_sets, "gene": expanded_genes})
    lists = {}
    genes = set(adata.raw.var_names if adata.raw is not None else adata.var_names)
    for cell_type, group in expanded.groupby("cell_type", sort=False):
        entries = group["gene"].astype(str).tolist()
        keep = sorted(set(entries).intersection(genes))
        if keep:
            lists[str(cell_type)] = keep
    if not lists:
        raise SystemExit("None of the supplied marker genes were present in the expression object.")
    # Store the exact usable genes so plotted marker panels match the genes
    # that actually contributed to the annotation score.
    adata.uns["codespring_marker_gene_sets"] = {name: list(genes_for_type) for name, genes_for_type in lists.items()}
    for name, genes_for_type in lists.items():
        # Use normalized log-expression, not the scaled/clipped matrix used
        # for PCA and neighbor construction.
        sc.tl.score_genes(adata, genes_for_type, score_name=f"marker_score__{name}", use_raw=True)
    clusters = adata.obs["cluster"].astype(str)
    score_columns = [f"marker_score__{x}" for x in lists]
    means = adata.obs.assign(cluster=clusters).groupby("cluster", observed=True)[score_columns].mean()
    winners = means.idxmax(axis=1).str.replace("marker_score__", "", regex=False)
    adata.obs[annotation_name] = pd.Categorical(clusters.map(winners).fillna("Unassigned"))
    adata.obs[f"annotation_source__{annotation_name}"] = "marker-list cluster scoring"
    adata.obs["annotation_source"] = f"marker-list cluster scoring ({annotation_name})"
    means.to_csv(tables_dir / f"marker_annotation_cluster_scores__{annotation_name}.tsv", sep="\t")
    means.to_csv(tables_dir / "marker_annotation_cluster_scores.tsv", sep="\t")


def marker_list_panels(marker_sets, max_cell_types=6, max_genes=32):
    panels, current, current_genes = [], {}, 0
    for label, genes in marker_sets.items():
        unique_genes = list(dict.fromkeys(str(gene) for gene in genes if str(gene)))
        blocks = [unique_genes[i:i + max_genes] for i in range(0, len(unique_genes), max_genes)]
        for index, block in enumerate(blocks, start=1):
            block_label = f"{label} (part {index})" if len(blocks) > 1 else str(label)
            if current and (len(current) >= max_cell_types or current_genes + len(block) > max_genes):
                panels.append(current); current, current_genes = {}, 0
            current[block_label] = block; current_genes += len(block)
    if current:
        panels.append(current)
    return panels


def save_marker_annotation_panels(adata, marker_sets, figures: Path, annotation_name="cell_type"):
    if not marker_sets:
        return
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from matplotlib.lines import Line2D
    expression = adata.raw if adata.raw is not None else adata
    available = set(expression.var_names)
    for panel_index, panel in enumerate(marker_list_panels(marker_sets), start=1):
        marker_rows, seen = [], set()
        for marker_set, values in panel.items():
            for gene in values:
                if gene in available and gene not in seen:
                    marker_rows.append((marker_set, gene)); seen.add(gene)
        genes = [gene for _, gene in marker_rows]
        if not genes:
            continue
        row_labels = [f"{marker_set} | {gene}" for marker_set, gene in marker_rows]
        row_breaks = [index - 0.5 for index in range(1, len(marker_rows)) if marker_rows[index][0] != marker_rows[index - 1][0]]
        title = "Marker-list annotation: " + "; ".join(panel.keys())
        groupings = [("cluster", "Cluster")]
        if annotation_name in adata.obs.columns:
            groupings.append((annotation_name, "Cell type"))
        for group_column, group_label in groupings:
            groups = adata.obs[group_column].astype(str)
            group_levels = sorted(groups.unique(), key=lambda value: (not str(value).replace(".", "", 1).isdigit(), float(value) if group_column == "cluster" and str(value).replace(".", "", 1).isdigit() else str(value)))
            means = np.zeros((len(genes), len(group_levels)), dtype=float); fractions = np.zeros((len(genes), len(group_levels)), dtype=float)
            for group_index, group in enumerate(group_levels):
                subset = expression[groups.to_numpy() == group, genes]
                values = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
                means[:, group_index] = np.asarray(values.mean(axis=0)).ravel(); fractions[:, group_index] = np.asarray((values > 0).mean(axis=0)).ravel()
            z_scores = np.clip((means - means.mean(axis=1, keepdims=True)) / np.maximum(means.std(axis=1, keepdims=True), 1e-8), -2.5, 2.5)
            width, height = max(8.0, 4.4 + 0.62 * len(group_levels)), max(5.5, 2.5 + 0.24 * len(genes))
            fig, ax = plt.subplots(figsize=(width, height))
            x_coords, y_coords = np.meshgrid(np.arange(len(group_levels)), np.arange(len(genes)))
            dots = ax.scatter(x_coords.ravel(), y_coords.ravel(), s=12 + 260 * fractions.ravel(), c=z_scores.ravel(), cmap="RdBu_r", vmin=-2.5, vmax=2.5, edgecolors="#374151", linewidths=0.18)
            ax.set_xticks(np.arange(len(group_levels)), group_levels, rotation=45, ha="right"); ax.set_yticks(np.arange(len(genes)), row_labels)
            for row_break in row_breaks: ax.axhline(row_break, color="#4B5563", linewidth=0.7)
            ax.invert_yaxis(); ax.set_xlabel(group_label); ax.set_ylabel("Marker set | gene")
            ax.set_title(f"Marker expression by {group_label}", fontweight="bold", loc="left", fontsize=13)
            ax.text(0, 1.01, title + ". Color is scaled within each gene; size is the fraction of positive cells.", transform=ax.transAxes, fontsize=9, va="bottom")
            ax.xaxis.grid(True, color="#E5E7EB", linewidth=0.45); ax.set_axisbelow(True)
            colorbar = fig.colorbar(dots, ax=ax, pad=0.02); colorbar.set_label("Relative mean expression\n(row Z-score)")
            size_handles = [Line2D([], [], marker="o", linestyle="", markerfacecolor="#D1D5DB", markeredgecolor="#374151", markersize=np.sqrt(12 + 260 * value / 100), label=f"{value}%") for value in (0, 25, 50, 75, 100)]
            ax.legend(handles=size_handles, title="Cells expressing", bbox_to_anchor=(1.18, 0.33), loc="center left", frameon=False)
            fig.tight_layout(); fig.savefig(figures / f"07_marker_annotation_dotplot_by_{'cluster' if group_column == 'cluster' else 'cell_type'}_panel_{panel_index:02d}.png", dpi=220, bbox_inches="tight"); plt.close(fig)
            fig, ax = plt.subplots(figsize=(width, height))
            image = ax.imshow(z_scores, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5))
            ax.set_xticks(np.arange(len(group_levels)), group_levels, rotation=45, ha="right"); ax.set_yticks(np.arange(len(genes)), row_labels)
            for row_break in row_breaks: ax.axhline(row_break, color="#4B5563", linewidth=0.7)
            ax.set_xlabel(group_label); ax.set_ylabel("Marker set | gene"); ax.set_title(f"Marker expression by {group_label}", fontweight="bold", loc="left", fontsize=13)
            ax.text(0, 1.01, title + ". Mean normalized expression, scaled within each gene.", transform=ax.transAxes, fontsize=9, va="bottom")
            colorbar = fig.colorbar(image, ax=ax, pad=0.02); colorbar.set_label("Row Z-score")
            fig.tight_layout(); fig.savefig(figures / f"07_marker_annotation_heatmap_by_{'cluster' if group_column == 'cluster' else 'cell_type'}_panel_{panel_index:02d}.png", dpi=220, bbox_inches="tight"); plt.close(fig)


def save_cluster_marker_heatmaps(adata, markers, figures: Path, genes_per_cluster=10, clusters_per_panel=4):
    if markers.empty or "group" not in markers.columns or "names" not in markers.columns:
        return
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    ranked = []
    for cluster, group in markers.groupby("group", sort=False, observed=True):
        if "pvals_adj" in group.columns:
            group = group.sort_values("pvals_adj", kind="stable")
        ranked.append(group.head(genes_per_cluster).assign(source_cluster=str(cluster)))
    top = pd.concat(ranked, ignore_index=True) if ranked else pd.DataFrame()
    if top.empty:
        return
    cluster_levels = sorted(adata.obs["cluster"].astype(str).unique(), key=lambda value: (not str(value).replace(".", "", 1).isdigit(), float(value) if str(value).replace(".", "", 1).isdigit() else str(value)))
    expression = adata.raw if adata.raw is not None else adata
    available = set(expression.var_names)
    source_clusters = [str(value) for value in top["source_cluster"].drop_duplicates()]
    panels = [source_clusters[index:index + clusters_per_panel] for index in range(0, len(source_clusters), clusters_per_panel)]
    group_labels = adata.obs["cluster"].astype(str).to_numpy()
    for panel_index, panel_clusters in enumerate(panels, start=1):
        panel_top = top[top["source_cluster"].isin(panel_clusters) & top["names"].astype(str).isin(available)].copy()
        if panel_top.empty:
            continue
        genes = panel_top["names"].astype(str).tolist()
        raw_labels = [f"{gene}  [cluster {cluster}]" for gene, cluster in zip(genes, panel_top["source_cluster"].astype(str))]
        seen, labels = {}, []
        for label in raw_labels:
            seen[label] = seen.get(label, 0) + 1
            labels.append(label if seen[label] == 1 else f"{label} ({seen[label]})")
        means = np.zeros((len(genes), len(cluster_levels)), dtype=float)
        for cluster_index, cluster in enumerate(cluster_levels):
            subset = expression[group_labels == cluster, genes]
            values = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
            means[:, cluster_index] = np.asarray(values.mean(axis=0)).ravel()
        z_scores = (means - means.mean(axis=1, keepdims=True)) / np.maximum(means.std(axis=1, keepdims=True), 1e-8)
        z_scores = np.clip(z_scores, -2.5, 2.5)
        width, height = max(8.0, 3.8 + 0.58 * len(cluster_levels)), max(6.0, 2.8 + 0.22 * len(genes))
        fig, ax = plt.subplots(figsize=(width, height))
        image = ax.imshow(z_scores, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5))
        ax.set_xticks(np.arange(len(cluster_levels)), cluster_levels, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(labels)), labels)
        ax.set_xlabel("Cluster"); ax.set_ylabel("Marker gene [source cluster]")
        ax.set_title("Top cluster markers: " + ", ".join(panel_clusters) + "\nTop 10 ranked markers per selected cluster; mean normalized expression, row-scaled", fontweight="bold", loc="left", fontsize=11)
        colorbar = fig.colorbar(image, ax=ax, pad=0.02); colorbar.set_label("Row Z-score")
        fig.tight_layout(); fig.savefig(figures / f"08_cluster_marker_heatmap_panel_{panel_index:02d}.png", dpi=180, bbox_inches="tight"); plt.close(fig)


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


def write_cell_metadata_outputs(adata, tables_dir: Path):
    cell_metadata = adata.obs.copy()
    if "cell" in cell_metadata.columns:
        cell_metadata = cell_metadata.rename(columns={"cell": "input_cell"})
    cell_metadata.insert(0, "cell", cell_metadata.index.astype(str))
    cell_metadata.to_csv(tables_dir / "cell_metadata.tsv", sep="\t", index=False)
    save_umap_coordinates(adata, tables_dir)


def write_composition_outputs(adata, annotation_name: str, tables_dir: Path):
    if annotation_name not in adata.obs.columns:
        raise SystemExit(f"Annotation metadata field was not created: {annotation_name}")
    current = adata.obs.groupby(["sample_id", annotation_name], observed=True).size().reset_index(name="cells")
    current = current.rename(columns={annotation_name: "annotation_label"})
    current.insert(1, "annotation_field", annotation_name)
    current["proportion_within_sample"] = current["cells"] / current.groupby("sample_id", observed=True)["cells"].transform("sum")
    long_path = tables_dir / "composition_by_sample.tsv"
    if long_path.exists():
        previous = pd.read_csv(long_path, sep="\t")
        if "annotation_field" in previous.columns:
            previous = previous[previous["annotation_field"].astype(str) != annotation_name]
            current = pd.concat([previous, current], ignore_index=True)
    current.to_csv(long_path, sep="\t", index=False)
    legacy = current[current["annotation_field"] == annotation_name].rename(columns={"annotation_label": "cell_type"})
    legacy[["sample_id", "cell_type", "cells", "proportion_within_sample"]].to_csv(tables_dir / "cell_type_by_sample.tsv", sep="\t", index=False)
    cluster_sizes = adata.obs.groupby(["cluster", annotation_name], observed=True).size().reset_index(name="cells")
    cluster_sizes = cluster_sizes.rename(columns={annotation_name: "annotation_label"})
    cluster_sizes.insert(1, "annotation_field", annotation_name)
    cluster_sizes.to_csv(tables_dir / "cluster_annotation_sizes.tsv", sep="\t", index=False)
    cluster_sizes.rename(columns={"annotation_label": "cell_type"})[["cluster", "cell_type", "cells"]].to_csv(tables_dir / "cluster_cell_type_sizes.tsv", sep="\t", index=False)


def score_signatures(adata, path: Path, tables_dir: Path, annotation_name: str, signature_species="same", ortholog_path=None):
    signatures = read_table(path)
    if not {"signature", "gene"}.issubset(signatures.columns):
        raise SystemExit("Signature list needs signature and gene columns.")
    available = set(adata.raw.var_names if adata.raw is not None else adata.var_names)
    expanded_genes, expanded_sets = [], []
    for signature, group in signatures.groupby("signature", sort=False):
        for value in group["gene"].astype(str):
            parsed = [x for x in re.split(r"[,;\s]+", value) if x]
            expanded_genes.extend(parsed); expanded_sets.extend([str(signature)] * len(parsed))
    if signature_species != "same":
        expanded_genes, expanded_sets = map_cross_species_genes(expanded_genes, expanded_sets, signature_species, Path(ortholog_path).expanduser(), available, tables_dir / "signature_ortholog_mapping.tsv")
    signatures = pd.DataFrame({"signature": expanded_sets, "gene": expanded_genes})
    score_columns = []
    coverage = []
    for signature, group in signatures.groupby("signature", sort=False):
        requested = []
        for value in group["gene"].astype(str):
            requested.extend([x for x in re.split(r"[,;\s]+", value) if x])
        requested = list(dict.fromkeys(requested))
        present = [x for x in requested if x in available]
        if not present:
            coverage.append({"signature": signature, "metadata_column": "", "genes_requested": len(requested), "genes_found": 0, "coverage": 0})
            continue
        column = "signature__" + safe_metadata_name(signature, "score")
        sc.tl.score_genes(adata, present, score_name=column, use_raw=True, random_state=1234)
        score_columns.append(column)
        coverage.append({"signature": signature, "metadata_column": column, "genes_requested": len(requested), "genes_found": len(present), "coverage": len(present) / max(1, len(requested))})
    pd.DataFrame(coverage).to_csv(tables_dir / "signature_gene_coverage.tsv", sep="\t", index=False)
    if not score_columns:
        raise SystemExit("None of the supplied signature genes were found in the normalized expression object.")
    per_cell = adata.obs[["sample_id"] + ([annotation_name] if annotation_name in adata.obs.columns else []) + score_columns].copy()
    per_cell.insert(0, "cell", adata.obs_names.astype(str))
    per_cell.to_csv(tables_dir / "signature_scores_per_cell.tsv", sep="\t", index=False)
    group_columns = ["sample_id"] + ([annotation_name] if annotation_name in adata.obs.columns else [])
    summary = adata.obs.groupby(group_columns, observed=True)[score_columns].agg(["mean", "median", "count"])
    summary.columns = [f"{score}__{stat}" for score, stat in summary.columns]
    summary.reset_index().to_csv(tables_dir / "signature_scores_summary.tsv", sep="\t", index=False)
    return score_columns


def export_differential_inputs(adata, p, tables_dir: Path, out_dir: Path):
    group_column = p["de_group_column"]
    reference, comparison = p["de_reference"], p["de_comparison"]
    if group_column not in adata.obs.columns:
        raise SystemExit(f"Differential-expression group field is absent from cell metadata: {group_column}")
    if not reference or not comparison or reference == comparison:
        raise SystemExit("Choose two different reference and comparison groups for differential expression.")
    annotation_column = p["de_annotation_column"]
    targets = p["de_annotation_values"] or ["all"]
    targets = list(dict.fromkeys(targets))
    if any(str(x).lower() not in {"all", "all cells"} for x in targets):
        if not annotation_column or annotation_column not in adata.obs.columns:
            raise SystemExit("Choose an annotation metadata field before requesting individual-population differential expression.")
    contrast_slug = f"{safe_metadata_name(comparison, 'comparison')}_vs_{safe_metadata_name(reference, 'reference')}"
    root = out_dir / "differential_expression" / contrast_slug
    manifest_rows, completed = [], 0
    for annotation_value in targets:
        is_global = str(annotation_value).lower() in {"all", "all cells"}
        population_slug = "global" if is_global else f"{safe_metadata_name(annotation_column, 'annotation')}__{safe_metadata_name(annotation_value, 'population')}"
        job_dir = root / population_slug
        job_dir.mkdir(parents=True, exist_ok=True)
        file_slug = f"{contrast_slug}__{population_slug}"
        row = {
            "job_id": file_slug, "comparison": comparison, "reference": reference,
            "annotation_column": "" if is_global else annotation_column,
            "annotation_value": "all" if is_global else str(annotation_value),
            "output_dir": str(job_dir), "file_slug": file_slug, "status": "prepared", "note": "",
            "count_file": "", "design_file": "", "cell_result_file": "",
        }
        try:
            keep = adata.obs[group_column].astype(str).isin([reference, comparison])
            if not is_global:
                keep &= adata.obs[annotation_column].astype(str) == str(annotation_value)
            subset = adata[keep].copy()
            if subset.n_obs == 0:
                raise ValueError("No cells remain for this population.")
            if not {reference, comparison}.issubset(set(subset.obs[group_column].astype(str))):
                raise ValueError("Both comparison groups must contain cells in this population.")
            if p["de_method"] in {"both", "pseudobulk"}:
                counts = subset.layers.get("counts")
                if counts is None:
                    raise ValueError("Raw counts are unavailable in layers['counts'] for pseudobulk DE.")
                design_rows, matrices = [], []
                for sample_id, cell_index in subset.obs.groupby("sample_id", observed=True).groups.items():
                    positions = subset.obs_names.get_indexer(cell_index)
                    sample_groups = subset.obs.iloc[positions][group_column].astype(str).unique()
                    if len(sample_groups) != 1:
                        raise ValueError(f"Biological replicate {sample_id} contains multiple values of {group_column}.")
                    matrices.append(np.asarray(counts[positions].sum(axis=0)).ravel())
                    design_row = {"sample_id": str(sample_id), "group": sample_groups[0], "cells": len(positions)}
                    for covariate in p["de_covariates"]:
                        if covariate not in subset.obs.columns:
                            raise ValueError(f"Requested covariate is absent: {covariate}")
                        values = subset.obs.iloc[positions][covariate].astype(str).unique()
                        if len(values) != 1:
                            raise ValueError(f"Covariate {covariate} is not constant within replicate {sample_id}.")
                        design_row[covariate] = values[0]
                    design_rows.append(design_row)
                count_file, design_file = job_dir / f"pseudobulk_counts__{file_slug}.tsv", job_dir / f"pseudobulk_design__{file_slug}.tsv"
                count_table = pd.DataFrame(np.vstack(matrices).T, index=subset.var_names, columns=[x["sample_id"] for x in design_rows])
                count_table.insert(0, "gene", count_table.index.astype(str)); count_table.to_csv(count_file, sep="\t", index=False)
                pd.DataFrame(design_rows).to_csv(design_file, sep="\t", index=False)
                row["count_file"], row["design_file"] = str(count_file), str(design_file)
            if p["de_method"] in {"both", "cell", "cell_level"}:
                sc.tl.rank_genes_groups(subset, groupby=group_column, groups=[comparison], reference=reference, method="wilcoxon", use_raw=True, pts=True)
                result = sc.get.rank_genes_groups_df(subset, group=comparison)
                result.insert(0, "comparison", f"{comparison}_vs_{reference}")
                result.insert(1, "population", "global" if is_global else str(annotation_value))
                result["analysis_level"] = "cell-level exploratory Wilcoxon; cells are not biological replicates"
                cell_file = job_dir / f"cell_level_Wilcoxon__{file_slug}.tsv"
                result.to_csv(cell_file, sep="\t", index=False); row["cell_result_file"] = str(cell_file)
                if is_global:
                    result.to_csv(tables_dir / "cell_level_differential_expression.tsv", sep="\t", index=False)
            completed += 1
        except Exception as exc:
            row["status"], row["note"] = "skipped", str(exc)
        manifest_rows.append(row)
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(tables_dir / "differential_jobs.tsv", sep="\t", index=False)
    if completed == 0:
        raise SystemExit("No requested differential-expression population could be prepared. Review tables/differential_jobs.tsv.")


def main():
    if len(sys.argv) not in {4, 5}:
        raise SystemExit("Usage: scrna_pipeline_scanpy.py <samples.tsv> <out_dir> <params.tsv> [inspect|qc|preprocess|cluster|annotate|score|differential|all]")
    samples_path, out_dir, params_path = map(Path, sys.argv[1:4])
    stage = sys.argv[4].lower() if len(sys.argv) == 5 else "all"
    stages = {"inspect", "qc", "preprocess", "cluster", "annotate", "score", "differential", "all"}
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
    if not 2 <= p["n_neighbors"] <= 200:
        raise SystemExit("UMAP neighbours must be between 2 and 200.")
    if not 0 <= p["umap_min_dist"] <= 2 or not 0.1 <= p["umap_spread"] <= 10:
        raise SystemExit("UMAP minimum distance must be 0–2 and spread must be 0.1–10.")
    if p["umap_metric"] not in {"euclidean", "cosine", "manhattan", "correlation"}:
        raise SystemExit("Unsupported UMAP distance metric.")
    if p["umap_init_pos"] not in {"spectral", "random"}:
        raise SystemExit("UMAP initialization must be spectral or random.")
    np.random.seed(p["seed"])
    input_checkpoint = checkpoints / "01_input_scanpy.h5ad"
    qc_checkpoint = checkpoints / "02_qc_scanpy.h5ad"
    preprocess_checkpoint = checkpoints / "03_preprocessed_scanpy.h5ad"
    cluster_checkpoint = checkpoints / "04_clustered_scanpy.h5ad"
    processed_object = objects / "processed_scanpy.h5ad"

    def mark_complete(name):
        (out_dir / f"_STAGE_{name.upper()}_COMPLETE").write_text("complete\n")

    def require_checkpoint(path: Path, prior: str):
        if not path.exists():
            raise SystemExit(f"The {prior} stage has not completed. Run {prior} before this stage.")
        return sc.read_h5ad(path)

    def discard_checkpoint(path: Path):
        """Remove only a checkpoint that has been superseded successfully."""
        try:
            path.unlink()
        except FileNotFoundError:
            pass

    # Validate the manifest once.  Each submitted stage then opens just the
    # one checkpoint it needs; separate SLURM jobs cannot share memory, but
    # they also must not serially reload older checkpoints on the way there.
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
    if "capture_id" not in samples.columns:
        samples["capture_id"] = samples["sample_id"]
    samples["capture_id"] = samples["capture_id"].fillna("").astype(str).str.strip()
    samples.loc[samples["capture_id"] == "", "capture_id"] = samples.loc[samples["capture_id"] == "", "sample_id"]
    if "expected_doublet_rate" in samples.columns:
        raw_rates = samples["expected_doublet_rate"].replace({"": np.nan, "nan": np.nan, "NA": np.nan, "None": np.nan})
        rates = pd.to_numeric(raw_rates, errors="coerce")
        supplied = raw_rates.notna()
        if ((supplied & (rates.isna() | (rates <= 0) | (rates >= 1))).any()):
            raise SystemExit("expected_doublet_rate must be blank for automatic estimation or greater than 0 and less than 1.")
        for capture_id, values in rates[supplied].groupby(samples.loc[supplied, "capture_id"]):
            if values.nunique() > 1:
                raise SystemExit(f"expected_doublet_rate must be constant within capture_id: {capture_id}")

    # Stage 1: inspect each source object and preserve a raw-count checkpoint.
    if stage in {"inspect", "all"}:
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
        recommendations = qc_recommendations(adata)
        recommendations.to_csv(tables / "qc_recommended_thresholds.tsv", sep="\t", index=False)
        global_recommendation = recommendations[recommendations["sample_id"] == "Recommended global"]
        suggested_cutoffs = global_recommendation.iloc[0][["min_features", "min_counts", "max_features", "max_percent_mt"]].to_dict() if len(global_recommendation) else None
        save_qc_plots(adata, figures, prefix="00_qc_pre_filter", cutoffs=suggested_cutoffs, state_label="Unfiltered cells — suggested cutoffs")
        # Continue a single processed AnnData object without rebuilding a
        # valid embedding. This writes only a project-local copy; the source
        # H5AD remains read-only and unchanged.
        detected = pd.DataFrame([item[1] for item in inputs])
        reusable_processed_input = (
            len(samples) == 1
            and str(detected.iloc[0].get("input_kind", "")) == "scanpy_h5ad"
            and bool(detected.iloc[0].get("umap_detected", False))
            and bool(detected.iloc[0].get("clusters_detected", False))
        )
        if reusable_processed_input:
            if "cluster" not in adata.obs.columns:
                for cluster_source in ("leiden", "louvain", "seurat_clusters"):
                    if cluster_source in adata.obs.columns:
                        adata.obs["cluster"] = adata.obs[cluster_source].astype(str)
                        break
            adata.uns["codespring_integration"] = "existing input object"
            adata.uns["codespring_doublets_removed"] = 0
            write_h5ad_checkpoint(adata, processed_object)
        write_h5ad_checkpoint(adata, input_checkpoint)
        mark_complete("inspect")
        if stage == "inspect":
            return
    elif stage == "qc":
        adata = require_checkpoint(input_checkpoint, "input inspection")
    elif stage == "preprocess":
        adata = require_checkpoint(qc_checkpoint, "QC and doublet handling")
    elif stage == "cluster":
        adata = require_checkpoint(preprocess_checkpoint, "normalization and PCA")
    elif stage == "annotate":
        annotation_input = cluster_checkpoint if cluster_checkpoint.exists() else processed_object
        adata = require_checkpoint(annotation_input, "UMAP and clustering")
    elif stage == "score":
        score_input = processed_object if processed_object.exists() and (not cluster_checkpoint.exists() or processed_object.stat().st_mtime >= cluster_checkpoint.stat().st_mtime) else cluster_checkpoint
        adata = require_checkpoint(score_input, "UMAP and clustering")
    elif stage == "differential":
        adata = require_checkpoint(processed_object, "annotation")

    # Stage 2: minimally prefilter, call doublets per capture, then perform the
    # complete user-reviewed QC while the data are still raw.
    if stage in {"qc", "all"}:
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        # Replace the inspection preview with the exact same unfiltered cells,
        # now annotated with the cutoffs the user chose for this QC run.
        save_qc_plots(adata, figures, prefix="00_qc_pre_filter", cutoffs=p, state_label="Before QC filtering")
        pd.DataFrame([{
            "min_features": p["min_features"], "min_counts": p["min_counts"],
            "max_features": p["max_features"], "max_percent_mt": p["max_percent_mt"],
            "remove_predicted_doublets": p["remove_doublets"],
        }]).to_csv(tables / "qc_cutoffs_applied.tsv", sep="\t", index=False)
        cells_before_qc = adata.obs.groupby("sample_id", observed=True).size().rename("cells_input")
        minimal_keep = adata.obs["total_counts"] >= 200
        adata = adata[minimal_keep].copy()
        if adata.n_obs == 0:
            raise SystemExit("The minimal 200-count prefilter removed every cell.")
        adata = run_scrublet(adata, p, tables, figures)
        keep = (adata.obs["n_genes_by_counts"] >= p["min_features"]) & (adata.obs["total_counts"] >= p["min_counts"]) & (adata.obs["pct_counts_mt"] <= p["max_percent_mt"])
        if p["max_features"] > 0:
            keep &= adata.obs["n_genes_by_counts"] <= p["max_features"]
        adata = adata[keep].copy()
        if adata.n_obs == 0:
            raise SystemExit("QC thresholds removed every cell.")
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
        qc = adata.obs[["sample_id", "capture_id", "n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_score", "predicted_doublet"]].copy()
        qc.insert(0, "cell", adata.obs_names)
        qc.to_csv(tables / "qc_cell_metrics.tsv", sep="\t", index=False)
        qc_summary = qc.groupby("sample_id", observed=True).agg(cells_after_qc_and_doublets=("total_counts", "size"), median_umis=("total_counts", "median"), median_genes=("n_genes_by_counts", "median"), median_percent_mt=("pct_counts_mt", "median")).reset_index()
        qc_summary.insert(1, "cells_input", qc_summary["sample_id"].map(cells_before_qc).astype(int))
        qc_summary.to_csv(tables / "qc_summary_by_sample.tsv", sep="\t", index=False)
        save_qc_plots(adata, figures, prefix="01_qc_post_filter", cutoffs=p, state_label="After QC filtering")
        write_h5ad_checkpoint(adata, qc_checkpoint)
        mark_complete("qc")
        if stage == "qc":
            return

    # Stage 3: normalization, highly variable features, scaling, and PCA.
    if stage in {"preprocess", "all"}:
        # Keep only one raw-count matrix in the QC checkpoint.  The counts
        # layer is needed downstream, so create it here immediately before X
        # is normalized instead of storing a second full sparse matrix during
        # the QC review stage.
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        hvg_batch_key = p["batch_column"] if p["batch_column"] in adata.obs.columns and adata.obs[p["batch_column"]].astype(str).nunique() > 1 else None
        hvg_kwargs = {"n_top_genes": min(3000, adata.n_vars)}
        if hvg_batch_key is not None:
            hvg_kwargs["batch_key"] = hvg_batch_key
        # Official Scanpy-style HVG selection on log-normalized expression is
        # batch-aware and needs no optional compiled dependency.
        sc.pp.highly_variable_genes(adata, flavor="seurat", **hvg_kwargs)
        hvg_table = adata.var.copy(); hvg_table.insert(0, "gene", hvg_table.index.astype(str))
        hvg_columns = [c for c in ["gene", "highly_variable", "means", "variances", "variances_norm"] if c in hvg_table.columns]
        hvg_table.loc[:, hvg_columns].to_csv(tables / "highly_variable_genes.tsv", sep="\t", index=False)
        hvg_count = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else 0
        use_hvg = hvg_count >= 2
        pca_features = hvg_count if use_hvg else adata.n_vars
        max_pcs = max(1, min(adata.n_obs - 1, pca_features - 1))
        n_pcs = min(p["n_pcs"], max_pcs)
        # Scanpy's sparse PCA supports implicit zero-centering, so this is
        # standard centered PCA without densifying the full expression matrix.
        sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", use_highly_variable=use_hvg, zero_center=True, random_state=p["seed"])
        pca_ratio = np.asarray(adata.uns["pca"]["variance_ratio"]).ravel()
        recommended_pcs = recommended_pcs_from_variance(pca_ratio)
        pd.DataFrame({"PC": np.arange(1, len(pca_ratio) + 1), "variance_explained": pca_ratio, "percent_variance_explained": 100 * pca_ratio}).to_csv(tables / "pca_variance_explained.tsv", sep="\t", index=False)
        pd.DataFrame([{"recommended_n_pcs": recommended_pcs, "basis": "PCA variance elbow (bounded to 10–50 PCs)"}]).to_csv(tables / "pca_recommended_parameters.tsv", sep="\t", index=False)
        save_pca_outputs(adata, figures, recommended_pcs=recommended_pcs)
        # The dashboard uses this full symbol list to request one gene at a
        # time from the post-UMAP H5AD's normalized `.raw` layer.
        pd.DataFrame({"gene": adata.raw.var_names.astype(str)}).to_csv(tables / "dashboard_all_genes.tsv", sep="\t", index=False)
        # Diagnostic embedding before any technical-batch correction. Keep a
        # separate copy so the integrated UMAP cannot overwrite this view.
        pre_n_pcs = min(p["n_pcs"], adata.obsm["X_pca"].shape[1])
        sc.pp.neighbors(adata, n_neighbors=min(p["n_neighbors"], max(2, adata.n_obs - 1)), n_pcs=pre_n_pcs, use_rep="X_pca", metric=p["umap_metric"])
        sc.tl.umap(adata, min_dist=p["umap_min_dist"], spread=p["umap_spread"], init_pos=p["umap_init_pos"], random_state=p["seed"])
        adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()
        pre = pd.DataFrame(adata.obsm["X_umap_unintegrated"][:, :2], index=adata.obs_names, columns=["UMAP_1", "UMAP_2"])
        pre_metadata = adata.obs.copy()
        if "cell" in pre_metadata.columns:
            pre_metadata = pre_metadata.rename(columns={"cell": "input_cell"})
        pre = pd.concat([pre, pre_metadata], axis=1)
        pre.insert(0, "cell", pre.index.astype(str))
        pre.to_csv(tables / "preintegration_umap_coordinates.tsv", sep="\t", index=False)
        save_umap(adata, "sample_id", figures / "02_preintegration_umap_sample.png", title="Before integration — sample")
        if p["batch_column"] in adata.obs.columns and adata.obs[p["batch_column"]].astype(str).nunique() > 1 and p["batch_column"] != "sample_id":
            save_umap(adata, p["batch_column"], figures / "02_preintegration_umap_batch.png", title=f"Before integration — {p['batch_column']}")
        write_h5ad_checkpoint(adata, preprocess_checkpoint)
        # Once normalized/PCA data are safely checkpointed, the QC object is
        # redundant.  The raw input checkpoint is retained so users can
        # revise cutoffs without rereading the original 10x folders.
        discard_checkpoint(qc_checkpoint)
        mark_complete("preprocess")
        if stage == "preprocess":
            return

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
                sce.pp.harmony_integrate(
                    adata, key=batch_key, basis="X_pca", adjusted_basis="X_harmony",
                    theta=p["harmony_theta"], lamb=p["harmony_lambda"],
                    max_iter_harmony=p["harmony_max_iter"], random_state=p["seed"],
                )
            except Exception as exc:
                raise SystemExit("Harmony integration was requested but harmonypy was unavailable or failed: " + str(exc)) from exc
            representation = "X_harmony"
        sc.pp.neighbors(adata, n_neighbors=min(p["n_neighbors"], max(2, adata.n_obs - 1)), n_pcs=min(n_pcs, adata.obsm[representation].shape[1]), use_rep=representation, metric=p["umap_metric"])
        sc.tl.umap(adata, min_dist=p["umap_min_dist"], spread=p["umap_spread"], init_pos=p["umap_init_pos"], random_state=p["seed"])
        sc.tl.leiden(adata, resolution=p["cluster_resolution"], key_added="cluster", random_state=p["seed"])
        adata.uns["codespring_integration"] = integration
        # These previews are available before annotation so users can assess
        # the UMAP and clustering result in the Run Pipeline step itself.
        save_umap(adata, "cluster", figures / "04_umap_clusters_pre_annotation.png", title="UMAP clusters")
        if adata.obs["sample_id"].astype(str).nunique() > 1:
            save_umap(adata, "sample_id", figures / "04_umap_samples_pre_annotation.png", title="UMAP by input sample")
        save_umap_coordinates(adata, tables)
        write_h5ad_checkpoint(adata, cluster_checkpoint)
        mark_complete("cluster")
        if stage == "cluster":
            return

    # Stage 5: annotation, markers, and the final portable results set.
    if stage in {"annotate", "all"}:
        integration = str(adata.uns.get("codespring_integration", p["integration"]))
        annotation_name = safe_metadata_name(p["annotation_name"])
        if p["celltype_file"] and p["celltype_file"].lower() != "none":
            apply_celltype_mapping(adata, Path(p["celltype_file"]).expanduser(), annotation_name)
        elif p["marker_file"] and p["marker_file"].lower() != "none":
            apply_marker_annotation(adata, Path(p["marker_file"]).expanduser(), tables, annotation_name, p["marker_species"], p["marker_ortholog_file"])
        elif not apply_existing_annotation(adata):
            adata.obs[annotation_name] = adata.obs["cluster"].astype("category")
            adata.obs[f"annotation_source__{annotation_name}"] = "cluster ID (no annotation supplied)"
            adata.obs["annotation_source"] = f"cluster ID ({annotation_name}; no annotation supplied)"
        elif annotation_name != "cell_type" and "cell_type" in adata.obs.columns:
            adata.obs[annotation_name] = adata.obs["cell_type"].copy()
            adata.obs[f"annotation_source__{annotation_name}"] = adata.obs["annotation_source"].astype(str)
        for col, name in [("sample_id", "03_umap_sample.png"), ("cluster", "04_umap_clusters.png"), (annotation_name, f"05_umap_{annotation_name}.png")]:
            save_umap(adata, col, figures / name)
        if "condition" in adata.obs.columns and adata.obs.condition.nunique() > 1:
            save_umap(adata, "condition", figures / "06_umap_condition.png")
        marker_sets = adata.uns.get("codespring_marker_gene_sets", {})
        if marker_sets:
            print("Rendering readable marker-list dot plots and heatmaps.", flush=True)
            save_marker_annotation_panels(adata, marker_sets, figures, annotation_name)
        try:
            sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon", use_raw=True)
            markers = sc.get.rank_genes_groups_df(adata, group=None)
            markers.to_csv(tables / "cluster_markers.tsv", sep="\t", index=False)
            markers.groupby("group", observed=True).head(10).to_csv(tables / "top10_markers_per_cluster.tsv", sep="\t", index=False)
            save_cluster_marker_heatmaps(adata, markers, figures)
        except Exception as exc:
            pd.DataFrame({"warning": [str(exc)]}).to_csv(tables / "cluster_markers.tsv", sep="\t", index=False)
        write_cell_metadata_outputs(adata, tables)
        write_composition_outputs(adata, annotation_name, tables)
        adata.uns["codespring_active_annotation"] = annotation_name
        write_h5ad_checkpoint(adata, processed_object)
        # The final H5AD includes the cluster representation, annotations,
        # counts, and normalized raw values, making the previous clustered
        # checkpoint unnecessary after successful annotation.
        discard_checkpoint(cluster_checkpoint)
        (out_dir / "run_summary.txt").write_text("\n".join([
            "engine: scanpy", "normalization: lognormalize", f"integration: {integration}", f"doublet_method: {p['doublet_method']}",
            f"doublets_removed: {int(adata.uns.get('codespring_doublets_removed', 0))}", "input_processing_inventory: tables/input_processing_detected.tsv",
            f"input_samples: {samples.shape[0]}", f"cells_after_qc: {adata.n_obs}", f"clusters: {adata.obs['cluster'].nunique()}", f"active_annotation: {annotation_name}", f"annotation_source: {adata.obs['annotation_source'].iloc[0]}",
        ]) + "\n")
        mark_complete("annotate")
        (out_dir / "_COMPLETE").write_text("complete\n")

    if stage == "score":
        if not p["signature_file"] or p["signature_file"].lower() == "none":
            raise SystemExit("Choose a signature TSV with signature and gene columns.")
        annotation_name = safe_metadata_name(str(adata.uns.get("codespring_active_annotation", p["annotation_name"])))
        score_signatures(adata, Path(p["signature_file"]).expanduser(), tables, annotation_name, p["signature_species"], p["signature_ortholog_file"])
        write_cell_metadata_outputs(adata, tables)
        write_h5ad_checkpoint(adata, processed_object)
        mark_complete("score")

    if stage == "differential":
        export_differential_inputs(adata, p, tables, out_dir)


if __name__ == "__main__":
    main()
