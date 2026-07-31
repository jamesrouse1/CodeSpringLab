#!/usr/bin/env python3
"""CodeSpringLab single-cell RNA-seq workflow (Scanpy engine).

Usage: python scrna_pipeline_scanpy.py <samples.tsv> <out_dir> <params.tsv>

The manifest and parameter files are shared with the Seurat engine.  This
script intentionally preserves raw counts in ``layers['counts']`` and uses a
scVI latent representation for multi-sample integration when requested.
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
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color=color, ax=ax, show=False, title="UMAP supplied with input object", frameon=False)
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
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sc.pl.umap(adata, color=color, ax=ax, show=False, title=title or str(color), frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def save_qc_plots(adata, figures: Path):
    """Write small, portable QC plots without an interactive display."""
    import matplotlib.pyplot as plt
    sc.pl.violin(
        adata,
        keys=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        groupby="sample_id",
        multi_panel=True,
        rotation=35,
        show=False,
    )
    plt.gcf().set_size_inches(11, 5)
    plt.gcf().tight_layout()
    plt.gcf().savefig(figures / "01_qc_violin.png", dpi=160)
    plt.close(plt.gcf())
    fig, ax = plt.subplots(figsize=(7, 5.5))
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", color="sample_id", ax=ax, show=False)
    fig.tight_layout()
    fig.savefig(figures / "02_qc_counts_vs_mt.png", dpi=160)
    plt.close(fig)


def save_doublet_plot(adata, figures: Path):
    """Show the score distribution used for a transparent doublet call."""
    import matplotlib.pyplot as plt
    if "doublet_score" not in adata.obs:
        return
    score = pd.to_numeric(adata.obs["doublet_score"], errors="coerce").dropna()
    if score.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(score, bins=min(50, max(10, int(np.sqrt(len(score))))), color="#5b7db1", edgecolor="white")
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
        report.append({"sample_id": "all", "method": "none", "cells_before": adata.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": "Doublet removal disabled"})
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
        # Scrublet is not informative on very small partitions. Do not invent
        # calls; retain these cells and record why in the summary.
        if sub.n_obs < 100 or sub.n_vars < 20:
            report.append({"sample_id": sample_id, "method": "scrublet", "cells_before": sub.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": "Skipped: fewer than 100 cells or 20 genes"})
            continue
        try:
            # The internal Scrublet feature filter can leave fewer than 30
            # dimensions, so choose a conservative PC count for small runs.
            scrublet_pcs = min(20, max(2, min(sub.n_obs, sub.n_vars) - 2))
            sc.pp.scrublet(sub, expected_doublet_rate=p["doublet_rate"], n_prin_comps=scrublet_pcs, random_state=p["seed"], verbose=False)
            adata.obs.loc[sub.obs_names, "doublet_score"] = pd.to_numeric(sub.obs["doublet_score"], errors="coerce").values
            predicted = sub.obs["predicted_doublet"].astype(bool)
            adata.obs.loc[sub.obs_names, "predicted_doublet"] = predicted.values
            predicted_n = int(predicted.sum())
            report.append({"sample_id": sample_id, "method": "scrublet", "cells_before": sub.n_obs, "predicted_doublets": predicted_n, "removed_doublets": predicted_n if p["remove_doublets"] else 0, "note": ""})
        except Exception as exc:
            if method == "scrublet":
                raise SystemExit(f"Scrublet failed for {sample_id}: {exc}") from exc
            report.append({"sample_id": sample_id, "method": "scrublet", "cells_before": sub.n_obs, "predicted_doublets": 0, "removed_doublets": 0, "note": f"Automatic Scrublet skipped after error: {exc}"})
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
    if len(sys.argv) != 4:
        raise SystemExit("Usage: scrna_pipeline_scanpy.py <samples.tsv> <out_dir> <params.tsv>")
    samples_path, out_dir, params_path = map(Path, sys.argv[1:])
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = out_dir / "tables"; figures = out_dir / "figures"; objects = out_dir / "objects"
    for d in (tables, figures, objects): d.mkdir(parents=True, exist_ok=True)
    p = params_from(params_path)
    if p["normalization"] not in {"lognormalize", "log1p"}:
        raise SystemExit("The Scanpy engine supports LogNormalize/log1p normalization. Use the Seurat engine for SCTransform.")
    if p["integration"] not in {"auto", "none", "scvi", "harmony"}:
        raise SystemExit("Scanpy integration must be auto, none, scvi, or harmony.")
    np.random.seed(p["seed"])
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
    obs = [item[0] for item in inputs]
    input_status = pd.DataFrame([item[1] for item in inputs])
    input_status.to_csv(tables / "input_processing_detected.tsv", sep="\t", index=False)
    adata = ad.concat(obs, join="outer", merge="same", index_unique=None, fill_value=0)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    cells_before_qc = adata.obs.groupby("sample_id", observed=True).size().rename("cells_input")
    keep = (adata.obs["n_genes_by_counts"] >= p["min_features"]) & (adata.obs["total_counts"] >= p["min_counts"]) & (adata.obs["pct_counts_mt"] <= p["max_percent_mt"])
    if p["max_features"] > 0:
        keep &= adata.obs["n_genes_by_counts"] <= p["max_features"]
    adata = adata[keep].copy()
    if adata.n_obs == 0:
        raise SystemExit("QC thresholds removed every cell.")
    # Run doublet detection on raw counts after initial QC and before every
    # expression transformation or integration step.
    adata = run_scrublet(adata, p, tables, figures)
    sc.pp.filter_genes(adata, min_cells=p["min_cells_per_gene"])
    adata.layers["counts"] = adata.X.copy()
    qc = adata.obs[["sample_id", "n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_score", "predicted_doublet"]].copy()
    qc.insert(0, "cell", adata.obs_names)
    qc.to_csv(tables / "qc_cell_metrics.tsv", sep="\t")
    qc_summary = qc.groupby("sample_id", observed=True).agg(cells_after_qc_and_doublets=("total_counts", "size"), median_umis=("total_counts", "median"), median_genes=("n_genes_by_counts", "median"), median_percent_mt=("pct_counts_mt", "median")).reset_index()
    qc_summary.insert(1, "cells_input", qc_summary["sample_id"].map(cells_before_qc).astype(int))
    qc_summary.to_csv(tables / "qc_summary_by_sample.tsv", sep="\t", index=False)
    save_qc_plots(adata, figures)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=min(3000, adata.n_vars), flavor="seurat_v3", layer="counts")
    except Exception:
        sc.pp.highly_variable_genes(adata, n_top_genes=min(3000, adata.n_vars), flavor="cell_ranger")
    hvg_table = adata.var.copy()
    hvg_table.insert(0, "gene", hvg_table.index.astype(str))
    hvg_columns = [c for c in ["gene", "highly_variable", "means", "variances", "variances_norm"] if c in hvg_table.columns]
    hvg_table.loc[:, hvg_columns].to_csv(tables / "highly_variable_genes.tsv", sep="\t", index=False)
    # Do not densify a real sparse count matrix: preserving sparsity is vital
    # for practical single-cell jobs with tens of thousands of genes.
    sc.pp.scale(adata, zero_center=False, max_value=10)
    hvg_count = int(adata.var["highly_variable"].sum()) if "highly_variable" in adata.var else 0
    use_hvg = hvg_count >= 2
    pca_features = hvg_count if use_hvg else adata.n_vars
    max_pcs = max(1, min(adata.n_obs - 1, pca_features - 1))
    n_pcs = min(p["n_pcs"], max_pcs)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", use_highly_variable=use_hvg, zero_center=False)
    pca_ratio = np.asarray(adata.uns["pca"]["variance_ratio"]).ravel()
    pd.DataFrame({"PC": np.arange(1, len(pca_ratio) + 1), "variance_explained": pca_ratio, "percent_variance_explained": 100 * pca_ratio}).to_csv(tables / "pca_variance_explained.tsv", sep="\t", index=False)
    integration = p["integration"]
    batch_key = p["batch_column"] if p["batch_column"] in adata.obs.columns else "sample_id"
    batch_count = int(adata.obs[batch_key].astype(str).nunique())
    if integration == "auto":
        integration = "scvi" if p["batch_column"] in adata.obs.columns and batch_count > 1 else "none"
    if integration in {"scvi", "harmony"} and batch_count < 2:
        raise SystemExit(f"{integration} integration requires at least two values in the selected batch column ({batch_key}). Choose none or supply the appropriate technical batch column.")
    representation = "X_pca"
    # A single AnnData object can still contain multiple technical batches.
    # Integrate based on the selected metadata values, not merely the number
    # of manifest rows, so imported atlas objects behave like multi-folder
    # inputs.
    if integration == "scvi" and batch_count > 1:
        try:
            import scvi
        except ImportError as exc:
            raise SystemExit("scVI integration was requested but scvi-tools is unavailable. Install scvi-tools or choose Harmony/none.") from exc
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
        model = scvi.model.SCVI(adata, n_latent=min(30, p["n_pcs"]))
        model.train(max_epochs=p["scvi_max_epochs"], early_stopping=True)
        adata.obsm["X_scVI"] = model.get_latent_representation()
        representation = "X_scVI"
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
    adata.obs.to_csv(tables / "cell_metadata.tsv", sep="\t")
    adata.obs.groupby(["cluster", "cell_type"], observed=True).size().reset_index(name="cells").to_csv(tables / "cluster_cell_type_sizes.tsv", sep="\t", index=False)
    adata.write_h5ad(objects / "processed_scanpy.h5ad", compression="gzip")
    (out_dir / "run_summary.txt").write_text("\n".join([
        "engine: scanpy", "normalization: lognormalize", f"integration: {integration}",
        f"doublet_method: {p['doublet_method']}", f"doublets_removed: {int(adata.uns.get('codespring_doublets_removed', 0))}",
        "input_processing_inventory: tables/input_processing_detected.tsv",
        f"input_samples: {samples.shape[0]}", f"cells_after_qc: {adata.n_obs}",
        f"clusters: {adata.obs['cluster'].nunique()}", f"annotation_source: {adata.obs['annotation_source'].iloc[0]}",
    ]) + "\n")
    (out_dir / "_COMPLETE").write_text("complete\n")


if __name__ == "__main__":
    main()
