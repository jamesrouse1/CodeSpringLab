#!/usr/bin/env python3
"""Run the complete Scanpy workflow on two small filtered 10x matrices."""

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


repo = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
pipeline = repo / "scripts_DoNotTouch" / "singleCellRNAseq" / "scrna_pipeline_scanpy.py"
assert pipeline.exists()
work = Path(tempfile.mkdtemp(prefix="codespring-scrna-scanpy-"))


def write_10x(sample_id: str, batch_shift: int) -> Path:
    path = work / sample_id
    path.mkdir()
    rng = np.random.default_rng(42 + batch_shift)
    counts = rng.poisson(0.8, size=(300, 120))
    counts[0:20, 0:60] += rng.poisson(2, size=(20, 60))
    counts[20:40, 60:120] += rng.poisson(2, size=(20, 60))
    counts[40:70, :] += batch_shift
    mmwrite(path / "matrix.mtx", sparse.coo_matrix(counts))
    genes = [f"Gene{i}" for i in range(296)] + ["MT-Co1", "MT-Co2", "Rpl1", "Rps1"]
    pd.DataFrame({0: [f"ENSG{i}" for i in range(300)], 1: genes, 2: "Gene Expression"}).to_csv(path / "features.tsv", sep="\t", header=False, index=False)
    pd.Series([f"BC{i}-1" for i in range(120)]).to_csv(path / "barcodes.tsv", sep="\t", header=False, index=False)
    return path


try:
    paths = [write_10x("sample_A", 0), write_10x("sample_B", 1)]
    manifest = pd.DataFrame({
        "sample_id": ["sample_A", "sample_B"],
        "capture_id": ["capture_1", "capture_1"],
        "input_path": [str(path) for path in paths],
        "condition": ["control", "treated"],
        "technical_batch": ["run_1", "run_2"],
    })
    manifest_path = work / "samples.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    params = {
        "normalization": "lognormalize", "integration": "harmony", "batch_column": "technical_batch",
        "cluster_resolution": "0.4", "min_features": "10", "min_counts": "0", "max_features": "0",
        "max_percent_mt": "100", "n_pcs": "10", "min_cells_per_gene": "1", "doublet_method": "none",
        "doublet_rate": "0", "remove_doublets": "false", "seed": "1234",
        "harmony_theta": "2", "harmony_lambda": "1", "harmony_max_iter": "20",
    }
    params_path = work / "params.tsv"
    pd.DataFrame({"key": list(params), "value": list(params.values())}).to_csv(params_path, sep="\t", index=False)
    out = work / "output"
    subprocess.run([sys.executable, str(pipeline), str(manifest_path), str(out), str(params_path), "all"], check=True)
    expected = [
        "_COMPLETE", "_STAGE_PREPROCESS_COMPLETE", "_STAGE_CLUSTER_COMPLETE",
        "figures/02_preintegration_umap_sample.png", "figures/02_preintegration_umap_batch.png",
        "tables/preintegration_umap_coordinates.tsv", "tables/umap_coordinates.tsv",
        "objects/processed_scanpy.h5ad",
    ]
    assert all((out / item).exists() for item in expected)
    pre = pd.read_csv(out / "tables/preintegration_umap_coordinates.tsv", sep="\t")
    final = pd.read_csv(out / "tables/umap_coordinates.tsv", sep="\t")
    assert len(pre) == len(final) == 240
    assert {"sample_id", "condition", "technical_batch"}.issubset(pre.columns)
    doublets = pd.read_csv(out / "tables/doublet_summary_by_capture.tsv", sep="\t")
    assert len(doublets) == 1 and doublets.loc[0, "capture_id"] == "capture_1"
    print("MULTI_MATRIX_SCANPY_HARMONY_OK")
finally:
    shutil.rmtree(work, ignore_errors=True)
