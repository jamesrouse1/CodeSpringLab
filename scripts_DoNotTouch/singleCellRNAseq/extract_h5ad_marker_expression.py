#!/usr/bin/env python3
"""Extract one normalized gene-expression vector from a Scanpy H5AD object."""
import argparse
from pathlib import Path

import anndata
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad")
    parser.add_argument("gene")
    parser.add_argument("output")
    args = parser.parse_args()
    adata = anndata.read_h5ad(args.h5ad, backed="r")
    if adata.raw is None:
        raise SystemExit("This H5AD has no normalized raw expression layer.")
    genes = pd.Index(adata.raw.var_names.astype(str))
    index = genes.get_indexer([args.gene])[0]
    if index < 0:
        raise SystemExit(f"Gene not found in normalized raw expression: {args.gene}")
    values = adata.raw.X[:, index]
    if hasattr(values, "toarray"):
        values = values.toarray()
    values = np.asarray(values).reshape(-1)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    pd.DataFrame({"cell": adata.obs_names.astype(str), "expression": values}).to_csv(temporary, sep="\t", index=False, compression="gzip")
    temporary.replace(output)


if __name__ == "__main__":
    main()
