#!/usr/bin/env python3
"""Write the complete normalized gene list from a processed Scanpy object."""
import argparse
from pathlib import Path

import anndata
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad")
    parser.add_argument("output")
    args = parser.parse_args()

    adata = anndata.read_h5ad(args.h5ad, backed="r")
    features = adata.raw.var_names if adata.raw is not None else adata.var_names
    genes = pd.Index(features.astype(str)).dropna().unique().sort_values()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": genes}).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()
