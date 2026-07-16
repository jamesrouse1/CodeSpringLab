# Bundled SEACR

CodeSpringLab bundles both `SEACR_1.3.sh` and its required companion `SEACR_1.3.R` for CUT&RUN peak calling so the workflow can run without an extra manual download step.

Peak calls are stored by analysis combination under `data/seacr/<normalization>_<stringency>/<sample>/`, for example `data/seacr/non_stringent/AKP_Creb-AA1/`. Running another combination therefore preserves earlier results. `norm` runs use raw fragment bedGraphs so SEACR performs target/control normalization internally; `non` runs use the already normalized Bowtie2 bedGraphs.

Source:

```text
https://github.com/FredHutch/SEACR
```

Citation requested by SEACR:

```text
Meers MP, Tenenbaum D, Henikoff S. Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics & Chromatin. 2019;12:42.
```

License:

```text
GPL-2.0
```

The bundled license is stored as `LICENSE_SEACR`.

To refresh the bundled copy:

```bash
bash ../scripts_DoNotTouch/SEACR/download_seacr.sh
```
