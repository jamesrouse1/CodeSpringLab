# CodeSpringLab RNA-Seq Shiny Server Bundle

This bundle lets the RNA-seq results explorer read directly from a completed CodeSpringLab run on the server.

Expected results layout:

`<results_root>/<project_name>/data`

For example:

`~/csl_results/example_dataset/data`

Files:

- `app_server.R`: config-driven Shiny app
- `shiny_results_config.R`: points the app to the completed project results

When launched by CodeSpringApp, the project config also supplies the genome species,
reference version, and exact GTF used for the RNA-seq analysis. The explorer uses that
GTF for gene-ID annotation so current M39/v50 projects are not silently annotated with
legacy M29/v42 resources.
- `run_rnaseq_results_explorer.sh`: loads modules and runs the app
- `sbatch_rnaseq_results_explorer.sh`: batch wrapper for the app
- `bulkRNAseq_shiny_notebook_chunk.py`: notebook cell snippet to print a config for the active run
