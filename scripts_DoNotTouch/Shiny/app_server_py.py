#!/usr/bin/env python3

import os
import re
from pathlib import Path

import pandas as pd
from shiny import App, render, ui


def parse_r_config(path: Path):
    cfg = {}
    if not path.exists():
        return cfg

    assign_re = re.compile(r"^\s*([A-Za-z0-9_]+)\s*<-\s*(.+?)\s*$")
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        m = assign_re.match(line)
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        if val.startswith('"') and val.endswith('"'):
            cfg[key] = val[1:-1]
        else:
            try:
                cfg[key] = int(val)
            except ValueError:
                cfg[key] = val
    return cfg


def safe_read_tsv(path: Path):
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def list_matching_files(path: Path, pattern: str):
    if not path.exists():
        return []
    reg = re.compile(pattern)
    return sorted([p.name for p in path.rglob("*") if p.is_file() and reg.search(p.name)])


default_config = Path(os.environ.get("RNASEQ_SHINY_CONFIG", "./shiny_results_config.R")).expanduser()
cfg = parse_r_config(default_config)

project_name = str(cfg.get("project_name", os.environ.get("RNASEQ_SHINY_PROJECT", "example_dataset")))
results_root = Path(
    str(cfg.get("results_root", os.environ.get("RNASEQ_SHINY_RESULTS_ROOT", "~/csl_results")))
).expanduser()
data_dir = Path(str(cfg.get("data_dir", results_root / project_name / "data"))).expanduser()

design_matrix_path = Path(
    str(cfg.get("design_matrix_path", data_dir / "design_matrix" / "design_matrix.txt"))
).expanduser()
star_summary_path = data_dir / "star_summary" / "summary_matrix.txt"
count_matrix_path = data_dir / "counts" / "count_matrix.txt"
featurecounts_summary_path = data_dir / "counts" / "featurecounts_summary.txt"
deseq2_dir = data_dir / "deseq2"
gseapy_dir = data_dir / "gseapy"

design_df = safe_read_tsv(design_matrix_path)
star_df = safe_read_tsv(star_summary_path)
count_df = safe_read_tsv(count_matrix_path)
featurecounts_df = safe_read_tsv(featurecounts_summary_path)

deseq2_files = list_matching_files(deseq2_dir, r"\.txt$|\.csv$|\.pdf$|\.png$")
gseapy_files = list_matching_files(gseapy_dir, r"\.csv$|\.pdf$|\.png$")

app_ui = ui.page_fluid(
    ui.h2("RNA-Seq Results Explorer (Py-Shiny)"),
    ui.markdown(
        f"""
**Project:** `{project_name}`  
**Data dir:** `{data_dir}`  
**Config:** `{default_config}`
"""
    ),
    ui.hr(),
    ui.navset_tab(
        ui.nav(
            "Overview",
            ui.h4("Detected resources"),
            ui.output_text_verbatim("overview_status"),
        ),
        ui.nav(
            "Design Matrix",
            ui.output_text_verbatim("design_tbl"),
        ),
        ui.nav(
            "STAR Summary",
            ui.output_text_verbatim("star_tbl"),
        ),
        ui.nav(
            "Count Matrix",
            ui.p("Showing first 500 rows for performance."),
            ui.output_text_verbatim("count_tbl"),
        ),
        ui.nav(
            "FeatureCounts Summary",
            ui.output_text_verbatim("featurecounts_tbl"),
        ),
        ui.nav(
            "DESeq2 Files",
            ui.output_text_verbatim("deseq2_files_txt"),
        ),
        ui.nav(
            "GSEA Files",
            ui.output_text_verbatim("gsea_files_txt"),
        ),
    ),
)


def server(input, output, session):
    def fmt_df(df, max_rows=500):
        if df is None or df.empty:
            return "<empty>"
        return df.head(max_rows).to_string(index=False)

    @output
    @render.text
    def overview_status():
        lines = [
            f"design_matrix: {'FOUND' if not design_df.empty else 'MISSING'} ({design_matrix_path})",
            f"star_summary: {'FOUND' if not star_df.empty else 'MISSING'} ({star_summary_path})",
            f"count_matrix: {'FOUND' if not count_df.empty else 'MISSING'} ({count_matrix_path})",
            f"featurecounts_summary: {'FOUND' if not featurecounts_df.empty else 'MISSING'} ({featurecounts_summary_path})",
            f"deseq2 files: {len(deseq2_files)}",
            f"gseapy files: {len(gseapy_files)}",
        ]
        return "\n".join(lines)

    @output
    @render.text
    def design_tbl():
        return fmt_df(design_df)

    @output
    @render.text
    def star_tbl():
        return fmt_df(star_df)

    @output
    @render.text
    def count_tbl():
        show_df = count_df.head(500) if not count_df.empty else count_df
        return fmt_df(show_df)

    @output
    @render.text
    def featurecounts_tbl():
        return fmt_df(featurecounts_df)

    @output
    @render.text
    def deseq2_files_txt():
        if not deseq2_files:
            return f"No DESeq2 output files found under: {deseq2_dir}"
        return "\n".join(deseq2_files)

    @output
    @render.text
    def gsea_files_txt():
        if not gseapy_files:
            return f"No GSEA output files found under: {gseapy_dir}"
        return "\n".join(gseapy_files)


app = App(app_ui, server)


if __name__ == "__main__":
    host = os.environ.get("RNASEQ_SHINY_HOST_OVERRIDE", str(cfg.get("host", "0.0.0.0")))
    port = int(os.environ.get("RNASEQ_SHINY_PORT_OVERRIDE", str(cfg.get("port", 3838))))
    app.run(host=host, port=port)
