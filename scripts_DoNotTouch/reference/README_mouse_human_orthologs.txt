Mouse-human ortholog table bundled for CodeSpringLab GSEA/pathway analysis.

Source:
http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt

Bundled clean table:
scripts_DoNotTouch/reference/mouse_human_orthologs_MGI.tsv

Rows:
24,585 mouse-human ortholog pairs from MGI/JAX homology classes.

Columns:
mgi_homology_class: MGI homology class key
mouse_gene_symbol: official mouse gene symbol from MGI row
mouse_entrez_id: mouse EntrezGene ID
mouse_mgi_id: MGI accession
human_gene_symbol: official human gene symbol from matching human row
human_entrez_id: human EntrezGene ID
human_hgnc_id: HGNC accession
human_omim_id: OMIM gene ID when present

Pipeline behavior:
bulkRNAseq.py uses this bundled table for mouse-to-human pathway gene mapping.
For GSEA, it keeps one-to-one and many-mouse-to-one-human mappings, while excluding mouse-to-many-human and many-to-many mappings.
When multiple retained mouse genes map to one human gene, their normalized expression values are averaged.
