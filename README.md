# Bulk RNA-seq Analysis Pipeline (DESeq2 + Enrichment)

## Overview

This repository provides a **production-quality bulk RNA-seq downstream analysis pipeline** written in **R**, built around **DESeq2**, **clusterProfiler**, and **enrichplot**.

The pipeline is designed for **PhD students and researchers** who want a workflow that is:

- biologically sound
- reproducible
- readable
- easily modifiable

It goes well beyond minimal demo scripts and is intended for **real research use**.

A full reference analysis using the **airway** dataset is included.

---

## What this pipeline does

Given RNA-seq count data, the pipeline performs:

### Differential expression
- DESeq2 model fitting
- Low-count prefiltering
- Log2 fold-change shrinkage (apeglm when available)
- Extraction of up- and down-regulated genes
- Export of full and shrunk result tables

### Quality control (QC)
- PCA (rlog / VST)
- MA plots
- Dispersion plots
- Cook’s distance outlier detection
- Sample-to-sample distance heatmaps
- Top-variable gene heatmaps
- Mean–SD plots

### Functional enrichment
- GO over-representation analysis (ORA)
- KEGG ORA
- GO GSEA
- KEGG GSEA
- Proper ranked gene list generation (ENTREZ IDs)
- Gene symbol mapping
- Network plots (emapplot, cnetplot)
- Direction-aware GSEA dotplots and ridgeplots

### Outputs
- Publication-ready PNG figures
- Structured CSV result tables
- One comprehensive PDF report containing all key plots

---

## Project structure

project/  
├── R/ – modular pipeline functions  
│   ├── 00_utils.R  
│   ├── 01_deseq2.R  
│   ├── 02_qc_plots.R  
│   ├── 03_enrichment.R  
│   ├── 04_reporting.R  
│   └── pipeline_functions.R  
├── scripts/ – runnable analysis scripts  
│   └── run_airway_pipeline.R  
├── plots/ – generated figures  
├── results/ – generated result tables  
├── _targets.R – optional targets workflow  
├── renv.lock – reproducible R environment  
└── README.md  

---

## Requirements

- R version 4.2 or higher
- Bioconductor version 3.17 or higher
- Internet connection (required for KEGG enrichment)

### Main R packages used
- DESeq2
- clusterProfiler
- enrichplot
- org.Hs.eg.db
- EnhancedVolcano
- renv (recommended)
- targets (optional)

---

## How to run the pipeline

The pipeline is designed to be run from the **project root directory**.

The main entry point is a single function:

**run_rnaseq_pipeline**

By default, the pipeline is run as a linear script, which is the most readable and beginner-friendly option:
source("scripts/run_airway_pipeline.R")

Running this function executes the complete workflow, including differential expression, QC, enrichment analyses, and reporting. All outputs are written to the `plots/` and `results/` directories.

All other functions in the `R/` directory are **internal helpers** and are not intended to be called directly by users.

**Optional: Using the targets pipeline**

For advanced users or repeated analyses, an optional targets workflow is provided.

The targets package enables:

* automatic caching of intermediate results
* skipping recomputation when code or inputs have not changed
* safer reruns during development or iterative analysis

install.packages("targets")
targets::tar_make()

---

## Outputs

After a successful run, the following outputs are generated.

### Results tables (`results/`)
- deseq2_results.csv
- deseq2_results_shrunk.csv
- Upregulated_genes.csv
- Downregulated_genes.csv
- GO_ORA_upregulated.csv
- KEGG_ORA.csv
- GO_GSEA.csv
- KEGG_GSEA.csv
- Ranked_Genes_for_GSEA_with_SYMBOL.csv
- pipeline.log

### Figures (`plots/`)
- PCA plots
- MA plots
- Volcano plot
- Sample distance heatmap
- Top-variable gene heatmap
- Enrichment dotplots
- GSEA ridgeplots
- Network plots
- Airway_DE_Analysis_All_Plots.pdf

---

## Using your own data

`run_rnaseq_pipeline()` uses the example `airway` dataset when called without data arguments. The original `run_airway_pipeline()` name remains available for compatibility. To analyze your own data, pass a raw count matrix and sample metadata table:

```r
source("R/pipeline_functions.R")

counts <- as.matrix(read.csv("data/counts.csv", row.names = 1, check.names = FALSE))
sample_data <- read.csv("data/samples.csv", row.names = 1)
sample_data$condition <- relevel(factor(sample_data$condition), ref = "control")

out <- run_rnaseq_pipeline(
  counts = counts,
  col_data = sample_data,
  design = ~ condition,
  group_variable = "condition"
)
```

The count matrix must have genes in rows, samples in columns, unique gene IDs as row names, and unique sample IDs as column names. Counts must be raw, finite, non-negative whole numbers; normalized or transformed values are not valid DESeq2 input. The sample metadata must have sample IDs as row names, with names matching the count matrix columns. Set the condition reference level explicitly so the reported fold-change direction is clear. For designs with more than one coefficient, pass `coef` to select the coefficient used for LFC shrinkage.

The current enrichment functions use human Ensembl IDs (`org.Hs.eg.db`) and human KEGG (`hsa`). Use human Ensembl gene IDs for custom data, with version suffixes removed from the count matrix row names. Non-human organisms and other identifier types need a matching annotation database and enrichment configuration before those results are meaningful. The combined PDF report is named with an `RNA-seq_` prefix for custom data and an `Airway_` prefix for the example; other output filenames are shared.

---

## targets workflow (optional)

An optional workflow using the **targets** package is provided.

This enables caching and avoids recomputation when inputs or code have not changed.

The script-based workflow is the **canonical and most readable version**.  
The targets workflow is provided as an execution wrapper, not a replacement.

---

## Reproducibility

This project uses **renv** to manage package versions.

The `renv.lock` file records the exact dependency versions used for the analysis, enabling full reproducibility across systems.

---

## Citation

If you use this pipeline in a publication or thesis, please cite the repository.
A `CITATION.cff` file is provided for convenience and can be accessed via the
" Cite this repository " button on GitHub.


---

## Philosophy

This pipeline prioritizes:
- clarity over cleverness
- biological correctness over minimal code
- reproducibility over speed
- extensibility without premature over-engineering

It is intended as a **real analysis template**, not a toy example.
