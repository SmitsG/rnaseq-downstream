# _targets.R
library(targets)
library(here)

tar_option_set(
  packages = c(
    "stats", 
    "utils",
    "airway",
    "DESeq2",
    "SummarizedExperiment",
    "ggplot2",
    "pheatmap",
    "RColorBrewer",
    "clusterProfiler",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "enrichplot",
    "DOSE",
    "EnhancedVolcano",
    "vsn",
    "matrixStats",
    "here"
  ),
  format = "rds"
)

source(here::here("R", "pipeline_functions.R"))

list(
  tar_target(
    pipeline_run,
    run_rnaseq_pipeline(alpha = 0.05, padj_cutoff = 0.05, shrink = TRUE)
  )
)
