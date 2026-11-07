# scripts/run_airway_pipeline.R
suppressPackageStartupMessages({
  library(here)
})

source(here::here("R", "pipeline_functions.R"))

out <- run_rnaseq_pipeline(alpha = 0.05, padj_cutoff = 0.05, shrink = TRUE)

message("Done.")
message("Results: ", out$dirs$results)
message("Plots:   ", out$dirs$plots)
message("Log:     ", out$log)
