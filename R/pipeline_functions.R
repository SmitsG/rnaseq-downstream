# R/pipeline_functions.R
suppressPackageStartupMessages({
  library(here)
})

source(here::here("R", "00_utils.R"))
source(here::here("R", "01_deseq2.R"))
source(here::here("R", "02_qc_plots.R"))
source(here::here("R", "03_enrichment.R"))
source(here::here("R", "04_reporting.R"))

suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
  library(SummarizedExperiment)
  library(ggplot2)
  library(EnhancedVolcano)
})

run_rnaseq_pipeline <- function(alpha = 0.05,
                                padj_cutoff = 0.05,
                                shrink = TRUE,
                                coef = NULL,
                                counts = NULL,
                                col_data = NULL,
                                design = ~ dex,
                                group_variable = NULL) {
  dirs <- init_output_dirs()
  log_path <- file.path(dirs$results, "pipeline.log")
  safe_remove(log_path)
  
  if (is.null(counts) && is.null(col_data)) {
    data("airway", package = "airway")
    se <- airway
    design <- ~ dex
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    counts <- SummarizedExperiment::assay(se)
    group_variable <- "dex"
    dataset_name <- "Airway"
  } else {
    dataset_name <- "RNA-seq"
    assert_that(!is.null(counts) && !is.null(col_data),
                "Provide both counts and col_data, or neither to use the airway example.")
    assert_that(is.matrix(counts) && is.numeric(counts),
                "counts must be a numeric matrix with genes in rows and samples in columns.")
    assert_that(nrow(counts) > 0 && ncol(counts) > 0,
                "counts must contain at least one gene and one sample.")
    assert_that(!is.null(rownames(counts)) && !anyDuplicated(rownames(counts)),
                "counts must have unique gene IDs as row names.")
    assert_that(!is.null(colnames(counts)) && !anyDuplicated(colnames(counts)),
                "counts must have unique sample IDs as column names.")
    assert_that(all(is.finite(counts)) && all(counts >= 0) &&
                  all(counts == floor(counts)),
                "counts must contain finite, non-negative whole numbers (raw counts).")
    assert_that(is.data.frame(col_data) && !is.null(rownames(col_data)),
                "col_data must be a data.frame with sample IDs as row names.")
    assert_that(!anyDuplicated(rownames(col_data)),
                "col_data must have unique sample IDs as row names.")
    assert_that(setequal(colnames(counts), rownames(col_data)),
                "Sample IDs in counts columns must match row names in col_data.")
    col_data <- col_data[colnames(counts), , drop = FALSE]
    design_variables <- all.vars(design)
    assert_that(length(design_variables) > 0 &&
                  all(design_variables %in% colnames(col_data)),
                "Every variable in design must be a column in col_data.")
    assert_that(all(vapply(col_data[design_variables], function(x) !anyNA(x), logical(1))),
                "Design variables in col_data cannot contain missing values.")
    if (is.null(group_variable)) group_variable <- tail(design_variables, 1)
    assert_that(group_variable %in% colnames(col_data),
                "group_variable must name a column in col_data.")
    col_data[[group_variable]] <- droplevels(factor(col_data[[group_variable]]))
    assert_that(nlevels(col_data[[group_variable]]) >= 2,
                "group_variable must have at least two observed groups.")
  }

  log_line(log_path, "Starting", dataset_name, "pipeline")
  log_line(log_path, "Loaded dataset | samples:", ncol(counts), "genes:", nrow(counts))

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = col_data,
    design = design
  )
  dds <- prefilter_dds(dds, min_total = 1)
  dds <- run_deseq2(dds)
  log_line(log_path, "DESeq2 fit complete | genes:", nrow(dds), "samples:", ncol(dds))
  
  res <- get_deseq_results(dds, alpha = alpha)
  save_csv(as.data.frame(res), file.path(dirs$results, "deseq2_results.csv"))
  
  if (isTRUE(shrink) && is.null(coef)) {
    available_coefs <- setdiff(DESeq2::resultsNames(dds), "Intercept")
    assert_that(length(available_coefs) == 1,
                paste0("Set coef explicitly for LFC shrinkage. Available coefficients: ",
                       paste(available_coefs, collapse = ", ")))
    coef <- available_coefs[[1]]
  }
  res_shrunk <- if (isTRUE(shrink)) shrink_lfc_safe(dds, coef = coef) else res
  res_shrunk <- res_shrunk[order(res_shrunk$pvalue), ]
  save_csv(as.data.frame(res_shrunk), file.path(dirs$results, "deseq2_results_shrunk.csv"))
  
  tx <- make_transforms(dds)
  vsd <- tx$vsd
  rld <- tx$rld
  ann <- as.data.frame(SummarizedExperiment::colData(dds)[, group_variable, drop = FALSE])
  
  pca <- plot_pca(rld, intgroup = group_variable,
                  title = paste("PCA (rlog) -", dataset_name))
  save_gg(pca, file.path(dirs$plots, "PCA_plot.png"))
  
  plot_ma_png(res, file.path(dirs$plots, "MA_plot.png"))
  plot_dispersion_png(dds, file.path(dirs$plots, "Dispersion_plot.png"))
  plot_cooks_png(dds, file.path(dirs$plots, "Cooks_distance_boxplot.png"))
  plot_sample_distance_heatmap(vsd, ann, file.path(dirs$plots, "Sample_distance_heatmap.png"))
  plot_top_variable_heatmap(vsd, ann, n = 20, path = file.path(dirs$plots, "Top20_VariableGenes_Heatmap.png"))
  
  ms <- plot_mean_sd(vsd, rld)
  save_gg(ms$vsd, file.path(dirs$plots, "MeanSD_VST.png"))
  save_gg(ms$rld, file.path(dirs$plots, "MeanSD_RLD.png"))
  
  volcano <- EnhancedVolcano::EnhancedVolcano(
    res_shrunk,
    lab = rownames(res_shrunk),
    x = "log2FoldChange",
    y = "pvalue",
    pCutoff = 0.05,
    FCcutoff = 1.5,
    title = paste("Volcano plot -", dataset_name, "(shrunk LFC if available)")
  )
  save_gg(volcano, file.path(dirs$plots, "Volcano_plot.png"), width = 7, height = 5)
  
  sig <- get_sig_genes(res, padj_cutoff = padj_cutoff)
  save_csv(data.frame(ENSEMBL = sig$up), file.path(dirs$results, "Upregulated_genes.csv"))
  save_csv(data.frame(ENSEMBL = sig$down), file.path(dirs$results, "Downregulated_genes.csv"))
  log_line(log_path, "Sig genes | up:", length(sig$up), "down:", length(sig$down))
  
  ego <- run_go_ora(sig$up, ont = "BP")
  save_csv(as.data.frame(ego), file.path(dirs$results, "GO_ORA_upregulated.csv"))
  
  go_dot <- safe_dotplot(ego, "GO BP ORA - Upregulated", showCategory = 10)
  if (!is.null(go_dot)) save_gg(go_dot, file.path(dirs$plots, "GO_Enrichment_Dotplot.png"), width = 8, height = 6)
  
  go_emap <- go_emapplot_safe(ego, showCategory = 10)
  if (!is.null(go_emap)) save_gg(go_emap, file.path(dirs$plots, "GO_emapplot.png"), width = 10, height = 8)
  
  go_cnet <- go_cnetplot_safe(ego, showCategory = 5)
  if (!is.null(go_cnet)) save_gg(go_cnet, file.path(dirs$plots, "GO_cnetplot.png"), width = 10, height = 8)
  
  rank_obj <- build_entrez_rank(res_shrunk)
  lfc_rank <- rank_obj$lfc_rank
  symbols <- add_symbols(names(lfc_rank))
  
  ranked_df <- data.frame(
    EntrezID = names(lfc_rank),
    Symbol = unname(symbols),
    log2FC = as.numeric(lfc_rank)
  )
  save_csv(ranked_df, file.path(dirs$results, "Ranked_Genes_for_GSEA_with_SYMBOL.csv"))
  log_line(log_path, "Ranked Entrez genes:", nrow(ranked_df))
  
  kegg_ora <- run_kegg_ora(names(lfc_rank), organism = "hsa")
  save_csv(as.data.frame(kegg_ora), file.path(dirs$results, "KEGG_ORA.csv"))
  
  kegg_dot <- safe_dotplot(kegg_ora, "KEGG ORA", showCategory = 10)
  if (!is.null(kegg_dot)) save_gg(kegg_dot, file.path(dirs$plots, "KEGG_Enrichment_Dotplot.png"), width = 8, height = 6)
  
  gse_go <- run_go_gsea(lfc_rank, ont = "BP")
  save_csv(as.data.frame(gse_go), file.path(dirs$results, "GO_GSEA.csv"))
  
  gse_kegg <- run_kegg_gsea(lfc_rank, organism = "hsa")
  save_csv(as.data.frame(gse_kegg), file.path(dirs$results, "KEGG_GSEA.csv"))
  
  ridge_go <- safe_ridgeplot(gse_go, "GSEA GO ridgeplot", showCategory = 15)
  ridge_kegg <- safe_ridgeplot(gse_kegg, "GSEA KEGG ridgeplot", showCategory = 15)
  if (!is.null(ridge_go)) save_gg(ridge_go, file.path(dirs$plots, "GSEA_GO_ridgeplot.png"), width = 9, height = 6)
  if (!is.null(ridge_kegg)) save_gg(ridge_kegg, file.path(dirs$plots, "GSEA_KEGG_ridgeplot.png"), width = 9, height = 6)
  
  dot_go_sign <- safe_dotplot_split(gse_go, "GSEA GO dotplot (split sign)", showCategory = 10)
  dot_kegg_sign <- safe_dotplot_split(gse_kegg, "GSEA KEGG dotplot (split sign)", showCategory = 10)
  if (!is.null(dot_go_sign)) save_gg(dot_go_sign, file.path(dirs$plots, "GSEA_GO_dotplot_sign.png"), width = 10, height = 6)
  if (!is.null(dot_kegg_sign)) save_gg(dot_kegg_sign, file.path(dirs$plots, "GSEA_KEGG_dotplot_sign.png"), width = 10, height = 6)
  
  pdf_path <- file.path(dirs$plots, paste0(dataset_name, "_DE_Analysis_All_Plots.pdf"))
  write_pdf_report(
    pdf_path,
    ggplots = list(
      pca, volcano, go_dot, kegg_dot,
      ms$vsd, ms$rld,
      ridge_go, ridge_kegg,
      dot_go_sign, dot_kegg_sign,
      go_emap, go_cnet
    ),
    base_plot_fns = list(function() DESeq2::plotMA(res, main = paste("MA plot -", dataset_name)))
  )
  
  log_line(log_path, "Saved PDF report:", pdf_path)
  log_line(log_path, "Finished airway pipeline")
  
  list(
    dds = dds,
    res = res,
    res_shrunk = res_shrunk,
    vsd = vsd,
    rld = rld,
    ego = ego,
    kegg_ora = kegg_ora,
    gse_go = gse_go,
    gse_kegg = gse_kegg,
    ranked = ranked_df,
    dirs = dirs,
    log = log_path
  )
}

# Retain the original function name for existing airway examples and scripts.
run_airway_pipeline <- run_rnaseq_pipeline
