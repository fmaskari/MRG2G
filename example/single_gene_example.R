source("MRG2G.R")
source("simulate_g2g_data.R")

d <- simulate_mrg2g_single_gene(n = 300, signal = 0.25, cov_structure = "AR1",
                                standardize = TRUE, seed = 1)

gs <- make_gene_structure(ld_blocks = d$design$ld_blocks,
                          cg_blocks = d$design$cg_blocks)

dat <- mrg2g_setup(Z = d$Z, X = d$X, y = d$y, U = NULL,
                   gene_structure = gs,
                   region_names   = d$design$block_names,
                   standardize    = TRUE)

fit <- mrg2g_fit(dat,
                 lambda       = 0.18,
                 a1           = 0.3,
                 a2           = 0.2,
                 stage2       = c("ols", "lasso", "en", "sgl"),
                 lambda_scale = 0.5)

cat("\nStage 1 selection\n")
cat(sprintf("  IVs selected:  %d (true %d)\n",
            sum(fit$IVs_selected), sum(d$truth$IVs_true)))
cat(sprintf("  CpGs selected: %d (true %d)\n",
            sum(fit$methyls_selected), sum(d$truth$CpG_select_true)))
cat(sprintf("  Regions kept:  %s\n", paste(fit$region_names_sel, collapse = ", ")))

cat("\nStage 2 MSE\n")
for (m in names(fit$estimates)) {
  cat(sprintf("  %-6s %.5f\n", m,
              mean((fit$estimates[[m]]$theta - d$truth$theta_std)^2)))
}

if (sum(fit$methyls_selected) >= 2) {
  wt <- mrg2g_wald(fit, corrections = "all")
  wt$null_block <- d$design$null_blocks[match(wt$block_name, d$design$block_names)]
  cat("\nRegional Wald test\n")
  print(wt[, c("block_name", "null_block", "wald_stat", "df", "pval_raw",
               "padj_bh", "reject_bh")])
}

