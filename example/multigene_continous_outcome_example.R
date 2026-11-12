source("MRG2G.R")
source("simulate_g2g_data.R")

# Matches MG_ar1_mod_withU.R: 3 genes (B1, B2, B3 designs), 3 confounders,
# AR1 methylation correlation (0.85 / 0.3 / 0.05, rho_gene = 0), signal 0.25
d <- simulate_mrg2g_example(n = 300, signal = 0.25, cov_structure = "AR1",
                standardize = TRUE, seed = 1, rho_gene = 0.03,
                rho_within = 0.25, rho_decay = 0.85, rho_block = 0.05)
gs <- make_gene_structure(
  ld_blocks = lapply(d$design$gene_structure, `[[`, "ld_blocks"),
  cg_blocks = lapply(d$design$gene_structure, `[[`, "cg_blocks"))

dat <- mrg2g_setup(Z = d$Z, X = d$X, y = d$y, U = d$U,
                   gene_structure = gs,
                   region_names   = d$design$block_names,
                   standardize    = TRUE)

# Stage 1 tuning: 5-fold CV predicting X from Z (y and truth never used).
# Structure (a1, a2) by CV minimum, lambda by the 1-SE rule, per gene.
# Set cores > 1 on Mac/Linux to run the fits in parallel.
tn <- mrg2g_tune(dat)
cat("\nTuned Stage 1 parameters\n")
print(tn$best_cis, row.names = FALSE)
cat(sprintf("  lambda_trans = %.3f\n", tn$lambda_trans))

fit <- mrg2g_fit(dat,
                 lambda       = tn$lambda,
                 a1           = tn$a1,
                 a2           = tn$a2,
                 lambda_trans = tn$lambda_trans,
                 stage2       = c("ols", "lasso", "en", "sgl"))

cat("\nStage 1 selection\n")
cat(sprintf("  IVs selected:  %d (true %d)\n",
            sum(fit$IVs_selected), sum(d$truth$IVs_true)))
cat(sprintf("  CpGs selected: %d (true %d)\n",
            sum(fit$methyls_selected), sum(d$truth$CpG_select_true)))
cat(sprintf("  Regions kept:  %s\n", paste(fit$region_names_sel, collapse = ", ")))
cat(sprintf("  SGL lambda_scale: %.1f (combined signal %.3f)\n",
            fit$lambda_scale_used, fit$combined_signal))

cat("\nStage 2 MSE vs standardized truth\n")
for (m in names(fit$estimates)) {
  cat(sprintf("  %-6s theta %.5f | gamma %.5f\n", m,
              mean((fit$estimates[[m]]$theta - d$truth$theta_std)^2),
              mean((fit$estimates[[m]]$gamma - d$truth$gamma_std)^2)))
}

if (sum(fit$methyls_selected) >= 2) {
  wt <- mrg2g_wald(fit, corrections = "all")
  # Simulation only: attach the known null status to check the test
  wt$null_block <- d$design$null_blocks[match(wt$block_name, d$design$block_names)]
  cat("\nRegional Wald test\n")
  print(wt[, c("block_name", "null_block", "wald_stat", "df", "pval_raw",
               "padj_bh", "reject_bh")])
}