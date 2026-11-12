# =============================================================================
# MRG2G.R
#
# MR-G2G: two-stage Mendelian randomization for CpG-level causal effects with
# gene/region structure (individual-level data, continuous outcome).
#
# Required package: glmnet, MASS, Matrix
# =============================================================================


setup_optimization <- function(Z, X, ld_blocks, cg_blocks, lambda, a1, a2) {
  
  if (lambda <= 0) stop("lambda must be positive")
  if (sum(ld_blocks) != ncol(Z)) stop("Sum of blocks must equal number of columns in Z")
  if (sum(cg_blocks) != ncol(X)) stop("Sum of h_blocks must equal number of columns in X")
  
  n <- nrow(X)
  h <- ncol(X)
  p <- ncol(Z) 
  L <- length(ld_blocks)
  C <- length(cg_blocks)
  
  if (nrow(Z) != n) stop("Z and X must have same number of rows")
  
  ld_group_ind <- rep(1:L, times = ld_blocks) 
  
  cg_group_ind <- rep(1:C, times = cg_blocks)
  
  cat("Precomputing matrix operations...\n")
  ZtZ <- crossprod(Z)  # p x p matrix
  ZtX <- crossprod(Z, X)  # p x h matrix
  
  ld_group_indices <- vector("list", L) 
  for (l in 1:L) {
    ld_group_indices[[l]] <- which(ld_group_ind == l)
  }
  
  cg_group_indices <- vector("list", C)
  for (c in 1:C) {
    cg_group_indices[[c]] <- which(cg_group_ind == c)
  }
  
  B <- matrix(0, nrow = p, ncol = h)
  
  cat("Setup complete: p =", p, ", h =", h, ", L =", L, ", C =", C, "\n")
  cat("SNP block sizes:", ld_blocks, "\n")
  cat("Methylation block sizes:", cg_blocks, "\n")
  
  list(
    B = B,
    Z = Z,
    X = X,
    ZtZ = ZtZ,
    ZtX = ZtX,
    ld_blocks = ld_blocks,
    cg_blocks = cg_blocks,
    ld_group_ind = ld_group_ind,
    cg_group_ind = cg_group_ind,
    ld_group_indices = ld_group_indices,
    cg_group_indices = cg_group_indices,
    lambda = lambda,
    a1 = a1,
    a2 = a2,
    L = L,
    C = C,
    h = h,
    n = n,
    p = p
  )
}



optimize_proximal_gradient <- function(setup, max_iter = 5000, tol = 1e-6, 
                                               step_size = NULL, verbose = TRUE,
                                               n_inner = 3, decay_schedule = FALSE) {
  
  B <- setup$B
  Z <- setup$Z
  X <- setup$X
  ZtZ <- setup$ZtZ
  ZtX <- setup$ZtX
  ld_group_indices <- setup$ld_group_indices
  cg_group_indices <- setup$cg_group_indices
  ld_blocks <- setup$ld_blocks
  cg_blocks <- setup$cg_blocks
  lambda <- setup$lambda
  a1 <- setup$a1
  a2 <- setup$a2
  L <- setup$L
  h <- setup$h
  p <- setup$p
  
  if (is.null(step_size)) {
    L_smooth <- max(eigen(ZtZ, only.values = TRUE)$values)
    step_size <- 0.9 / L_smooth
  }
  
  if (verbose) {
    cat("Starting proximal gradient with cyclic Dykstra\n")
    if (decay_schedule) {
      cat("  Inner stopping: decay schedule epsilon_k = 1/k^3 (Schmidt et al.)\n")
    } else {
      cat("  Inner Dykstra cycles:", n_inner, "\n")
    }
    cat("  Step size:", step_size, "\n")
    cat("  Weights: a1 =", a1, ", a2 =", a2, 
        ", a3 =", round(1 - a1 - a2, 3), "\n")
  }
  
  objective_values <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    B_old <- B
    
    grad <- ZtZ %*% B - ZtX
    V <- B - step_size * grad  
    
    B <- dykstra_proximal(V, step_size, setup, n_cycles = n_inner,
                                 decay_schedule = decay_schedule,
                                 outer_iter = iter)
	
	objective_values[iter] <- 0.5 * norm(setup$X - setup$Z %*% B, type = "F")^2 +
  		compute_penalty(B, setup)
    
    frobenius_diff <- norm(B - B_old, type = "F")
    relative_change <- frobenius_diff / (norm(B, type = "F") + 1e-8)
    
    if (!is.na(frobenius_diff) && !is.na(relative_change) && 
        (frobenius_diff < tol || relative_change < tol * 10)) {
      if (verbose) cat("Converged at iteration", iter, "with ||B_new - B_old||_F =", 
                       frobenius_diff, "\n")
      objective_values <- objective_values[1:iter]
      break
    }
    
    if (is.na(frobenius_diff) || is.nan(frobenius_diff) ||
        is.na(relative_change) || is.nan(relative_change)) {
      if (verbose) cat("WARNING: Numerical issue at iteration", iter, "\n")
      objective_values <- objective_values[1:iter]
      break
    }
    
    if (verbose && iter %% 50 == 0) {
      cat(sprintf("Iteration %d: Obj = %.6f, ||d^t|| = %.2e\n",
                  iter, objective_values[iter], frobenius_diff))
    }
  }
  
  threshold_tiny <- 1e-8
  B_clean <- B
  B_clean[abs(B_clean) < threshold_tiny] <- 0
  active_coeffs_clean <- (abs(B_clean) > 1e-10)
  n_active <- sum(active_coeffs_clean)
  sparsity <- 1 - n_active / length(B_clean)
  
  if (verbose) {
    cat("Proximal gradient (CYCLIC Dykstra) optimization complete:\n")
    cat("  Active coefficients:", n_active, "out of", length(B), "\n")
    cat("  Sparsity level:", round(sparsity * 100, 2), "%\n")
    cat("  Final objective:", tail(objective_values, 1), "\n")
  }
  
  list(
    B = B_clean,
    B_raw = B,
    active_coeffs = active_coeffs_clean,
    iterations = iter,
    objective_values = objective_values,
    converged = (iter < max_iter),
    sparsity = sparsity,
    final_objective = tail(objective_values, 1)
  )
}


dykstra_proximal <- function(V, step_size, setup, n_cycles = 3,
                                    decay_schedule = FALSE, outer_iter = 1) {
  
  lambda          <- setup$lambda
  a1              <- setup$a1
  a2              <- setup$a2
  a3              <- round(1 - a1 - a2, 3)
  ld_group_indices <- setup$ld_group_indices
  cg_group_indices <- setup$cg_group_indices
  ld_blocks       <- setup$ld_blocks
  cg_blocks       <- setup$cg_blocks
  L               <- setup$L
  C               <- setup$C
  h               <- setup$h
  p               <- setup$p
  
  threshold_l1 <- step_size * lambda * a1
  
  x  <- V
  d1 <- matrix(0, nrow = nrow(V), ncol = ncol(V)) 
  d2 <- matrix(0, nrow = nrow(V), ncol = ncol(V)) 
  d3 <- matrix(0, nrow = nrow(V), ncol = ncol(V)) 
  
  if (decay_schedule) {
    epsilon_k  <- 1 / (outer_iter^3)  
    max_cycles <- 50                    
  } else {
    epsilon_k  <- 0
    max_cycles <- n_cycles
  }
  
  for (cycle in 1:max_cycles) {
    x_prev <- x
    
    # ---- Step 1: L1 prox (finest structure) ----
    if (a1 > 0) {
      y1 <- x + d1
      x  <- sign(y1) * pmax(abs(y1) - threshold_l1, 0)
      d1 <- y1 - x
    }
    
    # ---- Step 2: LD-group prox ----
    if (a2 > 0) {
      y2 <- x + d2
      x_new <- y2
      for (l in 1:L) {
        ld_idx      <- ld_group_indices[[l]]
        p_l         <- ld_blocks[l]
        thresh_ld   <- step_size * lambda * a2 * sqrt(p_l)
        for (k in 1:h) {
          g_coeffs  <- y2[ld_idx, k]
          g_norm    <- sqrt(sum(g_coeffs^2))
          if (g_norm <= thresh_ld) {
            x_new[ld_idx, k] <- 0
          } else {
            x_new[ld_idx, k] <- g_coeffs * (1 - thresh_ld / g_norm)
          }
        }
      }
      d2 <- y2 - x_new
      x  <- x_new
    }
    
    # ---- Step 3: CG-group prox (coarsest structure, acts last) ----
    if (a3 > 0) {
      y3    <- x + d3
      x_new <- y3
      for (j in 1:p) {
        for (c_idx in 1:C) {
          cg_idx    <- cg_group_indices[[c_idx]]
          h_c       <- cg_blocks[c_idx]
          thresh_cg <- step_size * lambda * a3 * sqrt(h_c)
          g_coeffs  <- y3[j, cg_idx]
          g_norm    <- sqrt(sum(g_coeffs^2))
          if (g_norm <= thresh_cg) {
            x_new[j, cg_idx] <- 0
          } else {
            x_new[j, cg_idx] <- g_coeffs * (1 - thresh_cg / g_norm)
          }
        }
      }
      d3 <- y3 - x_new
      x  <- x_new
    }
    
    # Decay schedule stopping criterion: stop when inner change < epsilon_k = 1/k^3
    if (decay_schedule) {
      if (norm(x - x_prev, type = "F") < epsilon_k) break
    }
  }
  
  return(x)
}


two_stage_optimize_multigene <- function(Z, X, gene_structure,
                                         lambda_cis, a1_cis, a2_cis,
                                         lambda_trans,
                                         max_iter = 5000,
                                         tol = 1e-6,
                                         verbose = TRUE,
                                         n_inner = 3,
                                         decay_schedule = TRUE) {
  
  G <- length(gene_structure)
  if (length(lambda_cis) != G) stop("lambda_cis must have length G")
  if (length(a1_cis) != G) stop("a1_cis must have length G")
  if (length(a2_cis) != G) stop("a2_cis must have length G")
  
  n <- nrow(Z)
  p <- ncol(Z)
  h <- ncol(X)
  
  B <- matrix(0, nrow = p, ncol = h)
  
  if (verbose) {
    cat("\n")
    cat(rep("-", 70), "\n", sep="")
    cat("TWO-STAGE MULTI-GENE OPTIMIZATION\n")
    cat(rep("-", 70), "\n", sep="")
    cat("Problem size: n =", n, ", p =", p, ", h =", h, "\n")
    cat("Number of genes:", G, "\n")
    cat("Stage 1: Independent cis estimation\n")
    cat("Stage 2: Trans estimation with glmnet\n")
    cat("Stage 3: Refined cis estimation\n")
    cat(rep("-", 70), "\n\n", sep="")
  }
  
  if (verbose) {
    cat("STAGE 1: Estimating cis blocks independently\n")
    cat(rep("-", 70), "\n", sep="")
  }
  
  stage1_results <- list()
  
  for (g in 1:G) {
    if (verbose) {
      cat(sprintf("\nGene %d cis block:\n", g))
    }
    
    snp_idx <- gene_structure[[g]]$snp_range
    methyl_idx <- gene_structure[[g]]$methyl_range
    ld_blocks_g <- gene_structure[[g]]$ld_blocks
    cg_blocks_g <- gene_structure[[g]]$cg_blocks
    
    Z_g <- Z[, snp_idx, drop = FALSE]
    X_g <- X[, methyl_idx, drop = FALSE]
    
    if (verbose) {
      cat(sprintf("  Dimensions: Z_g (%d × %d), X_g (%d × %d)\n", 
                  nrow(Z_g), ncol(Z_g), nrow(X_g), ncol(X_g)))
      cat(sprintf("  Penalties: lambda = %.2f, a1 = %.2f, a2 = %.2f\n",
                  lambda_cis[g], a1_cis[g], a2_cis[g]))
    }
    
    setup_g <- setup_optimization(
      Z = Z_g,
      X = X_g,
      ld_blocks = ld_blocks_g,
      cg_blocks = cg_blocks_g,
      lambda = lambda_cis[g],
      a1 = a1_cis[g],
      a2 = a2_cis[g]
    )
    
    result_g <- optimize_proximal_gradient(
        setup          = setup_g,
        max_iter       = max_iter,
        tol            = tol,
        verbose        = verbose,
        n_inner        = n_inner,
        decay_schedule = decay_schedule
      )

    
    B[snp_idx, methyl_idx] <- result_g$B
    stage1_results[[g]] <- result_g
    
    if (verbose) {
      cat(sprintf("  Gene %d: %d active coefficients, objective = %.6e\n",
                  g, sum(abs(result_g$B) > 1e-8), result_g$final_objective))
    }
  }
  
  if (verbose) {
    cat("\n")
    cat("Stage 1 complete!\n")
    cat(sprintf("Total cis non-zeros: %d\n", sum(abs(B) > 1e-8)))
    cat("\n")
  }
  
  if (verbose) {
    cat("STAGE 2: Estimating trans blocks with glmnet\n")
    cat(rep("=", 70), "\n", sep="")
  }
  
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for Stage 2. Please install it with: install.packages('glmnet')")
  }
  
  stage2_results <- list()
  trans_nnz_total <- 0
  
  for (f in 1:G) {
    if (verbose) {
      cat(sprintf("\nTarget Gene %d methylation:\n", f))
    }
    
    methyl_idx_f <- gene_structure[[f]]$methyl_range
    snp_idx_f <- gene_structure[[f]]$snp_range
    
    R_f <- X[, methyl_idx_f, drop = FALSE] - 
      Z[, snp_idx_f, drop = FALSE] %*% B[snp_idx_f, methyl_idx_f]
    
    for (g in 1:G) {
      if (g == f) next  # Skip cis block
      
      snp_idx_g <- gene_structure[[g]]$snp_range
      Z_g <- Z[, snp_idx_g, drop = FALSE]
      
      if (verbose) {
        cat(sprintf("  Estimating Gene %d → Gene %d trans block...", g, f))
      }
      
      B_gf <- matrix(0, nrow = length(snp_idx_g), ncol = length(methyl_idx_f))
      
      for (k in 1:ncol(R_f)) {
        fit <- glmnet::glmnet(
          x = Z_g,
          y = R_f[, k],
          alpha = 1,  # Lasso
          lambda = lambda_trans / n,   
          standardize = FALSE,  
          intercept = FALSE   
        )
        
        B_gf[, k] <- as.numeric(coef(fit, s = lambda_trans / n))[-1]  
      }
      
      B[snp_idx_g, methyl_idx_f] <- B_gf
      nnz_gf <- sum(abs(B_gf) > 1e-8)
      trans_nnz_total <- trans_nnz_total + nnz_gf
      
      if (verbose) {
        cat(sprintf(" %d non-zeros\n", nnz_gf))
      }
      
      R_f <- R_f - Z_g %*% B_gf
      
      stage2_results[[paste0("G", g, "toG", f)]] <- list(
        B_gf = B_gf,
        nnz = nnz_gf
      )
    }
  }
  
  if (verbose) {
    cat("\n")
    cat("Stage 2 complete!\n")
    cat(sprintf("Total trans non-zeros: %d\n", trans_nnz_total))
    cat("\n")
  }
  
  if (verbose) {
    cat("STAGE 3: Refining cis blocks given trans estimates\n")
    cat(rep("=", 70), "\n", sep="")
  }
  
  stage3_results <- list()
  
  for (g in 1:G) {
    if (verbose) {
      cat(sprintf("\nRefining Gene %d cis block:\n", g))
    }
    
    snp_idx_g <- gene_structure[[g]]$snp_range
    methyl_idx_g <- gene_structure[[g]]$methyl_range
    ld_blocks_g <- gene_structure[[g]]$ld_blocks
    cg_blocks_g <- gene_structure[[g]]$cg_blocks
    
    Z_g <- Z[, snp_idx_g, drop = FALSE]
    
    R_g <- X[, methyl_idx_g, drop = FALSE]
    
    for (f in 1:G) {
      if (f == g) next  
      
      snp_idx_f <- gene_structure[[f]]$snp_range
      R_g <- R_g - Z[, snp_idx_f, drop = FALSE] %*% B[snp_idx_f, methyl_idx_g]
    }
    
    if (verbose) {
      cat(sprintf("  Residual norm after removing trans: %.4f\n", norm(R_g, "F")))
    }
    
    setup_g <- setup_optimization(
      Z = Z_g,
      X = R_g,   
      ld_blocks = ld_blocks_g,
      cg_blocks = cg_blocks_g,
      lambda = lambda_cis[g],
      a1 = a1_cis[g],
      a2 = a2_cis[g]
    )
    
      result_g <- optimize_proximal_gradient(
        setup          = setup_g,
        max_iter       = max_iter,
        tol            = tol,
        verbose        = verbose,
        n_inner        = n_inner,
        decay_schedule = decay_schedule
      )
    
    B[snp_idx_g, methyl_idx_g] <- result_g$B
    stage3_results[[g]] <- result_g
    
    if (verbose) {
      cat(sprintf("  Gene %d refined: %d active coefficients, objective = %.6e\n",
                  g, sum(abs(result_g$B) > 1e-8), result_g$final_objective))
    }
  }
  
  if (verbose) {
    cat("\n")
    cat("Stage 3 complete!\n")
    cat(sprintf("Final total non-zeros: %d\n", sum(abs(B) > 1e-8)))
    cat("\n")
    cat(rep("=", 70), "\n", sep="")
    cat("TWO-STAGE OPTIMIZATION COMPLETE\n")
    cat(rep("=", 70), "\n", sep="")
  }
  
  B_clean <- B
  B_clean[abs(B_clean) < 1e-8] <- 0
  
  list(
    B = B_clean,
    B_raw = B,
    stage1_results = stage1_results,
    stage2_results = stage2_results,
    stage3_results = stage3_results,
    total_nnz = sum(abs(B_clean) > 1e-8),
    cis_nnz = sum(sapply(1:G, function(g) {
      sum(abs(B_clean[gene_structure[[g]]$snp_range, 
                      gene_structure[[g]]$methyl_range]) > 1e-8)
    })),
    trans_nnz = trans_nnz_total
  )
}


fit_scaled_group_lasso <- function(y, X_hat, cg_blocks,
									U 			 = NULL, 	
                                    lambda_scale = 1.0,
                                    sigma_init   = NULL,
                                    max_outer    = 50,
                                    max_inner    = 2000,
                                    tol_outer    = 1e-4,
                                    tol_inner    = 1e-6,
                                    warm_start   = TRUE,
                                    verbose      = TRUE) {
  y <- as.numeric(y)
  n <- nrow(X_hat)
  h <- ncol(X_hat)
  B <- length(cg_blocks)
  q <- if (is.null(U)) 0L else ncol(U)


  stopifnot(sum(cg_blocks) == h)
  stopifnot(length(y) == n)
  stopifnot(lambda_scale > 0)
  if (q > 0) stopifnot(nrow(U) == n)

  X_fit      <- if (q > 0) cbind(X_hat, U) else X_hat
  blocks_fit <- c(cg_blocks, rep(1L, q))
  
  block_end   <- cumsum(cg_blocks)
  block_start <- c(1L, head(block_end, -1L) + 1L)

  if (verbose) {
    cat("=== Scaled Group Lasso  ===\n")
    cat(sprintf("n = %d, h = %d, B = %d blocks\n", n, h, B))
    cat(sprintf("lambda_scale = %.4f\n", lambda_scale))
    cat(sprintf("Block sizes: %s\n",
                paste(cg_blocks, collapse = ", ")))
  }

  omega_base <- rep(NA_real_, h)
  for (b in seq_len(B)) {
    idx_b   <- block_start[b]:block_end[b]
    d_b     <- cg_blocks[b]
    omega_b <- sqrt(d_b / n) + sqrt(2.0 * log(B) / n)
    omega_base[idx_b] <- omega_b
  }

  if (verbose) {
    cat(sprintf("Base omega range: [%.4f, %.4f]\n",
                min(omega_base), max(omega_base)))
  }
  
  omega_base <- c(omega_base, rep(0, q))

  if (is.null(sigma_init)) {
    sigma_hat <- sd(y)
  } else {
    sigma_hat <- sigma_init
  }

  if (verbose) cat(sprintf("Initial sigma = %.4f\n\n", sigma_hat))

  theta_hat    <- rep(0.0, h + q)
  sigma_path   <- numeric(max_outer)
  inner_results <- vector("list", max_outer)

  for (outer_iter in seq_len(max_outer)) {

    if (verbose) cat(sprintf("Outer iter %d/%d: sigma = %.6f\n",
                              outer_iter, max_outer, sigma_hat))

    omega_scaled <- lambda_scale * sigma_hat * omega_base

    inner_res <- solve_group_lasso_proxgrad(
      y          = y,
      X_hat      = X_fit,
      omega      = omega_scaled,
      cg_blocks  = blocks_fit,
      max_iter   = max_inner,
      tol        = tol_inner,
      theta_init = if (warm_start) theta_hat else NULL,
      verbose    = verbose && (outer_iter == 1) 
    )

    theta_hat          <- inner_res$theta
    inner_results[[outer_iter]] <- inner_res

    resid_outer  <- y - X_fit %*% theta_hat
    sigma_new    <- sqrt(sum(resid_outer^2) / n)

    sigma_path[outer_iter] <- sigma_new

    if (verbose) {
      p_active <- sum(abs(theta_hat[seq_len(h)]) > 1e-10)
      cat(sprintf("  sigma updated: %.6f -> %.6f | active theta: %d/%d\n",
                  sigma_hat, sigma_new, p_active, h))
    }


    rel_change_sigma <- abs(sigma_new - sigma_hat) / (abs(sigma_hat) + 1e-10)

    sigma_hat <- sigma_new

    if (rel_change_sigma < tol_outer) {
      if (verbose) cat(sprintf("\nConverged at outer iter %d (rel delta sigma = %.2e)\n",
                                outer_iter, rel_change_sigma))
      sigma_path <- sigma_path[seq_len(outer_iter)]
      inner_results <- inner_results[seq_len(outer_iter)]

      theta_hat[abs(theta_hat) < 1e-10] <- 0.0
      theta_conf <- theta_hat[h + seq_len(q)]
      theta_hat  <- theta_hat[seq_len(h)]
      
      n_active     <- sum(abs(theta_hat) > 1e-10)
      active_blocks <- which(sapply(seq_len(B), function(b) {
        idx_b <- block_start[b]:block_end[b]
        any(abs(theta_hat[idx_b]) > 1e-10)
      }))

      if (verbose) {
        cat(sprintf("Final sigma_hat = %.6f\n", sigma_hat))
        cat(sprintf("Active theta entries: %d / %d\n", n_active, h))
        cat(sprintf("Active blocks: %s\n",
                    if (length(active_blocks) > 0)
                      paste(active_blocks, collapse = ", ")
                    else "none"))
      }

      return(list(
      	theta_hat     = theta_hat,
	    theta_conf    = theta_conf,
    	sigma_hat     = sigma_hat,
	    omega_final   = omega_scaled,
    	omega_base    = omega_base,
	    sigma_path    = sigma_path,
    	outer_iters   = outer_iter,
	    converged     = TRUE,
    	inner_results = inner_results,
    	n_active      = n_active,
    	active_blocks = active_blocks
      	))
    }
  }

  warning(sprintf("fit_scaled_group_lasso: max_outer (%d) reached without convergence. ",
                  max_outer),
          sprintf("Final rel change in sigma = %.2e. Consider increasing max_outer.",
                  rel_change_sigma))

  theta_hat[abs(theta_hat) < 1e-10] <- 0.0
  theta_conf <- theta_hat[h + seq_len(q)]
  theta_hat  <- theta_hat[seq_len(h)]
  n_active <- sum(abs(theta_hat) > 1e-10)
  active_blocks <- which(sapply(seq_len(B), function(b) {
    idx_b <- block_start[b]:block_end[b]
    any(abs(theta_hat[idx_b]) > 1e-10)
  }))

  return(list(
    theta_hat     = theta_hat,
    theta_conf	  = theta_conf, 	
    sigma_hat     = sigma_hat,
    omega_final   = omega_scaled,
    omega_base    = omega_base,
    sigma_path    = sigma_path,
    outer_iters   = max_outer,
    converged     = FALSE,
    inner_results = inner_results,
    n_active      = n_active,
    active_blocks = active_blocks
  ))
}



solve_group_lasso_proxgrad <- function(y, X_hat, omega, cg_blocks,
                                        max_iter  = 2000,
                                        tol       = 1e-6,
                                        theta_init = NULL,
                                        verbose   = FALSE) {

  n <- nrow(X_hat)
  h <- ncol(X_hat)
  B <- length(cg_blocks)

  stopifnot(sum(cg_blocks) == h)
  stopifnot(length(omega) == h)

  block_end   <- cumsum(cg_blocks)
  block_start <- c(1L, head(block_end, -1L) + 1L)

  XtX <- crossprod(X_hat)         
  Xty <- crossprod(X_hat, y)        


  step_size <- 1.0 / (max(eigen(XtX, symmetric = TRUE,
                                only.values = TRUE)$values) / n + 1e-10)

  if (is.null(theta_init)) {
    theta <- rep(0.0, h)
  } else {
    theta <- theta_init
  }

  objective_values <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {

    grad    <- -(Xty - XtX %*% theta) / n    
    theta_g <- theta - step_size * grad     

    theta_new <- theta_g  

    for (b in seq_len(B)) {
      idx_b   <- block_start[b]:block_end[b]
      v_b     <- theta_g[idx_b]
      norm_vb <- sqrt(sum(v_b^2))
      omega_b <- omega[block_start[b]]   

      threshold_b <- step_size * omega_b

      if (norm_vb <= threshold_b) {
        theta_new[idx_b] <- 0.0
      } else {
        theta_new[idx_b] <- v_b * (1.0 - threshold_b / norm_vb)
      }
    }


    resid_iter  <- y - X_hat %*% theta_new
    smooth_part <- sum(resid_iter^2) / (2.0 * n)
    penalty_part <- sum(sapply(seq_len(B), function(b) {
      idx_b   <- block_start[b]:block_end[b]
      omega_b <- omega[block_start[b]]
      omega_b * sqrt(sum(theta_new[idx_b]^2))
    }))
    objective_values[iter] <- smooth_part + penalty_part

    delta_theta <- sqrt(sum((theta_new - theta)^2))
    theta       <- theta_new

    if (delta_theta < tol) {
      if (verbose) cat(sprintf("    [inner prox] Converged at iter %d, delta=%.2e\n",
                                iter, delta_theta))
      objective_values <- objective_values[seq_len(iter)]
      return(list(
        theta            = theta,
        objective_values = objective_values,
        iterations       = iter,
        converged        = TRUE
      ))
    }
  }

  if (verbose) cat(sprintf("    [inner prox] Max iter (%d) reached\n", max_iter))
  return(list(
    theta            = theta,
    objective_values = objective_values,
    iterations       = max_iter,
    converged        = FALSE
  ))
}


.build_na_result <- function(block_names, cg_blocks,
                              sigma_hat,
                              estimator_used, en_collapsed,
                              active_corrections, fdr_level, B) {
  out <- data.frame(
    block_name        = block_names,
    n_selected        = cg_blocks,
    wald_stat         = NA_real_,
    df                = NA_integer_,
    pval_raw          = NA_real_,
    lambda_scale_used = NA_real_,
    sigma_hat         = sigma_hat,
    estimator_used    = estimator_used,
    en_collapsed      = en_collapsed,
    stringsAsFactors  = FALSE
  )
  if ("none"        %in% active_corrections) {
    out$padj_none         <- NA_real_
    out$reject_none       <- NA_integer_
  }
  if ("bh"          %in% active_corrections) {
    out$padj_bh           <- NA_real_
    out$reject_bh         <- NA_integer_
  }
  if ("bonferroni"  %in% active_corrections) {
    out$padj_bonferroni   <- NA_real_
    out$reject_bonferroni <- NA_integer_
  }
  out
}


compute_penalty <- function(B, setup) {

  lambda           <- setup$lambda
  a1               <- setup$a1
  a2               <- setup$a2
  a3               <- round(1 - a1 - a2, 3)
  ld_group_indices <- setup$ld_group_indices
  cg_group_indices <- setup$cg_group_indices
  ld_blocks        <- setup$ld_blocks
  cg_blocks        <- setup$cg_blocks
  L                <- setup$L
  C                <- setup$C

  l1_pen <- a1 * lambda * sum(abs(B))

  ld_pen <- 0
  for (l in 1:L) {
    idx <- ld_group_indices[[l]]
    p_l <- ld_blocks[l]
    for (k in 1:ncol(B)) {
      ld_pen <- ld_pen + sqrt(p_l) * sqrt(sum(B[idx, k]^2))
    }
  }
  ld_pen <- a2 * lambda * ld_pen

  cg_pen <- 0
  for (j in 1:nrow(B)) {
    for (c in 1:C) {
      idx <- cg_group_indices[[c]]
      h_c <- cg_blocks[c]
      cg_pen <- cg_pen + sqrt(h_c) * sqrt(sum(B[j, idx]^2))
    }
  }
  cg_pen <- a3 * lambda * cg_pen

  l1_pen + ld_pen + cg_pen
}


construct_score_matrix <- function(X_G, X_neg_G, lambda = NULL,
                                   cond_threshold = 30) {

  n    <- nrow(X_G)
  d_b  <- ncol(X_G)
  p_neg <- ncol(X_neg_G)

  Gamma_hat <- matrix(0, nrow = p_neg, ncol = d_b)

  if (p_neg > 0) {
    if (is.null(lambda)) {
      lambda <- sqrt(log(p_neg) / n)
    }

    for (j in 1:d_b) {
      fit_j <- glmnet::glmnet(x = X_neg_G, y = X_G[, j],
                      alpha = 1,
                      lambda = lambda,
                      intercept = FALSE,
                      standardize = FALSE)
      Gamma_hat[, j] <- as.numeric(coef(fit_j)[-1])
    }
  }

  if (p_neg > 0) {
    Z_G <- X_G - X_neg_G %*% Gamma_hat
  } else {
    Z_G <- X_G
  }

  svd_ZG <- svd(Z_G)
  singular_values <- svd_ZG$d

  nonzero_sv <- singular_values[singular_values > 1e-10]
  cond_num   <- max(nonzero_sv) / min(nonzero_sv)

  svd_truncated <- FALSE

  if (cond_num > cond_threshold) {
    sv_threshold <- singular_values[1] / sqrt(cond_threshold)
    r <- sum(singular_values > sv_threshold)
    r <- max(r, 1)
    svd_truncated <- TRUE
    cat(sprintf("  [score matrix] Block collinearity detected (cond = %.1f). ",
                cond_num))
    cat(sprintf("Truncating to r = %d of %d directions.\n", r, d_b))
  } else {
    r <- length(nonzero_sv)
  }

  U_r <- svd_ZG$u[, 1:r, drop = FALSE]
  P_G <- U_r %*% t(U_r)

  return(list(
    Z_G           = Z_G,
    P_G           = P_G,
    U_r           = U_r,
    k_G           = r,
    Gamma_hat     = Gamma_hat,
    singular_values = singular_values,
    cond_num      = cond_num,
    svd_truncated = svd_truncated
  ))
}



mz_group_wald_mr <- function(y, X_hat, cg_blocks,
                              U                   = NULL,
                              lambda_score        = NULL,
                              fdr_level           = 0.05,
                              cond_threshold      = 30,
                              theta_init_override = NULL,
                              gamma_init_override = NULL,
                              sigma_hat_override  = NULL,
                              verbose             = TRUE) {

  n <- nrow(X_hat)
  h <- ncol(X_hat)
  q <- if (is.null(U)) 0L else ncol(U)
  B <- length(cg_blocks)

  stopifnot(sum(cg_blocks) == h)
  stopifnot(length(y) == n)
  if (q > 0) stopifnot(nrow(U) == n)

  if (is.null(theta_init_override) || is.null(sigma_hat_override))
    stop("mz_group_wald_mr: theta_init_override and sigma_hat_override are required.")
  if (q > 0 && is.null(gamma_init_override))
    stop("mz_group_wald_mr: gamma_init_override is required when U is supplied.")

  stopifnot(length(theta_init_override) == h)
  if (q > 0) stopifnot(length(gamma_init_override) == q)
  stopifnot(is.numeric(sigma_hat_override) && length(sigma_hat_override) == 1)
  stopifnot(sigma_hat_override > 0)

  theta_init <- theta_init_override
  if (!is.null(colnames(X_hat))) names(theta_init) <- colnames(X_hat)
  gamma_init <- if (q > 0) gamma_init_override else numeric(0)
  sigma_hat  <- sigma_hat_override

  block_end   <- cumsum(cg_blocks)
  block_start <- c(1, head(block_end, -1) + 1)

  if (verbose) {
    cat("=== Mitra & Zhang Group Wald Test (MR Version) ===\n")
    cat(sprintf("n = %d, h = %d methylation sites, q = %d confounders, B = %d blocks\n",
                n, h, q, B))
    cat(sprintf("  sigma_hat = %.4f, active exposures = %d\n",
                sigma_hat, sum(abs(theta_init) > 1e-10)))
    cat("\nStep 3: Computing block-level Wald statistics (MR)...\n")
  }

  results    <- data.frame(block      = seq_len(B),
                           block_size = cg_blocks,
                           wald_stat  = NA_real_,
                           df         = NA_integer_,
                           pval       = NA_real_,
                           padj       = NA_real_,
                           cond_num   = NA_real_,
                           svd_trunc  = NA)
  score_info <- vector("list", B)

  for (b in seq_len(B)) {

    if (verbose) cat(sprintf("  Block %d/%d (size %d)... ", b, B, cg_blocks[b]))

    if (cg_blocks[b] == 0) {
      if (verbose) cat("EMPTY (no sites selected) -- skipped\n")
      next
    }

    idx_b     <- block_start[b]:block_end[b]
    idx_neg_b <- setdiff(seq_len(h), idx_b)

    X_hat_b     <- X_hat[, idx_b,     drop = FALSE]
    X_hat_neg_b <- X_hat[, idx_neg_b, drop = FALSE]
    X_neg_b_aug <- if (q > 0) cbind(X_hat_neg_b, U) else X_hat_neg_b

    sm_mr <- construct_score_matrix(X_hat_b, X_neg_b_aug,
                                     lambda         = lambda_score,
                                     cond_threshold = cond_threshold)
    score_info[[b]] <- sm_mr

    P_b <- sm_mr$P_G
    k_b <- sm_mr$k_G

    theta_neg_b_mr  <- theta_init[idx_neg_b]
    partial_resid_b <- y - X_hat_neg_b %*% theta_neg_b_mr
    if (q > 0) partial_resid_b <- partial_resid_b - U %*% gamma_init

    proj_b_mr   <- P_b %*% partial_resid_b
    wald_stat_b <- sum(proj_b_mr^2) / sigma_hat^2
    pval_b      <- pchisq(wald_stat_b, df = k_b, lower.tail = FALSE)

    results$wald_stat[b] <- wald_stat_b
    results$df[b]        <- k_b
    results$pval[b]      <- pval_b
    results$cond_num[b]  <- sm_mr$cond_num
    results$svd_trunc[b] <- sm_mr$svd_truncated

    if (verbose) cat(sprintf("W = %.3f, df = %d, p = %.4f\n",
                              wald_stat_b, k_b, pval_b))
  }

  results$padj <- p.adjust(results$pval, method = "BH")

  if (verbose) {
    cat(sprintf("\nBH correction at FDR = %.2f:\n", fdr_level))
    cat(sprintf("  Blocks rejected: %d / %d\n",
                sum(results$padj < fdr_level, na.rm = TRUE), B))
  }

  return(list(
    results    = results,
    theta_init = theta_init,
    gamma_init = gamma_init,
    sigma_hat  = sigma_hat,
    score_info = score_info
  ))
}


mz_group_wald_test <- function(y,
                                X_hat,
                                cg_blocks,
                                U                = NULL,
                                block_names      = NULL,
                                null_blocks      = NULL,
                                estimator        = "sgl",
                                corrections      = "bh",
                                lambda_scale     = 0.5,
                                en_fallback      = "warn",
                                fdr_level        = 0.05,
                                cond_threshold   = 30,
                                verbose          = FALSE) {

  y   <- as.numeric(y)
  n   <- nrow(X_hat)
  h   <- ncol(X_hat)
  q   <- if (is.null(U)) 0L else ncol(U)
  B   <- length(cg_blocks)

  if (length(y) != n)
    stop("mz_group_wald_test: length(y) must equal nrow(X_hat).")

  if (q > 0 && nrow(U) != n)
    stop("mz_group_wald_test: nrow(U) must equal nrow(X_hat).")

  if (sum(cg_blocks) != h)
    stop(sprintf(
      "mz_group_wald_test: sum(cg_blocks) = %d but ncol(X_hat) = %d. ",
      sum(cg_blocks), h))

  if (!estimator %in% c("sgl", "en"))
    stop("mz_group_wald_test: estimator must be 'sgl' or 'en'.")

  if (!all(corrections %in% c("none", "bh", "bonferroni", "all")))
    stop("mz_group_wald_test: corrections must be 'none', 'bh', 'bonferroni', or 'all'.")

  if (!en_fallback %in% c("warn", "error", "sgl"))
    stop("mz_group_wald_test: en_fallback must be 'warn', 'error', or 'sgl'.")

  if (is.null(block_names)) {
    block_names <- paste0("Block", seq_len(B))
  } else if (length(block_names) != B) {
    stop("mz_group_wald_test: block_names must have the same length as cg_blocks.")
  }

  if ("all" %in% corrections) {
    active_corrections <- c("none", "bh", "bonferroni")
  } else {
    active_corrections <- unique(corrections)
  }

  lambda_scale_used <- as.numeric(lambda_scale)
  if (verbose) cat(sprintf("  Fixed lambda_scale = %.1f\n", lambda_scale_used))

  estimator_used <- estimator
  en_collapsed   <- NA

  if (estimator == "sgl") {

    sgl_fit <- tryCatch({
      fit_scaled_group_lasso(
        y            = y,
        X_hat        = X_hat,
        cg_blocks    = cg_blocks,
        U            = U,
        lambda_scale = lambda_scale_used,
        sigma_init   = sd(y),
        max_outer    = 50,
        max_inner    = 2000,
        tol_outer    = 1e-4,
        tol_inner    = 1e-6,
        warm_start   = TRUE,
        verbose      = FALSE
      )
    }, error = function(e) {
      stop(sprintf("mz_group_wald_test: fit_scaled_group_lasso failed: %s", e$message))
    })

    theta_init <- sgl_fit$theta_hat
    gamma_init <- sgl_fit$theta_conf
    sigma_hat  <- sgl_fit$sigma_hat

    if (verbose) {
      cat(sprintf("  SGL: sigma=%.4f | active=%d / %d | q=%d confounders\n",
                  sigma_hat, sgl_fit$n_active, h, q))
    }

  } else {

    if (!requireNamespace("glmnet", quietly = TRUE))
      stop("mz_group_wald_test: package 'glmnet' is required for estimator = 'en'.")

    X_en           <- if (q > 0) cbind(X_hat, U) else X_hat
    penalty_factor <- c(rep(1, h), rep(0, q))

    en_fit <- tryCatch({
      alpha_grid <- seq(0.1, 0.9, by = 0.1)
      cv_mse     <- numeric(length(alpha_grid))
      for (i in seq_along(alpha_grid)) {
        cv_tmp    <- glmnet::cv.glmnet(X_en, y, alpha = alpha_grid[i],
                                       penalty.factor = penalty_factor)
        cv_mse[i] <- min(cv_tmp$cvm)
      }
      best_alpha <- alpha_grid[which.min(cv_mse)]
      cv_best    <- glmnet::cv.glmnet(X_en, y, alpha = best_alpha,
                                      penalty.factor = penalty_factor)
      coefs_en   <- as.numeric(coef(cv_best, s = "lambda.min")[-1])
      list(theta_hat  = coefs_en[seq_len(h)],
           gamma_hat  = coefs_en[h + seq_len(q)],
           best_alpha = best_alpha)
    }, error = function(e) {
      stop(sprintf("mz_group_wald_test: elastic net fitting failed: %s", e$message))
    })

    theta_en     <- en_fit$theta_hat
    gamma_en     <- en_fit$gamma_hat
    en_collapsed <- sum(abs(theta_en)) < 1e-10

    if (en_collapsed) {
      if (en_fallback == "error") {
        stop("mz_group_wald_test: EN collapsed to zero vector (en_fallback = 'error').")
      } else if (en_fallback == "warn") {
        warning("mz_group_wald_test: EN collapsed to zero vector. Returning NA Wald stats.")
        out <- .build_na_result(
          block_names, cg_blocks, 
          NA_real_, "en", TRUE,
          active_corrections, fdr_level, B
        )
        return(out)
      } else {
        warning("mz_group_wald_test: EN collapsed to zero vector. Falling back to SGL.")
        sgl_fb <- tryCatch({
          fit_scaled_group_lasso(
            y            = y,
            X_hat        = X_hat,
            cg_blocks    = cg_blocks,
            U            = U,
            lambda_scale = lambda_scale_used,
            sigma_init   = sd(y),
            max_outer    = 50,
            max_inner    = 2000,
            tol_outer    = 1e-4,
            tol_inner    = 1e-6,
            warm_start   = TRUE,
            verbose      = FALSE
          )
        }, error = function(e) {
          stop(sprintf(
            "mz_group_wald_test: SGL fallback also failed: %s", e$message))
        })
        theta_init     <- sgl_fb$theta_hat
        gamma_init     <- sgl_fb$theta_conf
        sigma_hat      <- sgl_fb$sigma_hat
        estimator_used <- "sgl"
      }
    } else {
      theta_init <- theta_en
      gamma_init <- gamma_en
      resid_en   <- y - X_hat %*% theta_en
      if (q > 0) {
        resid_en  <- resid_en - U %*% gamma_en
        sigma_hat <- sqrt(sum(resid_en^2) / max(n - sum(abs(theta_en) > 1e-10) - q, 1))
      } else {
        sigma_hat <- sqrt(sum(resid_en^2) / n)
      }
    }

    if (verbose) {
      cat(sprintf("  EN: alpha=%.1f | active=%d / %d | collapsed=%s | sigma=%.4f\n",
                  en_fit$best_alpha, sum(abs(theta_init) > 1e-10),
                  h, en_collapsed, sigma_hat))
    }
  }

  wald_out <- tryCatch({
    mz_group_wald_mr(
      y                   = y,
      X_hat               = X_hat,
      cg_blocks           = cg_blocks,
      U                   = U,
      theta_init_override = theta_init,
      gamma_init_override = gamma_init,
      sigma_hat_override  = sigma_hat,
      fdr_level           = fdr_level,
      cond_threshold      = cond_threshold,
      verbose             = FALSE
    )
  }, error = function(e) {
    stop(sprintf("mz_group_wald_test: mz_group_wald_mr failed: %s", e$message))
  })

  wr        <- wald_out$results
  raw_pvals <- wr$pval

  if (verbose) {
    cat(sprintf("  Wald stats: [%s]\n",
                paste(round(wr$wald_stat, 2), collapse = ", ")))
    cat(sprintf("  Raw p-vals: [%s]\n",
                paste(round(raw_pvals, 4), collapse = ", ")))
  }

  out <- data.frame(
    block_name        = block_names,
    n_selected        = cg_blocks,
    wald_stat         = wr$wald_stat,
    df                = wr$df,
    pval_raw          = raw_pvals,
    lambda_scale_used = if (estimator_used == "en") NA_real_ else lambda_scale_used,
    sigma_hat         = sigma_hat,
    estimator_used    = estimator_used,
    en_collapsed      = en_collapsed,
    stringsAsFactors  = FALSE
  )



  if ("none" %in% active_corrections) {
    out$padj_none    <- raw_pvals
    out$reject_none  <- as.integer(!is.na(raw_pvals) & raw_pvals < fdr_level)
  }

  if ("bh" %in% active_corrections) {
    padj_bh         <- p.adjust(raw_pvals, method = "BH")
    out$padj_bh     <- padj_bh
    out$reject_bh   <- as.integer(!is.na(padj_bh) & padj_bh < fdr_level)
  }

  if ("bonferroni" %in% active_corrections) {
    padj_bonf             <- p.adjust(raw_pvals, method = "bonferroni")
    out$padj_bonferroni   <- padj_bonf
    out$reject_bonferroni <- as.integer(!is.na(padj_bonf) & padj_bonf < fdr_level)
  }

  return(out)
}


make_gene_structure <- function(ld_blocks, cg_blocks) {
  if (!is.list(ld_blocks)) ld_blocks <- list(ld_blocks)
  if (!is.list(cg_blocks)) cg_blocks <- list(cg_blocks)
  if (length(ld_blocks) != length(cg_blocks))
    stop("make_gene_structure: ld_blocks and cg_blocks must have one entry per gene.")

  snp_end    <- cumsum(vapply(ld_blocks, sum, numeric(1)))
  methyl_end <- cumsum(vapply(cg_blocks, sum, numeric(1)))
  snp_start    <- c(1, head(snp_end, -1) + 1)
  methyl_start <- c(1, head(methyl_end, -1) + 1)

  lapply(seq_along(ld_blocks), function(g) list(
    snp_range    = snp_start[g]:snp_end[g],
    methyl_range = methyl_start[g]:methyl_end[g],
    ld_blocks    = ld_blocks[[g]],
    cg_blocks    = cg_blocks[[g]]
  ))
}


.mrg2g_standardize <- function(M) {
  M  <- as.matrix(M)
  mu <- colMeans(M)
  s  <- apply(M, 2, sd)
  s[s == 0] <- 1
  out <- sweep(sweep(M, 2, mu, "-"), 2, s, "/")
  attr(out, "center") <- mu
  attr(out, "scale")  <- s
  out
}


mrg2g_setup <- function(Z, X, y, U = NULL, gene_structure,
                        region_names = NULL, standardize = TRUE) {

  Z <- as.matrix(Z)
  X <- as.matrix(X)
  y <- as.numeric(y)
  if (!is.null(U)) U <- as.matrix(U)

  n <- nrow(Z)
  if (nrow(X) != n || length(y) != n || (!is.null(U) && nrow(U) != n))
    stop("mrg2g_setup: Z, X, y and U must have the same number of observations.")

  ld_all <- unlist(lapply(gene_structure, `[[`, "ld_blocks"))
  cg_all <- unlist(lapply(gene_structure, `[[`, "cg_blocks"))
  if (sum(ld_all) != ncol(Z))
    stop("mrg2g_setup: total ld_blocks in gene_structure must equal ncol(Z).")
  if (sum(cg_all) != ncol(X))
    stop("mrg2g_setup: total cg_blocks in gene_structure must equal ncol(X).")
  for (g in seq_along(gene_structure)) {
    gs <- gene_structure[[g]]
    if (length(gs$snp_range) != sum(gs$ld_blocks) ||
        length(gs$methyl_range) != sum(gs$cg_blocks))
      stop(sprintf("mrg2g_setup: gene %d ranges do not match its block sizes.", g))
  }

  if (is.null(colnames(Z))) colnames(Z) <- paste0("snp", seq_len(ncol(Z)))
  if (is.null(colnames(X))) colnames(X) <- paste0("cg", seq_len(ncol(X)))
  if (!is.null(U) && is.null(colnames(U))) colnames(U) <- paste0("U", seq_len(ncol(U)))

  if (is.null(region_names)) {
    region_names <- paste0("Region", seq_along(cg_all))
  } else if (length(region_names) != length(cg_all)) {
    stop("mrg2g_setup: region_names must have one entry per CpG region (cg_blocks).")
  }

  scales <- NULL
  if (standardize) {
    Z <- .mrg2g_standardize(Z)
    X <- .mrg2g_standardize(X)
    if (!is.null(U)) U <- .mrg2g_standardize(U)
    y_s <- .mrg2g_standardize(y)
    scales <- list(Z = attr(Z, "scale"), X = attr(X, "scale"),
                   U = if (!is.null(U)) attr(U, "scale") else NULL,
                   y = attr(y_s, "scale"))
    y <- as.numeric(y_s)
  }

  structure(list(
    Z              = Z,
    X              = X,
    y              = y,
    U              = U,
    gene_structure = gene_structure,
    ld_blocks      = ld_all,
    cg_blocks      = cg_all,
    site_region_id = rep(seq_along(cg_all), times = cg_all),
    region_names   = region_names,
    standardized   = standardize,
    scales         = scales
  ), class = "mrg2g_data")
}


mrg2g_fit <- function(dat,
                      lambda, a1, a2,
                      lambda_trans   = NULL,
                      stage2         = "en",
                      lambda_scale   = 0.5,
                      max_iter       = 5000,
                      tol            = 1e-6,
                      n_inner        = 3,
                      decay_schedule = TRUE,
                      verbose        = FALSE) {

  stopifnot(inherits(dat, "mrg2g_data"))
  stage2 <- unique(stage2)
  if (!all(stage2 %in% c("ols", "lasso", "en", "sgl")))
    stop("mrg2g_fit: stage2 must be any of 'ols', 'lasso', 'en', 'sgl'.")

  Z <- dat$Z; X <- dat$X; y <- dat$y; U <- dat$U
  gs <- dat$gene_structure
  G  <- length(gs)
  q  <- if (is.null(U)) 0L else ncol(U)
  h  <- ncol(X)
  n_obs  <- nrow(Z)
  lambda <- rep_len(lambda, G) * n_obs  
  a1     <- rep_len(a1, G)
  a2     <- rep_len(a2, G)

  if (G == 1) {
    setup1 <- setup_optimization(
      Z = Z[, gs[[1]]$snp_range, drop = FALSE],
      X = X[, gs[[1]]$methyl_range, drop = FALSE],
      ld_blocks = gs[[1]]$ld_blocks, cg_blocks = gs[[1]]$cg_blocks,
      lambda = lambda, a1 = a1, a2 = a2)
    stage1 <- optimize_proximal_gradient(
      setup1, max_iter = max_iter, tol = tol, verbose = verbose,
      n_inner = n_inner, decay_schedule = decay_schedule)
    B_hat <- matrix(0, ncol(Z), h, dimnames = list(colnames(Z), colnames(X)))
    B_hat[gs[[1]]$snp_range, gs[[1]]$methyl_range] <- stage1$B
  } else {
    if (is.null(lambda_trans))
      stop("mrg2g_fit: lambda_trans is required when there is more than one gene.")
    stage1 <- two_stage_optimize_multigene(
      Z, X, gs,
      lambda_cis = lambda, a1_cis = a1, a2_cis = a2,
      lambda_trans = lambda_trans * n_obs,
      max_iter = max_iter, tol = tol, verbose = verbose,
      n_inner = n_inner, decay_schedule = decay_schedule)
    B_hat <- stage1$B
    dimnames(B_hat) <- list(colnames(Z), colnames(X))
  }

  IVs_selected     <- apply(B_hat, 1, function(r) sqrt(sum(r^2)) > 0)
  methyls_selected <- apply(B_hat, 2, function(cl) sqrt(sum(cl^2)) > 0)

  region_counts <- as.integer(table(factor(dat$site_region_id[methyls_selected],
                                           levels = seq_along(dat$cg_blocks))))
  regions_kept  <- region_counts > 0
  cg_blocks_sel <- region_counts[regions_kept]

  out <- list(
    B_hat            = B_hat,
    IVs_selected     = IVs_selected,
    methyls_selected = methyls_selected,
    cg_blocks_sel    = cg_blocks_sel,
    regions_kept     = regions_kept,
    region_names_sel = dat$region_names[regions_kept],
    X_hat            = NULL,
    estimates        = list(),
    stage2           = stage2,
    lambda           = lambda / n_obs,
    lambda_trans     = lambda_trans,
    lambda_scale     = lambda_scale,
    stage1           = stage1,
    data             = dat
  )

  if (sum(IVs_selected) == 0 || sum(methyls_selected) == 0) {
    warning("mrg2g_fit: Stage 1 selected no instruments or no CpG sites; Stage 2 skipped.")
    return(structure(out, class = "mrg2g_fit"))
  }

  Z_sel <- Z[, IVs_selected, drop = FALSE]
  X_sel <- X[, methyls_selected, drop = FALSE]
  X_hat <- matrix(0, nrow = nrow(X_sel), ncol = ncol(X_sel),
                  dimnames = list(NULL, colnames(X_sel)))
  for (j in seq_len(ncol(X_sel))) {
    X_hat[, j] <- fitted(lm(X_sel[, j] ~ Z_sel))
  }
  out$X_hat <- X_hat
  h_sel <- ncol(X_hat)

  D  <- if (q > 0) cbind(X_hat, U) else X_hat
  pf <- c(rep(1, h_sel), rep(0, q))

  to_full <- function(coefs) {
    theta <- setNames(rep(0, h), colnames(X))
    theta[colnames(X_hat)] <- coefs[seq_len(h_sel)]
    gamma <- if (q > 0) setNames(coefs[h_sel + seq_len(q)], colnames(U)) else numeric(0)
    list(theta = theta, gamma = gamma)
  }

  if ("ols" %in% stage2) {
    ols_fit <- lm(y ~ D)
    out$estimates$ols <- c(to_full(as.numeric(coef(ols_fit))[-1]), list(fit = ols_fit))
  }

  if ("lasso" %in% stage2) {
    lasso_fit <- glmnet::cv.glmnet(D, y, alpha = 1, penalty.factor = pf)
    out$estimates$lasso <- c(
      to_full(as.numeric(coef(lasso_fit, s = "lambda.min"))[-1]),
      list(fit = lasso_fit))
  }

  if ("en" %in% stage2) {
    alpha_grid  <- seq(0.1, 0.9, by = 0.1)
    cv_mse_vals <- sapply(alpha_grid, function(a)
      min(glmnet::cv.glmnet(D, y, alpha = a, penalty.factor = pf)$cvm))
    best_alpha  <- alpha_grid[which.min(cv_mse_vals)]
    en_fit      <- glmnet::cv.glmnet(D, y, alpha = best_alpha, penalty.factor = pf)
    out$estimates$en <- c(
      to_full(as.numeric(coef(en_fit, s = "lambda.min"))[-1]),
      list(fit = en_fit, alpha = best_alpha))
  }

  if ("sgl" %in% stage2) {
    sgl_fit <- fit_scaled_group_lasso(
      y = y, X_hat = X_hat, cg_blocks = cg_blocks_sel, U = U,
      lambda_scale = lambda_scale, sigma_init = sd(y),
      max_outer = 50, max_inner = 2000, tol_outer = 1e-4, tol_inner = 1e-6,
      warm_start = TRUE, verbose = FALSE)
    out$estimates$sgl <- c(
      to_full(c(sgl_fit$theta_hat, sgl_fit$theta_conf)),
      list(fit = sgl_fit, sigma_hat = sgl_fit$sigma_hat))
  }

  if (dat$standardized) {
    for (m in names(out$estimates)) {
      out$estimates[[m]]$theta_raw <-
        out$estimates[[m]]$theta * dat$scales$y / dat$scales$X
    }
  }

  structure(out, class = "mrg2g_fit")
}


mrg2g_wald <- function(fit,
                       estimator      = "sgl",
                       corrections    = "bh",
                       lambda_scale   = fit$lambda_scale,
                       fdr_level      = 0.05,
                       cond_threshold = 30,
                       verbose        = FALSE) {

  stopifnot(inherits(fit, "mrg2g_fit"))
  if (is.null(fit$X_hat) || ncol(fit$X_hat) < 2)
    stop("mrg2g_wald: at least two CpG sites must be selected in Stage 1.")

  mz_group_wald_test(
    y              = fit$data$y,
    X_hat          = fit$X_hat,
    cg_blocks      = fit$cg_blocks_sel,
    U              = fit$data$U,
    block_names    = fit$region_names_sel,
    estimator      = estimator,
    corrections    = corrections,
    lambda_scale   = lambda_scale,
    fdr_level      = fdr_level,
    cond_threshold = cond_threshold,
    verbose        = verbose
  )
}


.mrg2g_quiet <- function(expr) {
  utils::capture.output(val <- expr)
  val
}


.mrg2g_cis_fit <- function(Z, X, ld_blocks, cg_blocks, lambda, a1, a2,
                           max_iter, tol, n_inner, decay_schedule) {
  .mrg2g_quiet({
    s <- setup_optimization(Z, X, ld_blocks = ld_blocks, cg_blocks = cg_blocks,
                            lambda = lambda * nrow(Z), a1 = a1, a2 = a2)
    optimize_proximal_gradient(s, max_iter = max_iter, tol = tol, verbose = FALSE,
                               n_inner = n_inner, decay_schedule = decay_schedule)
  })$B
}


.mrg2g_boundary_check <- function(value, grid, what) {
  if (value %in% range(grid))
    warning(sprintf("mrg2g_tune: %s at grid boundary (%.3f); widen the grid.", what, value),
            call. = FALSE)
}


mrg2g_tune <- function(dat,
                       lambda_grid       = seq(0.01, 0.40, by = 0.03),
                       a_grid            = NULL,
                       lambda_trans_grid = seq(0.05, 1.00, by = 0.05),
                       K                 = 5,
                       seed              = 2026,
                       cores             = 1,
                       max_iter          = 5000,
                       tol               = 1e-6,
                       n_inner           = 3,
                       decay_schedule    = TRUE,
                       verbose           = TRUE) {

  stopifnot(inherits(dat, "mrg2g_data"))

  if (is.null(a_grid)) {
    a_grid <- expand.grid(a1 = c(0.1, 0.2, 0.4, 0.6),
                          a2 = seq(0, 0.8, by = 0.2))
    a_grid <- a_grid[a_grid$a1 + a_grid$a2 <= 1 + 1e-9, , drop = FALSE]
  }
  if (any(a_grid$a1 >= 1))
    stop("mrg2g_tune: a_grid must not contain a1 = 1 (no structured penalty).")

  Z  <- dat$Z
  X  <- dat$X
  gs <- dat$gene_structure
  G  <- length(gs)
  n  <- nrow(Z)

  old_seed <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  fold <- sample(rep_len(seq_len(K), n))
  if (is.null(old_seed)) rm(".Random.seed", envir = globalenv())
  else assign(".Random.seed", old_seed, envir = globalenv())

  run <- function(jobs, FUN) {
    if (cores > 1) parallel::mclapply(jobs, FUN, mc.cores = cores)
    else lapply(jobs, FUN)
  }


  cis_tables <- vector("list", G)
  best_cis   <- data.frame(gene = seq_len(G), lambda = NA_real_,
                           a1 = NA_real_, a2 = NA_real_)

  for (g in seq_len(G)) {
    Zg <- Z[, gs[[g]]$snp_range,    drop = FALSE]
    Xg <- X[, gs[[g]]$methyl_range, drop = FALSE]

    tab <- expand.grid(lambda = lambda_grid, row = seq_len(nrow(a_grid)))
    tab$a1  <- a_grid$a1[tab$row]
    tab$a2  <- a_grid$a2[tab$row]
    tab$row <- NULL

    jobs <- expand.grid(i = seq_len(nrow(tab)), k = seq_len(K))
    if (verbose) cat(sprintf("Gene %d: %d cis fits (%d settings x %d folds)\n",
                             g, nrow(jobs), nrow(tab), K))

    err <- run(seq_len(nrow(jobs)), function(r) {
      i  <- jobs$i[r]; k <- jobs$k[r]
      tr <- fold != k
      B  <- .mrg2g_cis_fit(Zg[tr, , drop = FALSE], Xg[tr, , drop = FALSE],
                           gs[[g]]$ld_blocks, gs[[g]]$cg_blocks,
                           tab$lambda[i], tab$a1[i], tab$a2[i],
                           max_iter, tol, n_inner, decay_schedule)
      mean((Xg[!tr, , drop = FALSE] - Zg[!tr, , drop = FALSE] %*% B)^2)
    })
    E <- matrix(unlist(err), nrow = nrow(tab), ncol = K)

    tab$cvm <- rowMeans(E)
    tab$se  <- apply(E, 1, sd) / sqrt(K)

    j_min  <- which.min(tab$cvm)
    same_a <- which(tab$a1 == tab$a1[j_min] & tab$a2 == tab$a2[j_min])
    elig   <- same_a[tab$cvm[same_a] <= tab$cvm[j_min] + tab$se[j_min]]
    j      <- elig[which.max(tab$lambda[elig])]

    best_cis[g, c("lambda", "a1", "a2")] <- tab[j, c("lambda", "a1", "a2")]
    cis_tables[[g]] <- tab[order(tab$cvm), ]

    .mrg2g_boundary_check(tab$lambda[j_min], lambda_grid, sprintf("gene %d CV-min lambda", g))
    .mrg2g_boundary_check(tab$lambda[j],     lambda_grid, sprintf("gene %d 1-SE lambda", g))
    if (verbose) cat(sprintf("  a1=%.2f a2=%.2f lambda=%.3f\n",
                             tab$a1[j], tab$a2[j], tab$lambda[j]))
  }

  lambda_trans <- NULL
  trans_table  <- NULL

  if (G > 1) {
    jobs <- expand.grid(i = seq_along(lambda_trans_grid), k = seq_len(K))
    if (verbose) cat(sprintf("Trans: %d multigene fits\n", nrow(jobs)))

    err <- run(seq_len(nrow(jobs)), function(r) {
      i    <- jobs$i[r]; k <- jobs$k[r]
      tr   <- fold != k
      n_tr <- sum(tr)
      res  <- .mrg2g_quiet(two_stage_optimize_multigene(
        Z[tr, , drop = FALSE], X[tr, , drop = FALSE], gs,
        lambda_cis     = best_cis$lambda * n_tr,
        a1_cis         = best_cis$a1,
        a2_cis         = best_cis$a2,
        lambda_trans   = lambda_trans_grid[i] * n_tr,
        max_iter       = max_iter,
        tol            = tol,
        verbose        = FALSE,
        n_inner        = n_inner,
        decay_schedule = decay_schedule))
      mean((X[!tr, , drop = FALSE] - Z[!tr, , drop = FALSE] %*% res$B)^2)
    })
    E <- matrix(unlist(err), nrow = length(lambda_trans_grid), ncol = K)

    trans_table <- data.frame(lambda_trans = lambda_trans_grid,
                              cvm = rowMeans(E),
                              se  = apply(E, 1, sd) / sqrt(K))
    j_min <- which.min(trans_table$cvm)
    elig  <- which(trans_table$cvm <= trans_table$cvm[j_min] + trans_table$se[j_min])
    j     <- elig[which.max(trans_table$lambda_trans[elig])]
    lambda_trans <- trans_table$lambda_trans[j]

    .mrg2g_boundary_check(trans_table$lambda_trans[j_min], lambda_trans_grid, "CV-min lambda_trans")
    .mrg2g_boundary_check(lambda_trans,                    lambda_trans_grid, "1-SE lambda_trans")
    if (verbose) cat(sprintf("  lambda_trans=%.3f\n", lambda_trans))
  }

  list(lambda       = best_cis$lambda,
       a1           = best_cis$a1,
       a2           = best_cis$a2,
       lambda_trans = lambda_trans,
       best_cis     = best_cis,
       cis_tables   = cis_tables,
       trans_table  = trans_table,
       folds        = fold)
}