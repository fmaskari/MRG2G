

#Setup optimization environment for multi-gene scenario
# @Z: Matrix of SNP data (n x p), values in {0, 1, 2}
# @X: Matrix of methylation data (n x h)
# @blocks: Vector of block sizes for SNP LD groups
# @h_blocks: Vector of block sizes for methylation groups
# @lambda: Penalty parameter (> 0)
# @a1: L1 penalty weight (0 <= a1 <= 1)
# @a2: Group L2 penalty weight (0 <= a2 <= 1, a1 + a2 <= 1)
setup_optimization <- function(Z, X, ld_blocks, cg_blocks, lambda, a1, a2) {
  
  # Validation
  if (lambda <= 0) stop("lambda must be positive")
  #if (a1 < 0 || a2 < 0 || (a1 + a2) > 1) stop("Invalid penalty weights")
  if (sum(ld_blocks) != ncol(Z)) stop("Sum of blocks must equal number of columns in Z")
  if (sum(cg_blocks) != ncol(X)) stop("Sum of h_blocks must equal number of columns in X")
  
  n <- nrow(X)
  h <- ncol(X)
  p <- ncol(Z) 
  L <- length(ld_blocks)
  C <- length(cg_blocks)
  
  # Check dimensions
  if (nrow(Z) != n) stop("Z and X must have same number of rows")
  
  # Create group indicator for SNP LD blocks
  ld_group_ind <- rep(1:L, times = ld_blocks) #each snp marked with which LD group it belongs to
  
  # Create group indicator for methylation blocks
  cg_group_ind <- rep(1:C, times = cg_blocks)
  
  # Precompute expensive operations
  cat("Precomputing matrix operations...\n")
  ZtZ <- crossprod(Z)  # p x p matrix
  ZtX <- crossprod(Z, X)  # p x h matrix
  
  # Create group index lists for SNP LD blocks
  ld_group_indices <- vector("list", L) #each sublist includes list of SNPs belonging to that group
  for (l in 1:L) {
    ld_group_indices[[l]] <- which(ld_group_ind == l)
  }
  
  # Create group index lists for methylation blocks
  cg_group_indices <- vector("list", C)
  for (c in 1:C) {
    cg_group_indices[[c]] <- which(cg_group_ind == c)
  }
  
  # Initialize coefficient matrix (all zeros)
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



#' Soft thresholding operator
soft_threshold <- function(x, threshold) {
  sign(x) * pmax(abs(x) - threshold, 0)
}

#' Compute the full objective function
compute_objective <- function(B, Z, X, lambda, a1, a2, 
						ld_group_indices, ld_blocks,
						cg_group_indices, cg_blocks) {
  
  p <- nrow(B)
  h <- ncol(B)
  L <- length(ld_group_indices)
  C <- length(cg_group_indices)
  
  # Prediction error term
  X_pred <- Z %*% B
  mse_term <- 0.5 * norm(X - X_pred, type = "F")^2
  
  # L1 penalty term
  l1_penalty <- 0
  for (j in 1:nrow(B)) {
    l1_penalty <- l1_penalty + sum(abs(B[j, ]))
  }
  l1_penalty <- a1 * lambda * l1_penalty
  
  # LD Group L2 penalty term
  ld_group_penalty <- 0
  for (l in 1:L) {
    ld_group_idx <- ld_group_indices[[l]]
    p_l <- ld_blocks[l]
    for (k in 1:ncol(B)) {
      ld_group_penalty <- ld_group_penalty + sqrt(p_l) * norm(B[ld_group_idx, k], type = "2")
    }
  }
  ld_group_penalty <- a2 * lambda * ld_group_penalty
  
  # Functional CG Block L2 penalty term
  cg_block_penalty <- 0
  for (j in 1:nrow(B)) {
  	for (c in 1:C){
  		cg_idx <- cg_group_indices[[c]]
  		h_c <- cg_blocks[c]
  		cg_block_penalty <- cg_block_penalty + sqrt(h_c) * norm(B[j, cg_idx], type = "2")
  	}
  	}
  cg_block_penalty <- (1 - a1 - a2) * lambda * cg_block_penalty
  
  return(mse_term + l1_penalty + ld_group_penalty + cg_block_penalty)
}


# ============================================================================
# PROXIMAL GRADIENT DESCENT ALGORITHM (ALTERNATIVE TO COORDINATE DESCENT)
# ============================================================================

#' Proximal gradient descent optimization
#' @param setup Setup object from setup_optimization()
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param step_size Fixed step size (or use backtracking)
#' @param verbose Print progress information
optimize_proximal_gradient <- function(setup, max_iter = 5000, tol = 1e-6, 
                                      step_size = NULL, verbose = TRUE) {
  
  # Extract setup components
  B <- setup$B
  Z <- setup$Z
  X <- setup$X
  ZtZ <- setup$ZtZ
  ZtX <- setup$ZtX
  ld_group_indices <- setup$ld_group_indices  #list the size of L; each sublist includes list of SNPs belonging to that ld group
  cg_group_indices <- setup$cg_group_indices
  ld_blocks <- setup$ld_blocks
  cg_blocks <- setup$cg_blocks
  lambda <- setup$lambda
  a1 <- setup$a1
  a2 <- setup$a2
  L <- setup$L
  h <- setup$h
  p <- setup$p
  
  # Compute step size if not provided (based on Lipschitz constant)
  if (is.null(step_size)) {
    L_smooth <- max(eigen(ZtZ, only.values = TRUE)$values)
    step_size <- 0.9 / L_smooth  # Conservative step size
  }
  
  if (verbose) cat("Starting proximal gradient optimization with step size:", step_size, "\n")
  
  objective_values <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    B_old <- B
    
    # Step 1: Gradient step on smooth part
    grad <- ZtZ %*% B - ZtX  # Gradient of (1/2)||X - ZB||²_F
    B_temp <- B - step_size * grad
    
    # Step 2: Apply composite proximal operator
    B <- apply_composite_proximal(B_temp, step_size, setup)
    
    # Compute objective
    objective_values[iter] <- compute_objective(B, Z, X, lambda, a1, a2, 
    							ld_group_indices, ld_blocks,
    							cg_group_indices, cg_blocks)
    			 					
    
    # Check convergence
    frobenius_diff <- norm(B - B_old, type = "F")
    relative_change <- frobenius_diff / (norm(B, type = "F") + 1e-8)
    
    if (frobenius_diff < tol || relative_change < tol * 10) {
      if (verbose) cat("Converged at iteration", iter, "with ||B_new - B_old||_F =", 
                      frobenius_diff, "\n")
      objective_values <- objective_values[1:iter]
      break
    }
    
    if (verbose && iter %% 50 == 0) {
      cat("Iteration", iter, ": Objective =", objective_values[iter], 
          ", ||Delta B||_F =", frobenius_diff, "\n")
    }
  }
  
  # Final results
  threshold_tiny <- 1e-8
  B_clean <- B
  B_clean[abs(B_clean) < threshold_tiny] <- 0
  active_coeffs_clean <- (abs(B_clean) > 1e-10)
  n_active <- sum(active_coeffs_clean)
  sparsity <- 1 - n_active / length(B_clean)
  
  if (verbose) {
    cat("Proximal gradient optimization complete:\n")
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

#' Apply composite proximal operator for all penalties
apply_composite_proximal <- function(B, step_size, setup) {
  lambda <- setup$lambda
  a1 <- setup$a1
  a2 <- setup$a2
  ld_group_indices <- setup$ld_group_indices
  cg_group_indices <- setup$cg_group_indices
  ld_blocks <- setup$ld_blocks
  cg_blocks <- setup$cg_blocks
  L <- setup$L
  C <- setup$C
  h <- setup$h
  p <- setup$p
  
  # Step 1: L1 proximal operator (soft thresholding)
  threshold_l1 <- step_size * a1 * lambda
  B1 <- sign(B) * pmax(abs(B) - threshold_l1, 0)
  
  # Step 2: LD Group L2 proximal operator
  B2 <- B1
  for (l in 1:L) {
    ld_group_idx <- ld_group_indices[[l]]
    p_l <- ld_blocks[l]
    threshold_ld_group <- step_size * a2 * lambda * sqrt(p_l)
    
    for (k in 1:h) {
      ld_group_coeffs <- B1[ld_group_idx, k]
      ld_col_group_norm <- norm(ld_group_coeffs, type = "2")
      
      if (ld_col_group_norm <= threshold_ld_group) {
        B2[ld_group_idx, k] <- 0
      } else {
        B2[ld_group_idx, k] <- ld_group_coeffs * (1 - threshold_ld_group / ld_col_group_norm)
      }
    }
  }
  
  # Step 3: Functional CG Block L2 proximal operator
  B3 <- B2
  a3 <- round(1 - a1 - a2, 3)
  
  for (j in 1:p) {
  	for (c in 1:C){
  	  cg_idx <- cg_group_indices[[c]]
      h_c <- cg_blocks[c]
      threshold_cg_group <- step_size * a3 * lambda * sqrt(h_c)
      cg_group_coeffs <- B2[j, cg_idx]
      cg_row_group_norm <- norm(cg_group_coeffs, type = "2")
      
      if (cg_row_group_norm <= threshold_cg_group) {
        B3[j, cg_idx] <- 0
      } else {
        B3[j, cg_idx] <- cg_group_coeffs * (1 - threshold_cg_group / cg_row_group_norm)
      }
  	}
    
  }
  
  return(B3)
}

#' Two-Stage Optimization for Multi-Gene Sparse Regression
#' 
#' Stage 1: Estimate cis blocks independently using single-gene optimization
#' Stage 2: Estimate trans blocks using glmnet given fixed cis blocks
#' Stage 3: Refine cis blocks given trans estimates
#' 
#' @param Z Matrix of SNP data (n x p)
#' @param X Matrix of methylation data (n x h)
#' @param gene_structure List of gene structures (snp_range, methyl_range, ld_blocks, cg_blocks)
#' @param lambda_cis Vector of lambda values for each gene's cis block
#' @param a1_cis Vector of a1 values for each gene
#' @param a2_cis Vector of a2 values for each gene
#' @param lambda_trans Penalty parameter for trans blocks (single value)
#' @param max_iter Maximum iterations for single-gene optimization
#' @param tol Convergence tolerance
#' @param verbose Print progress
#' 
#' @return List containing optimized B matrix and diagnostics
two_stage_optimize_multigene <- function(Z, X, gene_structure,
                                         lambda_cis, a1_cis, a2_cis,
                                         lambda_trans,
                                         max_iter = 5000,
                                         tol = 1e-6,
                                         verbose = TRUE,
                                         use_dykstra = TRUE,
                                         n_inner = 3,
                                         decay_schedule = TRUE) {
  
  # Validate inputs
  G <- length(gene_structure)
  if (length(lambda_cis) != G) stop("lambda_cis must have length G")
  if (length(a1_cis) != G) stop("a1_cis must have length G")
  if (length(a2_cis) != G) stop("a2_cis must have length G")
  
  n <- nrow(Z)
  p <- ncol(Z)
  h <- ncol(X)
  
  # Initialize full B matrix
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
    cat("Cis solver: ", if (use_dykstra) paste0("Cyclic Dykstra (n_inner = ", n_inner, 
                                                ", decay_schedule = ", decay_schedule, ")") else "Sequential proximal", "\n")
    cat(rep("-", 70), "\n\n", sep="")
  }
  
  # ========================================================================
  # STAGE 1: INDEPENDENT CIS OPTIMIZATION
  # ========================================================================
  
  if (verbose) {
    cat("STAGE 1: Estimating cis blocks independently\n")
    cat(rep("-", 70), "\n", sep="")
  }
  
  stage1_results <- list()
  
  for (g in 1:G) {
    if (verbose) {
      cat(sprintf("\nGene %d cis block:\n", g))
    }
    
    # Extract gene-specific data
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
    
    # Setup single-gene optimization
    setup_g <- setup_optimization(
      Z = Z_g,
      X = X_g,
      ld_blocks = ld_blocks_g,
      cg_blocks = cg_blocks_g,
      lambda = lambda_cis[g],
      a1 = a1_cis[g],
      a2 = a2_cis[g]
    )
    
    # Optimize — use cyclic Dykstra or sequential proximal based on flag
    if (use_dykstra) {
      result_g <- optimize_proximal_gradient_dykstra(
        setup          = setup_g,
        max_iter       = max_iter,
        tol            = tol,
        verbose        = verbose,
        n_inner        = n_inner,
        decay_schedule = decay_schedule
      )
    } else {
      result_g <- optimize_proximal_gradient(
        setup    = setup_g,
        max_iter = max_iter,
        tol      = tol,
        verbose  = verbose
      )
    }
    
    # Store in full B matrix
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
  
  # ========================================================================
  # STAGE 2: TRANS ESTIMATION WITH GLMNET
  # ========================================================================
  
  if (verbose) {
    cat("STAGE 2: Estimating trans blocks with glmnet\n")
    cat(rep("=", 70), "\n", sep="")
  }
  
  # Check if glmnet is available
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for Stage 2. Please install it with: install.packages('glmnet')")
  }
  
  stage2_results <- list()
  trans_nnz_total <- 0
  
  # For each target gene
  for (f in 1:G) {
    if (verbose) {
      cat(sprintf("\nTarget Gene %d methylation:\n", f))
    }
    
    methyl_idx_f <- gene_structure[[f]]$methyl_range
    snp_idx_f <- gene_structure[[f]]$snp_range
    
    # Compute initial residual (remove cis contribution)
    R_f <- X[, methyl_idx_f, drop = FALSE] - 
      Z[, snp_idx_f, drop = FALSE] %*% B[snp_idx_f, methyl_idx_f]
    
    # For each source gene
    for (g in 1:G) {
      if (g == f) next  # Skip cis block
      
      snp_idx_g <- gene_structure[[g]]$snp_range
      Z_g <- Z[, snp_idx_g, drop = FALSE]
      
      if (verbose) {
        cat(sprintf("  Estimating Gene %d → Gene %d trans block...", g, f))
      }
      
      # Solve Lasso for each methylation site in gene f
      B_gf <- matrix(0, nrow = length(snp_idx_g), ncol = length(methyl_idx_f))
      
      for (k in 1:ncol(R_f)) {
        # Fit Lasso with glmnet
        # Note: glmnet uses alpha = 1 for Lasso, and lambda is the penalty
        # glmnet parameterization: (1/2n)||y - Xb||^2 + lambda||b||_1
        # Our parameterization: (1/2)||y - Xb||^2 + lambda||b||_1
        # So we need to scale lambda by n
        
        fit <- glmnet::glmnet(
          x = Z_g,
          y = R_f[, k],
          alpha = 1,  # Lasso
          lambda = lambda_trans / n,  # Scale by sample size
          standardize = FALSE,  # Don't standardize (we want raw coefficients)
          intercept = FALSE  # No intercept
        )
        
        # Extract coefficients (glmnet returns sparse matrix)
        B_gf[, k] <- as.numeric(coef(fit, s = lambda_trans / n))[-1]  # Remove intercept
      }
      
      # Store trans block
      B[snp_idx_g, methyl_idx_f] <- B_gf
      nnz_gf <- sum(abs(B_gf) > 1e-8)
      trans_nnz_total <- trans_nnz_total + nnz_gf
      
      if (verbose) {
        cat(sprintf(" %d non-zeros\n", nnz_gf))
      }
      
      # Update residual (sequential updates)
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
  
  # ========================================================================
  # STAGE 3: REFINE CIS BLOCKS
  # ========================================================================
  
  if (verbose) {
    cat("STAGE 3: Refining cis blocks given trans estimates\n")
    cat(rep("=", 70), "\n", sep="")
  }
  
  stage3_results <- list()
  
  for (g in 1:G) {
    if (verbose) {
      cat(sprintf("\nRefining Gene %d cis block:\n", g))
    }
    
    # Extract gene-specific indices
    snp_idx_g <- gene_structure[[g]]$snp_range
    methyl_idx_g <- gene_structure[[g]]$methyl_range
    ld_blocks_g <- gene_structure[[g]]$ld_blocks
    cg_blocks_g <- gene_structure[[g]]$cg_blocks
    
    Z_g <- Z[, snp_idx_g, drop = FALSE]
    
    # Compute residual: X_g minus all trans contributions TO gene g
    R_g <- X[, methyl_idx_g, drop = FALSE]
    
    for (f in 1:G) {
      if (f == g) next  # Skip cis
      
      snp_idx_f <- gene_structure[[f]]$snp_range
      R_g <- R_g - Z[, snp_idx_f, drop = FALSE] %*% B[snp_idx_f, methyl_idx_g]
    }
    
    if (verbose) {
      cat(sprintf("  Residual norm after removing trans: %.4f\n", norm(R_g, "F")))
    }
    
    # Setup single-gene optimization with residual
    setup_g <- setup_optimization(
      Z = Z_g,
      X = R_g,  # Use residual instead of raw X_g
      ld_blocks = ld_blocks_g,
      cg_blocks = cg_blocks_g,
      lambda = lambda_cis[g],
      a1 = a1_cis[g],
      a2 = a2_cis[g]
    )
    
    # Re-optimize — use cyclic Dykstra or sequential proximal based on flag
    if (use_dykstra) {
      result_g <- optimize_proximal_gradient_dykstra(
        setup          = setup_g,
        max_iter       = max_iter,
        tol            = tol,
        verbose        = verbose,
        n_inner        = n_inner,
        decay_schedule = decay_schedule
      )
    } else {
      result_g <- optimize_proximal_gradient(
        setup    = setup_g,
        max_iter = max_iter,
        tol      = tol,
        verbose  = verbose
      )
    }
    
    # Update B with refined cis block
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
  
  # Clean up final B
  B_clean <- B
  B_clean[abs(B_clean) < 1e-8] <- 0
  
  # Return results
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



# ============================================================================
# PROXIMAL GRADIENT WITH CYCLIC DYKSTRA INNER SOLVER (Algorithm 10.31)
# Following Combettes & Pesquet (2011), Proposition 10.32
# ============================================================================

# Proximal gradient descent with CYCLIC Dykstra inner solver
# @param setup Setup object from setup_optimization()
# @param max_iter Maximum number of outer iterations
# @param tol Convergence tolerance
# @param step_size Fixed step size (or computed from Lipschitz constant)
# @param verbose Print progress information
# @param n_inner Number of Dykstra inner iterations per outer iteration
optimize_proximal_gradient_dykstra <- function(setup, max_iter = 5000, tol = 1e-6, 
                                               step_size = NULL, verbose = TRUE,
                                               n_inner = 3, decay_schedule = FALSE) {
  
  # Extract setup components
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
  
  # Compute step size if not provided (based on Lipschitz constant)
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
    
    # Step 1: Gradient step on smooth part
    grad <- ZtZ %*% B - ZtX
    V <- B - step_size * grad  # V = B^(t) - gamma * grad f(B^(t))
    
    # Step 2: Approximate prox_{gamma * u}(V) via cyclic Dykstra
    B <- cyclic_dykstra_proximal(V, step_size, setup, n_cycles = n_inner,
                                 decay_schedule = decay_schedule,
                                 outer_iter = iter)
    
    # Compute objective
    objective_values[iter] <- compute_objective(B, Z, X, lambda, a1, a2, 
                                                ld_group_indices, ld_blocks,
                                                cg_group_indices, cg_blocks)
    
    # Check convergence
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
  
  # Final results
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

# Cyclic Dykstra proximal algorithm
#
# Solves: min_x f_1(x) + f_2(x) + f_3(x) + (1/2gamma)||x - V||^2
# where f_1 = a1*lambda*||.||_1, f_2 = a2*lambda*||.||_LD, f_3 = a3*lambda*||.||_CG
#
# Cyclic ordering following Yu (2013) Corollary 3:
#   finer structure first, coarser structure last: L1 -> LD -> CG
# Each prox feeds its output directly into the next (sequential composition).
# Dykstra increment variables (d1, d2, d3) carry the dual corrections
# across cycles, exactly as in classical Dykstra (Boyle & Dykstra 1986).
#
# @param V Input matrix (gradient step output)
# @param step_size Outer step size gamma
# @param setup Optimization setup object
# @param n_cycles Number of cyclic passes
# @return Approximate prox_{gamma*(f1+f2+f3)}(V)
cyclic_dykstra_proximal <- function(V, step_size, setup, n_cycles = 3,
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
  
  # Thresholds: gamma * lambda * a_i * (group weight)
  threshold_l1 <- step_size * lambda * a1
  
  # Current iterate and Dykstra increment variables (one per penalty)
  x  <- V
  d1 <- matrix(0, nrow = nrow(V), ncol = ncol(V))  # L1 increment
  d2 <- matrix(0, nrow = nrow(V), ncol = ncol(V))  # LD-group increment
  d3 <- matrix(0, nrow = nrow(V), ncol = ncol(V))  # CG-group increment
  
  # Determine stopping rule
  if (decay_schedule) {
    epsilon_k  <- 1 / (outer_iter^3)  # Schmidt et al. theoretically grounded schedule
    max_cycles <- 50                   # safety cap to avoid infinite inner loop
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



# ============================================================================
# REGION WALD TEST HELPER FUNCTIONS
# ============================================================================

fit_scaled_group_lasso <- function(y, X_hat, cg_blocks,
                                    lambda_scale = 1.0,
                                    sigma_init   = NULL,
                                    max_outer    = 50,
                                    max_inner    = 2000,
                                    tol_outer    = 1e-4,
                                    tol_inner    = 1e-6,
                                    warm_start   = TRUE,
                                    verbose      = TRUE) {

  n <- nrow(X_hat)
  h <- ncol(X_hat)
  B <- length(cg_blocks)

  stopifnot(sum(cg_blocks) == h)
  stopifnot(length(y) == n)
  stopifnot(lambda_scale > 0)

  block_end   <- cumsum(cg_blocks)
  block_start <- c(1L, head(block_end, -1L) + 1L)

  if (verbose) {
    cat("=== Scaled Group Lasso (Mitra & Zhang eq 2.25 / Algorithm 3.14) ===\n")
    cat(sprintf("n = %d, h = %d, B = %d blocks\n", n, h, B))
    cat(sprintf("lambda_scale = %.4f\n", lambda_scale))
    cat(sprintf("Block sizes: %s\n",
                paste(cg_blocks, collapse = ", ")))
  }

  # ------------------------------------------------------------------
  # Compute BASE penalty weights: omega_base_b = sqrt(d_b/n) + sqrt(2 log B / n)
  #
  # These match the theoretically motivated weights in Theorem 7:
  #   omega_{*,j} proportional to (sqrt(|G_j|) + sqrt(2 log M)) / sqrt(n)
  #
  # The base weights are BEFORE sigma scaling. After each outer iteration,
  # we rescale by sigma_hat: omega_b = lambda_scale * sigma_hat * omega_base_b
  #
  # We store as an h-vector (one entry per site, same value within a block)
  # for compatibility with the inner solver.
  # ------------------------------------------------------------------
  omega_base <- rep(NA_real_, h)
  for (b in seq_len(B)) {
    idx_b   <- block_start[b]:block_end[b]
    d_b     <- cg_blocks[b]
    # Paper's omega_{*,j} = sqrt(d_b/n) + sqrt(2 log B / n)
    omega_b <- sqrt(d_b / n) + sqrt(2.0 * log(B) / n)
    omega_base[idx_b] <- omega_b
  }

  if (verbose) {
    cat(sprintf("Base omega range: [%.4f, %.4f]\n",
                min(omega_base), max(omega_base)))
  }

  # ------------------------------------------------------------------
  # Initialize sigma
  #
  # The paper suggests starting from a reasonable sigma estimate.
  # Using sd(y) is conservative but safe. Alternatively, use
  # ||y||_2 / sqrt(n) (oracle if theta* = 0) or residual sd from
  # an OLS fit (if p << n).
  # ------------------------------------------------------------------
  if (is.null(sigma_init)) {
    sigma_hat <- sd(y)
  } else {
    sigma_hat <- sigma_init
  }

  if (verbose) cat(sprintf("Initial sigma = %.4f\n\n", sigma_hat))

  # ------------------------------------------------------------------
  # Initialize theta at zero (or warm start will update this)
  # ------------------------------------------------------------------
  theta_hat    <- rep(0.0, h)
  sigma_path   <- numeric(max_outer)
  inner_results <- vector("list", max_outer)

  # ------------------------------------------------------------------
  # OUTER LOOP: Algorithm 3.14
  #
  # Each iteration:
  #   (a) Rescale weights: omega' = lambda_scale * sigma_hat * omega_base
  #   (b) Solve inner group lasso with fixed omega'
  #   (c) Update sigma: sigma_hat_new = ||y - X_hat theta_hat||_2 / sqrt(n)
  #   (d) Check convergence of sigma
  #
  # Lemma 2 of the paper proves this algorithm converges to the joint
  # minimizer of the scaled objective (3.12). The profile loss function
  # is convex in sigma (eq 3.15), ensuring the outer loop is well-behaved.
  # ------------------------------------------------------------------
  for (outer_iter in seq_len(max_outer)) {

    if (verbose) cat(sprintf("Outer iter %d/%d: sigma = %.6f\n",
                              outer_iter, max_outer, sigma_hat))

    # ----------------------------------------------------------------
    # Step (a): Rescale penalty weights by current sigma_hat
    #           omega' = lambda_scale * sigma_hat * omega_base
    #           (eq 3.14: omega' <- sigma_hat^{(k+1)} * omega)
    # ----------------------------------------------------------------
    omega_scaled <- lambda_scale * sigma_hat * omega_base

    # ----------------------------------------------------------------
    # Step (b): Solve inner group lasso with fixed omega'
    #           theta_hat <- argmin_theta { (1/2n)||y - X_hat theta||^2
    #                                       + sum_b omega'_b ||theta_b||_2 }
    #           (eq 3.14: beta_hat^{(k+1)} <- argmin_beta L_{omega'}(beta))
    # ----------------------------------------------------------------
    inner_res <- solve_group_lasso_proxgrad(
      y          = y,
      X_hat      = X_hat,
      omega      = omega_scaled,
      cg_blocks  = cg_blocks,
      max_iter   = max_inner,
      tol        = tol_inner,
      theta_init = if (warm_start) theta_hat else NULL,
      verbose    = verbose && (outer_iter == 1)  # only verbose on first iter
    )

    theta_hat          <- inner_res$theta
    inner_results[[outer_iter]] <- inner_res

    # ----------------------------------------------------------------
    # Step (c): Update sigma estimate
    #           sigma_hat <- ||y - X_hat theta_hat||_2 / sqrt(n)
    #           (eq 3.14: sigma_hat^{(k+1)} <- ||y - X beta_hat^{(k)}||_2 / sqrt(n))
    #
    # Note: the paper uses sqrt(n), not sqrt(n - p_active), because the
    # scaled lasso formulation absorbs the degrees-of-freedom adjustment
    # into the sigma/2 term in the joint objective (3.12). This is the
    # oracle noise level sigma* = ||y - X beta*||_2 / sqrt(n).
    # ----------------------------------------------------------------
    resid_outer  <- y - X_hat %*% theta_hat
    sigma_new    <- sqrt(sum(resid_outer^2) / n)

    sigma_path[outer_iter] <- sigma_new

    if (verbose) {
      p_active <- sum(abs(theta_hat) > 1e-10)
      cat(sprintf("  sigma updated: %.6f -> %.6f | active theta: %d/%d\n",
                  sigma_hat, sigma_new, p_active, h))
    }

    # ----------------------------------------------------------------
    # Step (d): Convergence check on sigma
    #           Stop when relative change in sigma is below tol_outer
    # ----------------------------------------------------------------
    rel_change_sigma <- abs(sigma_new - sigma_hat) / (abs(sigma_hat) + 1e-10)

    sigma_hat <- sigma_new

    if (rel_change_sigma < tol_outer) {
      if (verbose) cat(sprintf("\nConverged at outer iter %d (rel delta sigma = %.2e)\n",
                                outer_iter, rel_change_sigma))
      sigma_path <- sigma_path[seq_len(outer_iter)]
      inner_results <- inner_results[seq_len(outer_iter)]

      # Clean up tiny theta values (numerical zeros)
      theta_hat[abs(theta_hat) < 1e-10] <- 0.0

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

  # If we reach here, max_outer iterations exhausted
  warning(sprintf("fit_scaled_group_lasso: max_outer (%d) reached without convergence. ",
                  max_outer),
          sprintf("Final rel change in sigma = %.2e. Consider increasing max_outer.",
                  rel_change_sigma))

  theta_hat[abs(theta_hat) < 1e-10] <- 0.0
  n_active <- sum(abs(theta_hat) > 1e-10)
  active_blocks <- which(sapply(seq_len(B), function(b) {
    idx_b <- block_start[b]:block_end[b]
    any(abs(theta_hat[idx_b]) > 1e-10)
  }))

  return(list(
    theta_hat     = theta_hat,
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

mz_group_wald_mr <- function(y, X_hat, Z = NULL, B_hat = NULL,
                              cg_blocks,
                              lambda_init         = NULL,
                              lambda_score        = NULL,
                              alpha_sgl           = 0.05,
                              fdr_level           = 0.05,
                              cond_threshold      = 30,
                              # Override parameters: supply pre-computed estimates
                              # from fit_scaled_group_lasso() to bypass internal
                              # estimation and use the theoretically justified
                              # scaled group lasso estimator directly.
                              # When both are supplied, Steps 1 and 2 are skipped.
                              theta_init_override = NULL,   # h-vector from sgl_fit$theta_hat
                              sigma_hat_override  = NULL,   # scalar from sgl_fit$sigma_hat
                              verbose             = TRUE) {

  n <- nrow(X_hat)
  h <- ncol(X_hat)   # total number of methylation sites
  B <- length(cg_blocks)

  stopifnot(sum(cg_blocks) == h)
  stopifnot(length(y) == n)

  # Validate overrides if supplied
  if (!is.null(theta_init_override)) {
    stopifnot(length(theta_init_override) == h)
  }
  if (!is.null(sigma_hat_override)) {
    stopifnot(is.numeric(sigma_hat_override) && length(sigma_hat_override) == 1)
    stopifnot(sigma_hat_override > 0)
  }

  # Both overrides must be supplied together or not at all
  use_overrides <- !is.null(theta_init_override) && !is.null(sigma_hat_override)

  block_id    <- rep(seq_len(B), times = cg_blocks)
  block_end   <- cumsum(cg_blocks)
  block_start <- c(1, head(block_end, -1) + 1)

  if (verbose) {
    cat("=== Mitra & Zhang Group Wald Test (MR Version) ===\n")
    cat(sprintf("n = %d, h = %d methylation sites, B = %d blocks\n", n, h, B))
    if (!is.null(Z)) {
      cat(sprintf("q = %d SNPs (instruments)\n", ncol(Z)))
      cat("  Note: Wald test validity requires strong instruments (F >> 10)\n")
    }
    if (use_overrides) {
      cat("  Using pre-supplied theta_init and sigma_hat from fit_scaled_group_lasso()\n")
      cat("  (Steps 1 and 2 bypassed -- theoretically justified scaled group lasso estimates)\n")
    }
  }

  # --------------------------------------------------------------------
  # STEP 1: Initial estimator theta_hat_init
  #
  # If theta_init_override is supplied (from fit_scaled_group_lasso()),
  # use it directly -- this is the theoretically preferred path because
  # the scaled group lasso satisfies the working assumption (2.31) with
  # a proven rate (Theorem 7 of the paper).
  #
  # Otherwise, fall back to grpreg (requires grpreg package).
  # --------------------------------------------------------------------
  if (use_overrides) {

    # Use pre-computed estimate from fit_scaled_group_lasso() (eq 2.25)
    theta_init <- theta_init_override
    if (!is.null(colnames(X_hat))) names(theta_init) <- colnames(X_hat)

  } else {

    if (verbose) cat("\nStep 1: Fitting Stage 2 SGL estimator on X_hat...\n")

    group_vec <- block_id

    if (is.null(lambda_init)) {
      cv_fit_mr <- cv.grpreg(X = X_hat, y = y, group = group_vec,
                             penalty = "grLasso",
                             alpha   = alpha_sgl,
                             nfolds  = 5)
      lambda_init <- cv_fit_mr$lambda.min
      if (verbose) cat(sprintf("  CV-selected lambda_init = %.4f\n", lambda_init))
    }

    fit_init_mr <- grpreg(X = X_hat, y = y, group = group_vec,
                          penalty = "grLasso",
                          alpha   = alpha_sgl,
                          lambda  = lambda_init)

    theta_init <- as.numeric(coef(fit_init_mr))[-1]
    names(theta_init) <- colnames(X_hat)
  }

  # --------------------------------------------------------------------
  # STEP 2: Noise estimate sigma_hat
  #
  # If sigma_hat_override is supplied (from fit_scaled_group_lasso()),
  # use it directly -- this is the jointly optimized sigma_hat proven
  # consistent by Theorem 6 (eq 3.21), which is what the chi-squared
  # approximation in (2.24) requires.
  #
  # Otherwise, estimate from residuals with conservative df correction.
  # --------------------------------------------------------------------
  if (use_overrides) {

    # Use jointly estimated sigma from fit_scaled_group_lasso() (eq 3.14)
    sigma_hat <- sigma_hat_override
    if (verbose) {
      p_active <- sum(abs(theta_init) > 1e-10)
      cat(sprintf("  sigma_hat (override) = %.4f, active exposures = %d\n",
                  sigma_hat, p_active))
    }

  } else {

    resid_mr  <- y - X_hat %*% theta_init
    p_active  <- sum(abs(theta_init) > 1e-10)
    df_resid  <- max(n - p_active, 1)
    sigma_hat <- sqrt(sum(resid_mr^2) / df_resid)

    if (verbose) cat(sprintf("  sigma_hat = %.4f, active exposures = %d\n",
                              sigma_hat, p_active))
  }

  # MR-specific diagnostic: first-stage F-statistics
  if (!is.null(Z) && !is.null(B_hat)) {
    iv_strength <- compute_iv_strength(Z, X_hat, cg_blocks, verbose)
  } else {
    iv_strength <- NULL
  }

  # --------------------------------------------------------------------
  # STEP 3-5: Loop over blocks
  #
  # Identical logic to the standard case, but with X_hat replacing X.
  #
  # KEY POINT: In the MR setting, X_hat = Z B_hat lies in the column
  # space of Z. The score matrix Z_b (NOT to be confused with the SNP
  # matrix Z!) is constructed by regressing X_hat_b on X_hat_{-b}.
  # Under between-block independence of CpG blocks, X_hat_b and
  # X_hat_{-b} are approximately uncorrelated, so Z_b ≈ X_hat_b.
  # --------------------------------------------------------------------
  if (verbose) cat("\nStep 3: Computing block-level Wald statistics (MR)...\n")

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

    # ----------------------------------------------------------------
    # Guard: skip blocks with zero selected sites
    # This happens when Stage 1 IV selection selected no methylation
    # sites from this block. The Wald test is undefined for an empty
    # block -- leave all results as NA and continue.
    # ----------------------------------------------------------------
    if (cg_blocks[b] == 0) {
      if (verbose) cat("EMPTY (no sites selected) -- skipped\n")
      next
    }

    idx_b     <- block_start[b]:block_end[b]
    idx_neg_b <- setdiff(seq_len(h), idx_b)

    # X_hat partitioned into block b and complement
    X_hat_b     <- X_hat[, idx_b,     drop = FALSE]  # n x d_b
    X_hat_neg_b <- X_hat[, idx_neg_b, drop = FALSE]  # n x (h - d_b)

    # ----------------------------------------------------------------
    # STEP 3a-b: Score matrix and projection
    #
    # Z_b = X_hat_b - X_hat_{-b} Gamma_hat_b
    #
    # In MR, between-block independence of CpG blocks means
    # X_hat_b ⊥ X_hat_{-b} approximately, so Gamma_hat_b ≈ 0
    # and Z_b ≈ X_hat_b.
    #
    # Within-block collinearity of X_hat_b is handled by SVD truncation
    # (same as standard case -- high within-block CpG correlation leads
    # to near-singular Z_b^T Z_b, addressed by rank-r projection).
    # ----------------------------------------------------------------
    sm_mr <- construct_score_matrix(X_hat_b, X_hat_neg_b,
                                     lambda         = lambda_score,
                                     cond_threshold = cond_threshold)
    score_info[[b]] <- sm_mr

    P_b <- sm_mr$P_G
    k_b <- sm_mr$k_G   # effective degrees of freedom after SVD truncation

    # ----------------------------------------------------------------
    # STEP 3c: Partial residual in MR setting
    #
    # partial_resid_b = y - X_hat_{-b} %*% theta_init_{-b}
    #
    # This is y with the estimated causal contributions of all
    # OTHER blocks removed, leaving block b's causal signal + noise.
    #
    # Under H_0 (block b has no causal effect on y):
    #   partial_resid_b ≈ eps  (pure noise)
    # Under H_1:
    #   partial_resid_b ≈ X_hat_b theta*_b + eps  (signal + noise)
    # ----------------------------------------------------------------
    theta_neg_b_mr  <- theta_init[idx_neg_b]
    partial_resid_b <- y - X_hat_neg_b %*% theta_neg_b_mr

    # ----------------------------------------------------------------
    # STEP 3d-e: Projection and Wald statistic
    #
    # W_b = (1/sigma_hat^2) ||P_b * partial_resid_b||_2^2
    #
    # Under H_0: theta*_b = 0:
    #   W_b ~ chi^2(k_b)   asymptotically
    #
    # The MR exclusion restriction ensures that under H_0, the
    # instruments Z affect y only through X_hat, so partial_resid_b
    # under H_0 is indeed approximately eps (not correlated with Z).
    # ----------------------------------------------------------------
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

  # BH correction
  results$padj <- p.adjust(results$pval, method = "BH")

  if (verbose) {
    cat(sprintf("\nBH correction at FDR = %.2f:\n", fdr_level))
    cat(sprintf("  Blocks rejected: %d / %d\n",
                sum(results$padj < fdr_level, na.rm = TRUE), B))
  }

  return(list(
    results      = results,
    theta_init   = theta_init,
    sigma_hat    = sigma_hat,
    score_info   = score_info,
    iv_strength  = iv_strength
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

  # Block start/end indices
  block_end   <- cumsum(cg_blocks)
  block_start <- c(1L, head(block_end, -1L) + 1L)

  # Precompute X_hat^T X_hat and X_hat^T y  (used repeatedly)
  XtX <- crossprod(X_hat)          # h x h
  Xty <- crossprod(X_hat, y)       # h x 1

  # Lipschitz constant for gradient of smooth part:
  # L = lambda_max(X^T X) / n
  # Use power iteration for efficiency (avoid full eigendecomposition)
  step_size <- 1.0 / (max(eigen(XtX, symmetric = TRUE,
                                only.values = TRUE)$values) / n + 1e-10)

  # Initialize theta
  if (is.null(theta_init)) {
    theta <- rep(0.0, h)
  } else {
    theta <- theta_init
  }

  objective_values <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {

    # ------------------------------------------------------------------
    # Gradient step:
    #   gradient = -(1/n) X_hat^T (y - X_hat theta)
    #            = -(1/n)(Xty - XtX theta)
    # Gradient descent update: theta_grad = theta - step * gradient
    # ------------------------------------------------------------------
    grad    <- -(Xty - XtX %*% theta) / n   # h x 1
    theta_g <- theta - step_size * grad      # h x 1 gradient step

    # ------------------------------------------------------------------
    # Proximal step: block-wise soft thresholding
    # For each block b:
    #   v_b = theta_g[block b]
    #   theta_new[block b] = v_b * max(0, 1 - step*omega_b / ||v_b||_2)
    #
    # omega is an h-vector where omega[j] = omega_b for j in block b.
    # We use omega[block_start[b]] to get the block-level weight.
    # ------------------------------------------------------------------
    theta_new <- theta_g  # copy, then apply proximal in place

    for (b in seq_len(B)) {
      idx_b   <- block_start[b]:block_end[b]
      v_b     <- theta_g[idx_b]
      norm_vb <- sqrt(sum(v_b^2))
      omega_b <- omega[block_start[b]]  # block-level weight

      threshold_b <- step_size * omega_b

      if (norm_vb <= threshold_b) {
        # Entire block shrinks to zero (group-level sparsity)
        theta_new[idx_b] <- 0.0
      } else {
        # Block soft-thresholding: shrink toward zero by factor
        theta_new[idx_b] <- v_b * (1.0 - threshold_b / norm_vb)
      }
    }

    # ------------------------------------------------------------------
    # Compute objective for convergence monitoring:
    #   f(theta) = (1/2n)||y - X_hat theta||^2 + sum_b omega_b ||theta_b||_2
    # ------------------------------------------------------------------
    resid_iter  <- y - X_hat %*% theta_new
    smooth_part <- sum(resid_iter^2) / (2.0 * n)
    penalty_part <- sum(sapply(seq_len(B), function(b) {
      idx_b   <- block_start[b]:block_end[b]
      omega_b <- omega[block_start[b]]
      omega_b * sqrt(sum(theta_new[idx_b]^2))
    }))
    objective_values[iter] <- smooth_part + penalty_part

    # ------------------------------------------------------------------
    # Convergence check: relative change in theta
    # ------------------------------------------------------------------
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

# ============================================================================
# REGION WALD TEST FUNCTION (Framework: Mitra & Zhang (2016))
# ============================================================================

# mz_group_wald_test()
#
# Region-level group Wald test for G2G Mendelian Randomization
# Framework: Mitra & Zhang (2016), adapted to the two-stage MR setting
#
# DESCRIPTION:
#   Performs the Mitra & Zhang group Wald test for causal effects of
#   methylation blocks on an outcome, using instrumented exposures from
#   a two-stage IV pipeline. Tests H_0: theta_b = 0 for each CpG block b,
#   where theta_b is the vector of causal effects for block b.
#
#   The test statistic for block b is:
#     W_b = ||P_b r_b||^2 / sigma_hat^2  ~  chi^2(d_b)  under H_0
#   where:
#     P_b   = X_hat_b (X_hat_b^T X_hat_b)^dagger X_hat_b^T  (score matrix)
#           = P_Z X_b (X_b^T P_Z X_b)^dagger X_b^T P_Z      (2SLS annihilator)
#     r_b   = y - X_hat_{-b} theta_{-b}                      (partial residual)
#     d_b   = rank(P_b) after SVD truncation                  (degrees of freedom)
#
#   P_b has rank d_b = cg_blocks[b] by construction because
#   X_hat_b = Z_b B_ols already maps instrument space to methylation space.
#   This satisfies the rank condition of Mitra & Zhang Theorem 2.
#
# INITIAL ESTIMATOR:
#   "sgl" (default, recommended):
#     Scaled sparse group lasso via fit_scaled_group_lasso().
#     Theoretically justified by MZ Theorem 6.
#     lambda_scale is selected adaptively based on combined_signal:
#       High noise (combined_signal < lambda_threshold): lambda_noisy
#       Detectable signal (combined_signal >= lambda_threshold): lambda_signal
#     Rationale: in the high-noise regime, greater penalization is needed
#     to prevent the estimator from assigning spurious nonzero coefficients
#     to null CpG sites, which inflates Type I error via the partial residual.
#
#   "en" (experimental):
#     Tuned elastic net via cv.glmnet() with alpha grid search [0.1, 0.9].
#     Better point estimates when signal is present. Collapses to zero
#     vector in ~18-26% of applications at effect size 0.10. en_fallback
#     controls behaviour on collapse. Not covered by MZ Theorem 6.
#
# COMBINED SIGNAL:
#   combined_signal = mean(first_stage_R^2) * second_stage_R^2
#   Computed internally from X_hat and y. Measures IV-methylation-outcome
#   signal strength. Used only for adaptive lambda selection (SGL).
#
# CORRECTIONS:
#   Applied simultaneously across all B block p-values.
#   "none"       : raw p-values vs fdr_level
#   "bh"         : Benjamini-Hochberg FDR control (recommended for
#                  multi-gene discovery applications, e.g. B = 30 genes)
#   "bonferroni" : FWER control, threshold = fdr_level / B
#   "all"        : all three correction columns returned
#
# REQUIREMENTS:
#   fit_scaled_group_lasso()  [scaled_group_lasso_stage2.R]
#   mz_group_wald_mr()        [mitra_zhang_wald_test.R]
#   glmnet                    [CRAN, required only when estimator = "en"]
#
# @param y               numeric vector length n. Outcome.
# @param X_hat           numeric matrix n x h. Instrumented exposures.
#                        Columns must correspond exactly to cg_blocks order.
# @param cg_blocks       integer vector length B. Selected CpG sites per block.
#                        Entries of 0 are allowed (empty blocks → NA output).
# @param block_names     character vector length B. Labels for each block.
#                        Defaults to "Block1", "Block2", ...
# @param null_blocks     logical vector length B. TRUE = null block (for
#                        downstream power/Type I error calculation). Optional.
#                        If provided, a null_block column is added to output.
# @param estimator       "sgl" (default) or "en".
# @param corrections     "none", "bh", "bonferroni", or "all". Default "bh".
# @param lambda_scale    numeric scalar or "adaptive" (default).
#                        If numeric, used directly. If "adaptive", selected
#                        from lambda_noisy / lambda_signal based on
#                        combined_signal vs lambda_threshold.
#                        Ignored when estimator = "en".
# @param lambda_threshold numeric. Threshold on combined_signal for adaptive
#                        lambda selection. Default 0.10.
# @param lambda_noisy    numeric. lambda_scale when combined_signal <
#                        lambda_threshold (high noise regime). Default 2.0.
# @param lambda_signal   numeric. lambda_scale when combined_signal >=
#                        lambda_threshold (detectable signal regime). Default 1.0.
# @param en_fallback     "warn" (default), "error", or "sgl".
#                        Behaviour when EN collapses to zero vector:
#                        "warn"  — issue warning, return NA Wald stats
#                        "error" — stop with error
#                        "sgl"   — silently refit with SGL and use that instead
# @param fdr_level       numeric. Significance threshold. Default 0.05.
# @param cond_threshold  numeric. SVD condition number threshold for score
#                        matrix truncation in mz_group_wald_mr(). Default 30.
# @param verbose         logical. Print per-block diagnostics. Default FALSE.
#
# @return data frame with one row per block:
#   block_name        block label
#   n_selected        CpG sites selected for this block
#   wald_stat         W_b test statistic (NA for empty blocks)
#   df                degrees of freedom = rank(P_b) after SVD truncation
#   pval_raw          chi-squared p-value before any correction
#   combined_signal   IV signal strength (same value repeated for all rows)
#   lambda_scale_used lambda_scale applied (same for all rows; NA if estimator="en")
#   sigma_hat         sigma from initial estimator (same for all rows)
#   estimator_used    "sgl" or "en" (may differ from input if en_fallback="sgl")
#   en_collapsed      TRUE/FALSE if EN zeroed out (NA if estimator="sgl")
#   null_block        TRUE/FALSE (only present if null_blocks argument provided)
#   --- correction columns based on `corrections` argument ---
#   padj_none         raw p-value           (if corrections = "none" or "all")
#   reject_none       0/1 at fdr_level      (if corrections = "none" or "all")
#   padj_bh           BH-adjusted p-value   (if corrections = "bh"   or "all")
#   reject_bh         0/1 at fdr_level      (if corrections = "bh"   or "all")
#   padj_bonferroni   Bonferroni-adjusted   (if corrections = "bonferroni" or "all")
#   reject_bonferroni 0/1 at fdr_level      (if corrections = "bonferroni" or "all")
# =============================================================================

mz_group_wald_test <- function(y,
                                X_hat,
                                cg_blocks,
                                block_names      = NULL,
                                null_blocks      = NULL,
                                estimator        = "sgl",
                                corrections      = "bh",
                                lambda_scale     = "adaptive",
                                lambda_threshold = 0.10,
                                lambda_noisy     = 2.0,
                                lambda_signal    = 1.0,
                                en_fallback      = "warn",
                                fdr_level        = 0.05,
                                cond_threshold   = 30,
                                verbose          = FALSE) {

  # ---------------------------------------------------------------------------
  # INPUT VALIDATION
  # ---------------------------------------------------------------------------
  y   <- as.numeric(y)
  n   <- nrow(X_hat)
  h   <- ncol(X_hat)
  B   <- length(cg_blocks)

  if (length(y) != n)
    stop("mz_group_wald_test: length(y) must equal nrow(X_hat).")

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

  if (!is.null(null_blocks) && length(null_blocks) != B)
    stop("mz_group_wald_test: null_blocks must have the same length as cg_blocks.")

  # Resolve block names
  if (is.null(block_names)) {
    block_names <- paste0("Block", seq_len(B))
  } else if (length(block_names) != B) {
    stop("mz_group_wald_test: block_names must have the same length as cg_blocks.")
  }

  # Resolve which correction columns to produce
  if ("all" %in% corrections) {
    active_corrections <- c("none", "bh", "bonferroni")
  } else {
    active_corrections <- unique(corrections)
  }

  # ---------------------------------------------------------------------------
  # COMBINED SIGNAL
  #
  # combined_signal = mean(first_stage_R^2) * second_stage_R^2
  # Used to select lambda_scale in the adaptive rule (SGL only).
  # first_stage_R^2 per column is estimated by regressing each X_hat column
  # on y and taking R^2 — a proxy for how much of y variance is explained by
  # the instrumented exposure, which reflects IV strength * causal effect.
  #
  # Note: this is an approximation computed entirely from X_hat and y with
  # no reference to the true signal level. It is fully data-driven.
  # ---------------------------------------------------------------------------
  second_stage_r2  <- summary(lm(y ~ X_hat))$r.squared
  first_stage_r2   <- sapply(seq_len(h), function(j) {
    summary(lm(X_hat[, j] ~ y))$r.squared
  })
  combined_signal  <- mean(first_stage_r2) * second_stage_r2

  if (verbose) {
    cat(sprintf("  combined_signal = %.4f\n", combined_signal))
  }

  # ---------------------------------------------------------------------------
  # ADAPTIVE LAMBDA SELECTION (SGL only)
  # ---------------------------------------------------------------------------
  if (is.character(lambda_scale) && lambda_scale == "adaptive") {
    lambda_scale_used <- if (combined_signal < lambda_threshold)
      lambda_noisy else lambda_signal
    if (verbose) {
      cat(sprintf("  Adaptive lambda_scale = %.1f (%s SNR regime)\n",
                  lambda_scale_used,
                  if (combined_signal < lambda_threshold) "low" else "high"))
    }
  } else {
    lambda_scale_used <- as.numeric(lambda_scale)
    if (verbose) cat(sprintf("  Fixed lambda_scale = %.1f\n", lambda_scale_used))
  }

  # ---------------------------------------------------------------------------
  # INITIAL ESTIMATOR
  # ---------------------------------------------------------------------------
  estimator_used <- estimator
  en_collapsed   <- NA   # stays NA for SGL

  if (estimator == "sgl") {

    # --- Sparse group lasso ---
    sgl_fit <- tryCatch({
      fit_scaled_group_lasso(
        y            = y,
        X_hat        = X_hat,
        cg_blocks    = cg_blocks,
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
    sigma_hat  <- sgl_fit$sigma_hat

    if (verbose) {
      cat(sprintf("  SGL: sigma=%.4f | active=%d / %d\n",
                  sigma_hat, sgl_fit$n_active, h))
    }

  } else {

    # --- Elastic net ---
    if (!requireNamespace("glmnet", quietly = TRUE))
      stop("mz_group_wald_test: package 'glmnet' is required for estimator = 'en'.")

    en_fit <- tryCatch({
      alpha_grid <- seq(0.1, 0.9, by = 0.1)
      cv_mse     <- numeric(length(alpha_grid))
      for (i in seq_along(alpha_grid)) {
        cv_tmp    <- glmnet::cv.glmnet(X_hat, y, alpha = alpha_grid[i])
        cv_mse[i] <- min(cv_tmp$cvm)
      }
      best_alpha <- alpha_grid[which.min(cv_mse)]
      cv_best    <- glmnet::cv.glmnet(X_hat, y, alpha = best_alpha)
      theta_en   <- as.numeric(glmnet::coef(cv_best, s = "lambda.min")[-1])
      list(theta_hat = theta_en, best_alpha = best_alpha)
    }, error = function(e) {
      stop(sprintf("mz_group_wald_test: elastic net fitting failed: %s", e$message))
    })

    theta_en     <- en_fit$theta_hat
    en_collapsed <- sum(abs(theta_en)) < 1e-10

    if (en_collapsed) {
      if (en_fallback == "error") {
        stop("mz_group_wald_test: EN collapsed to zero vector (en_fallback = 'error').")
      } else if (en_fallback == "warn") {
        warning("mz_group_wald_test: EN collapsed to zero vector. Returning NA Wald stats.")
        # Build NA result and return early
        out <- .build_na_result(
          block_names, cg_blocks, null_blocks,
          combined_signal, NA_real_, "en", TRUE,
          active_corrections, fdr_level, B
        )
        return(out)
      } else {
        # en_fallback == "sgl": silently refit with SGL
        warning("mz_group_wald_test: EN collapsed to zero vector. Falling back to SGL.")
        sgl_fb <- tryCatch({
          fit_scaled_group_lasso(
            y            = y,
            X_hat        = X_hat,
            cg_blocks    = cg_blocks,
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
        sigma_hat      <- sgl_fb$sigma_hat
        estimator_used <- "sgl"
      }
    } else {
      # EN produced nonzero estimates — compute sigma from instrumented residuals
      # (same convention as SGL: ||y - X_hat theta||^2 / n)
      theta_init <- theta_en
      sigma_hat  <- sqrt(sum((y - X_hat %*% theta_en)^2) / n)
    }

    if (verbose) {
      cat(sprintf("  EN: alpha=%.1f | active=%d / %d | collapsed=%s | sigma=%.4f\n",
                  en_fit$best_alpha, sum(abs(theta_init) > 1e-10),
                  h, en_collapsed, sigma_hat))
    }
  }

  # ---------------------------------------------------------------------------
  # MZ GROUP WALD TEST
  # Delegates entirely to mz_group_wald_mr() — the MR-adapted implementation
  # that uses X_hat_b as the score matrix (not raw Z_b), ensuring rank(P_b)=d_b
  # as required by Mitra & Zhang Theorem 2.
  # ---------------------------------------------------------------------------
  wald_out <- tryCatch({
    mz_group_wald_mr(
      y                   = y,
      X_hat               = X_hat,
      cg_blocks           = cg_blocks,
      theta_init_override = theta_init,
      sigma_hat_override  = sigma_hat,
      fdr_level           = fdr_level,
      cond_threshold      = cond_threshold,
      verbose             = FALSE
    )
  }, error = function(e) {
    stop(sprintf("mz_group_wald_test: mz_group_wald_mr failed: %s", e$message))
  })

  wr        <- wald_out$results   # data.frame: block, block_size, wald_stat,
                                  #             df, pval, padj, cond_num, svd_trunc
  raw_pvals <- wr$pval

  if (verbose) {
    cat(sprintf("  Wald stats: [%s]\n",
                paste(round(wr$wald_stat, 2), collapse = ", ")))
    cat(sprintf("  Raw p-vals: [%s]\n",
                paste(round(raw_pvals, 4), collapse = ", ")))
  }

  # ---------------------------------------------------------------------------
  # BUILD OUTPUT DATA FRAME
  # ---------------------------------------------------------------------------
  out <- data.frame(
    block_name        = block_names,
    n_selected        = cg_blocks,
    wald_stat         = wr$wald_stat,
    df                = wr$df,
    pval_raw          = raw_pvals,
    combined_signal   = combined_signal,
    lambda_scale_used = if (estimator == "en") NA_real_ else lambda_scale_used,
    sigma_hat         = sigma_hat,
    estimator_used    = estimator_used,
    en_collapsed      = en_collapsed,
    stringsAsFactors  = FALSE
  )

  # Add null_block column if provided
  if (!is.null(null_blocks)) {
    out$null_block <- null_blocks
  }

  # ---------------------------------------------------------------------------
  # APPLY CORRECTIONS
  # ---------------------------------------------------------------------------
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
    padj_bonf            <- p.adjust(raw_pvals, method = "bonferroni")
    out$padj_bonferroni  <- padj_bonf
    out$reject_bonferroni <- as.integer(!is.na(padj_bonf) & padj_bonf < fdr_level)
  }

  return(out)
}

# =============================================================================
# INTERNAL HELPER: build NA result for EN collapse with en_fallback = "warn"
# Not exported — used only inside mz_group_wald_test()
# =============================================================================
.build_na_result <- function(block_names, cg_blocks, null_blocks,
                              combined_signal, sigma_hat,
                              estimator_used, en_collapsed,
                              active_corrections, fdr_level, B) {
  out <- data.frame(
    block_name        = block_names,
    n_selected        = cg_blocks,
    wald_stat         = NA_real_,
    df                = NA_integer_,
    pval_raw          = NA_real_,
    combined_signal   = combined_signal,
    lambda_scale_used = NA_real_,
    sigma_hat         = sigma_hat,
    estimator_used    = estimator_used,
    en_collapsed      = en_collapsed,
    stringsAsFactors  = FALSE
  )
  if (!is.null(null_blocks)) out$null_block <- null_blocks
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

# =============================================================================
# PART 2: SUMMARY-STATISTIC PIPELINE B FUNCTIONS
# =============================================================================

# ============================================================================
# PIPELINE B: SUMMARY-STATISTIC G2G-MR FUNCTIONS
# ============================================================================
# Companion to MR_Sim_functions_new_sim20.R (individual-level Pipeline A)
#
# Architecture:
#   Stage 1: Proximal gradient IV selection on (Sigma_ZZ, beta_X)
#            -- summary-stat reformulation of Pipeline A's Stage 1
#   Stage 2a: Penalized IVW regression of beta_Y on beta_X[S_hat, ]  (IVW)
#   Stage 2b: Same + unpenalized scalar Egger intercept               (Egger)
#
# Key identity: individual-level gradient  ZtZ %*% B - ZtX
#               = Sigma_ZZ %*% B - Sigma_ZZ %*% beta_X
#               = Sigma_ZZ %*% (B - beta_X)
# where Sigma_ZZ = Z'Z/n  (from LD reference panel)
#       beta_X   = Sigma_ZZ^{-1} * (Z'X/n)  (marginal SNP-methylation effects)
# So Sigma_ZZ %*% beta_X recovers Z'X/n exactly.
# ============================================================================


# ============================================================================
# STAGE 1: SUMMARY-STAT SETUP (analog of setup_optimization)
# ============================================================================

#' Setup summary-stat optimization environment
#'
#' @param beta_X   p x h matrix of marginal SNP-methylation effect estimates
#' @param Sigma_ZZ p x p LD matrix from reference panel (Z'Z/n scale)
#' @param ld_blocks Vector of SNP LD block sizes
#' @param cg_blocks Vector of methylation CpG block sizes
#' @param lambda    Penalty parameter (> 0). See scaling note below.
#' @param a1        L1 weight
#' @param a2        LD-group L2 weight
#' @param n         Sample size used for lambda scaling (default NULL = no scaling).
#'                  If provided, lambda is internally divided by n to match the
#'                  summary-stat objective scale. See details.
#'
#' @details
#' LAMBDA SCALING:
#' The individual-level objective is (1/2)||X - ZB||_F^2, with gradient ZtZ*B - ZtX.
#' The summary-stat objective is (1/2) tr[(B-beta_X)' Sigma_ZZ (B-beta_X)],
#' with gradient Sigma_ZZ*(B - beta_X), where Sigma_ZZ = Z'Z/n.
#' The summary-stat gradient is therefore 1/n times the individual-level gradient.
#' Consequently, a lambda that works for the individual-level pipeline must be
#' divided by n to have the same effective regularization strength here.
#' Pass n to have this scaling applied automatically, or pre-scale lambda yourself.
#'
#' Recovers Sigma_ZX = Sigma_ZZ %*% beta_X, which is the summary-stat
#' equivalent of Z'X/n. The gradient then becomes:
#'   grad = Sigma_ZZ %*% B - Sigma_ZX = Sigma_ZZ %*% (B - beta_X)
#' which is structurally identical to the individual-level gradient,
#' allowing the same proximal operators to apply unchanged.
#'
#' In simulation: Sigma_ZZ = Z'Z/n (exact), beta_X = Z'X/n solved per-column.
#' In real data:  Sigma_ZZ from reference panel (e.g., 1000 Genomes),
#'                beta_X   from mQTL summary statistics.
#'                Reference panel mismatch robustness is noted as future work.

setup_optimization_sumstat <- function(beta_X, Sigma_ZZ, ld_blocks, cg_blocks,
                                       lambda, a1, a2, n = NULL,
                                       Sigma_ZX = NULL) {

  # --- Input validation ---
  if (lambda <= 0) stop("lambda must be positive")
  if (nrow(beta_X) != nrow(Sigma_ZZ)) stop("beta_X rows must match Sigma_ZZ dimension")
  if (sum(ld_blocks) != nrow(beta_X)) stop("Sum of ld_blocks must equal nrow(beta_X)")
  if (sum(cg_blocks) != ncol(beta_X)) stop("Sum of cg_blocks must equal ncol(beta_X)")

  # --- Lambda scaling ---
  # Optional: if n is provided, lambda is divided by n to convert from
  # individual-level scale to summary-stat scale (Sigma_ZZ = Z'Z/n).
  # In normal use lambda is already pre-scaled so n is left NULL.
  if (!is.null(n)) {
    lambda <- lambda / n
  }

  p <- nrow(beta_X)
  h <- ncol(beta_X)
  L <- length(ld_blocks)
  C <- length(cg_blocks)

  # --- Group indicators (identical structure to setup_optimization) ---
  ld_group_ind <- rep(1:L, times = ld_blocks)
  cg_group_ind <- rep(1:C, times = cg_blocks)

  ld_group_indices <- vector("list", L)
  for (l in 1:L) ld_group_indices[[l]] <- which(ld_group_ind == l)

  cg_group_indices <- vector("list", C)
  for (c in 1:C) cg_group_indices[[c]] <- which(cg_group_ind == c)

  # --- Key precomputation ---
  # Sigma_ZX = Z'X/n — the correct gradient target in summary-stat space.
  # If provided directly (from compute_summary_stats), use it.
  # If not provided, fall back to Sigma_ZZ %*% beta_X ONLY if beta_X contains
  # the multivariate (joint) regression coefficients — NOT univariate marginals.
  # With univariate marginal betas, Sigma_ZZ %*% beta_X != Z'X/n due to LD.
  if (!is.null(Sigma_ZX)) {
    # Use directly — this is the preferred path when called via run_pipeline_B
    # which passes sumstats$Sigma_ZX from compute_summary_stats
  } else {
    # Fallback: only correct when beta_X = multivariate (joint) coefficients
    warning("Sigma_ZX not provided; computing as Sigma_ZZ %*% beta_X. ",
            "This is only correct if beta_X contains joint (not marginal) coefficients.")
    Sigma_ZX <- Sigma_ZZ %*% beta_X
  }

  # Lipschitz constant for step size (same formula as individual-level)
  L_smooth <- max(eigen(Sigma_ZZ, only.values = TRUE)$values)

  cat("Summary-stat setup complete: p =", p, ", h =", h, ", L =", L, ", C =", C, "\n")
  cat("SNP block sizes:", ld_blocks, "\n")
  cat("CpG block sizes:", cg_blocks, "\n")
  cat("Lipschitz constant (max eigenvalue Sigma_ZZ):", round(L_smooth, 4), "\n")

  list(
    B            = matrix(0, nrow = p, ncol = h),  # Initialization
    beta_X       = beta_X,
    Sigma_ZZ     = Sigma_ZZ,
    Sigma_ZX     = Sigma_ZX,          # Z'X/n analog — used in gradient
    ld_blocks    = ld_blocks,
    cg_blocks    = cg_blocks,
    ld_group_ind = ld_group_ind,
    cg_group_ind = cg_group_ind,
    ld_group_indices = ld_group_indices,
    cg_group_indices = cg_group_indices,
    lambda       = lambda,
    a1           = a1,
    a2           = a2,
    L            = L,
    C            = C,
    h            = h,
    p            = p,
    L_smooth     = L_smooth
  )
}


# ============================================================================
# STAGE 1: SUMMARY-STAT PROXIMAL GRADIENT (analog of optimize_proximal_gradient_dykstra)
# ============================================================================

#' Proximal gradient IV selection on summary statistics
#'
#' Structurally identical to optimize_proximal_gradient_dykstra(), except:
#'   - gradient uses Sigma_ZZ %*% (B - beta_X) instead of ZtZ %*% B - ZtX
#'   - objective uses summary-stat reconstruction loss (see below)
#'   - no individual-level Z, X matrices required
#'
#' Objective:
#'   (1/2) tr[(B - beta_X)' Sigma_ZZ (B - beta_X)] + penalties(B)
#'
#' This is the population-level analog of (1/2n)||X - ZB||_F^2,
#' recoverable from (Sigma_ZZ, beta_X) alone.

optimize_proximal_gradient_sumstat <- function(setup, max_iter = 5000, tol = 1e-6,
                                               step_size = NULL, verbose = TRUE,
                                               n_inner = 3, decay_schedule = FALSE) {

  # Extract setup components
  B            <- setup$B
  Sigma_ZZ     <- setup$Sigma_ZZ
  Sigma_ZX     <- setup$Sigma_ZX     # = Sigma_ZZ %*% beta_X
  ld_group_indices <- setup$ld_group_indices
  cg_group_indices <- setup$cg_group_indices
  ld_blocks    <- setup$ld_blocks
  cg_blocks    <- setup$cg_blocks
  lambda       <- setup$lambda
  a1           <- setup$a1
  a2           <- setup$a2
  L            <- setup$L
  h            <- setup$h
  p            <- setup$p

  if (is.null(step_size)) {
    step_size <- 0.9 / setup$L_smooth
  }

  if (verbose) {
    cat("Starting summary-stat proximal gradient\n")
    cat("  Step size:", round(step_size, 6), "\n")
    cat("  Weights: a1 =", a1, ", a2 =", a2, ", a3 =", round(1 - a1 - a2, 3), "\n")
  }

  objective_values <- numeric(max_iter)

  for (iter in 1:max_iter) {
    B_old <- B

    # --- Gradient step (summary-stat version) ---
    # Individual-level: grad = ZtZ %*% B - ZtX
    # Summary-stat:     grad = Sigma_ZZ %*% B - Sigma_ZX
    #                        = Sigma_ZZ %*% (B - beta_X)   [equivalent]
    grad <- Sigma_ZZ %*% B - Sigma_ZX
    V    <- B - step_size * grad

    # --- Proximal step (identical Dykstra structure as Pipeline A) ---
    # Re-use cyclic_dykstra_proximal() directly — it only needs setup,
    # step_size, and V, none of which depend on individual-level data.
    B <- cyclic_dykstra_proximal(V, step_size, setup,
                                 n_cycles       = n_inner,
                                 decay_schedule = decay_schedule,
                                 outer_iter     = iter)

    # --- Objective (summary-stat reconstruction loss) ---
    # (1/2) tr[(B - beta_X)' Sigma_ZZ (B - beta_X)] + penalties
    D   <- B - setup$beta_X
    obj_smooth <- 0.5 * sum(diag(t(D) %*% Sigma_ZZ %*% D))
    # NOTE: penalty terms identical to individual-level; reuse compute_objective_sumstat below
    objective_values[iter] <- obj_smooth +
      compute_penalty_sumstat(B, setup)

    # --- Convergence ---
    frobenius_diff   <- norm(B - B_old, type = "F")
    relative_change  <- frobenius_diff / (norm(B, type = "F") + 1e-8)

    if (!is.na(frobenius_diff) &&
        (frobenius_diff < tol || relative_change < tol * 10)) {
      if (verbose) cat("Converged at iteration", iter,
                       "with ||Delta B||_F =", frobenius_diff, "\n")
      objective_values <- objective_values[1:iter]
      break
    }

    if (verbose && iter %% 50 == 0) {
      cat(sprintf("Iteration %d: Obj = %.6f, ||dB|| = %.2e\n",
                  iter, objective_values[iter], frobenius_diff))
    }
  }

  B_clean <- B
  B_clean[abs(B_clean) < 1e-8] <- 0
  n_active <- sum(abs(B_clean) > 1e-10)
  sparsity <- 1 - n_active / length(B_clean)

  if (verbose) {
    cat("Summary-stat proximal gradient complete:\n")
    cat("  Active coefficients:", n_active, "out of", length(B), "\n")
    cat("  Sparsity:", round(sparsity * 100, 2), "%\n")
    cat("  Final objective:", tail(objective_values, 1), "\n")
  }

  list(
    B                = B_clean,
    B_raw            = B,
    active_coeffs    = (abs(B_clean) > 1e-10),
    iterations       = iter,
    objective_values = objective_values,
    converged        = (iter < max_iter),
    sparsity         = sparsity,
    final_objective  = tail(objective_values, 1)
  )
}


#' Penalty-only objective for summary-stat setting
#' (smooth loss is handled separately above; penalty structure is unchanged)
compute_penalty_sumstat <- function(B, setup) {

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


# ============================================================================
# HELPER: EXTRACT SELECTED IVs FROM SUMMARY-STAT B ESTIMATE
# ============================================================================

#' Extract IV selection indicators from estimated B matrix
#' (identical logic to individual-level pipeline)
extract_iv_selection_sumstat <- function(B_est) {
  IVs_selected  <- apply(B_est, 1, function(row) l2norm(row) > 0)
  methyls_selected <- apply(B_est, 2, function(col) l2norm(col) > 0)
  list(IVs_selected = IVs_selected, methyls_selected = methyls_selected,
       n_IVs = sum(IVs_selected), n_methyls = sum(methyls_selected))
}


# ============================================================================
# STAGE 2a: PENALIZED IVW (no Egger)
# ============================================================================

#' Summary-stat Stage 2: penalized IVW regression
#'
#' Solves: min_theta (1/2) ||W^{1/2}(beta_Y - beta_X_S theta)||^2 + lambda2 ||theta||_1
#'
#' where W = diag(se_Y^{-2}) is the IVW weight matrix.
#' This is a standard weighted lasso, fit via glmnet.
#'
#' @param beta_Y     p* vector of SNP-outcome marginal effects (selected IVs)
#' @param beta_X_S   p* x h* matrix of SNP-methylation effects (selected IVs x selected methyls)
#' @param se_Y       p* vector of standard errors for beta_Y
#' @param lambda2    Lasso penalty on theta (scalar or "cv" for cross-validation)
#' @param nfolds     Number of CV folds if lambda2 = "cv"
#'
#' @return List with theta_hat, selected methyls, fitted glmnet object

stage2_ivw_sumstat <- function(beta_Y, beta_X_S, se_Y, lambda2 = "cv",
                               alpha_en = 1, lambda_rule = "1se", theta_0,
                               nfolds = 5) {
  # @param alpha_en    glmnet alpha: 1 = lasso, 0.5 = elastic net, 0 = ridge
  # @param lambda_rule "1se" (default, more stable in small p* regimes) or "min"
  #
  # NOTE on lambda_rule:
  # With p* = 20-40 selected IVs, CV folds have only ~4-8 rows each.
  # lambda.min overfits to noisy CV error curves in this regime.
  # lambda.1se selects a sparser, more stable solution and is strongly
  # preferred here. Use lambda_rule = "min" only if p* is large (>100).

  p_star <- length(beta_Y)
  h_star <- ncol(beta_X_S)

  if (length(se_Y) != p_star) stop("se_Y must have same length as beta_Y")
  if (nrow(beta_X_S) != p_star) stop("beta_X_S must have p* rows")

  # IVW weights: w_j = 1 / se_Y[j]^2
  w <- 1 / se_Y^2
  w <- w / mean(w)   # normalize for numerical stability

  # --- Early exit: glmnet requires >= 2 columns ---
  # If only 1 methylation site selected, penalized regression is undefined.
  # Fall back to unpenalized WLS — equivalent to a single-instrument IV estimate.
  if (h_star < 2) {
    fit_ols   <- lm(beta_Y ~ beta_X_S - 1, weights = w)
    coefs <- coef(fit_ols)
    names(coefs) <- colnames(beta_X_S)
    theta_hat <- map_coefs_to_vector_benchmark(coefs, theta_0)
    return(list(
      theta_hat       = theta_hat,
  
      methyls_nonzero = (abs(theta_hat) > 0),
      n_selected      = h_star,
      lambda_used     = 0,
      fit             = fit_ols,
      alpha_en        = alpha_en,
      lambda_rule     = "WLS_fallback_single_site",
      method          = "IVW"
    ))
  }

  # --- Early exit: lambda2=0 — unpenalized IVW (analog of OLS/2SLS) ---
  if (is.numeric(lambda2) && lambda2 == 0) {
    fit_ols   <- lm(beta_Y ~ beta_X_S - 1, weights = w)
    coefs <- coef(fit_ols)
    names(coefs) <- colnames(beta_X_S)
    theta_hat <- map_coefs_to_vector_benchmark(coefs, theta_0)
    return(list(
      theta_hat       = theta_hat,
      methyls_nonzero = (abs(theta_hat) > 0),
      n_selected      = sum(abs(theta_hat) > 0),
      lambda_used     = 0,
      fit             = fit_ols,
      alpha_en        = alpha_en,
      lambda_rule     = "unpenalized_IVW",
      method          = "IVW"
    ))
  }

  cv_rule <- if (lambda_rule == "1se") "lambda.1se" else "lambda.min"

  if (lambda2 == "cv") {

    # Guard: cv.glmnet needs at least nfolds rows; reduce folds if needed
    nfolds_use <- min(nfolds, p_star)
    if (nfolds_use < 3) {
      # Too few IVs for CV — fall back to OLS (unpenalized)
      warning(sprintf("p* = %d too small for CV (nfolds=%d). Falling back to OLS.",
                      p_star, nfolds))
      fit_ols   <- lm(beta_Y ~ beta_X_S - 1, weights = w)
      coefs <- coef(fit_ols)
      names(coefs) <- colnames(beta_X_S)
      theta_hat <- map_coefs_to_vector_benchmark(coefs, theta_0)
      return(list(theta_hat = theta_hat,
                  methyls_nonzero = (abs(theta_hat) > 0),
                  lambda_used = 0, fit = fit_ols,
                  alpha_en = alpha_en, lambda_rule = "OLS_fallback",
                  method = "IVW"))
    }

    fit <- cv.glmnet(
      x           = beta_X_S,
      y           = beta_Y,
      weights     = w,
      alpha       = alpha_en,
      intercept   = FALSE,
      nfolds      = nfolds_use,
      standardize = FALSE
    )
    theta_hat   <- as.numeric(coef(fit, s = cv_rule))[-1]
    names(theta_hat) <- colnames(beta_X_S)
    theta_hat <- map_coefs_to_vector_benchmark(theta_hat, theta_0)
    lambda_used  <- fit[[cv_rule]]

  } else {
    fit <- glmnet(
      x           = beta_X_S,
      y           = beta_Y,
      weights     = w,
      alpha       = alpha_en,
      intercept   = FALSE,
      lambda      = lambda2,
      standardize = FALSE
    )
    theta_hat   <- as.numeric(coef(fit))[-1]
    names(theta_hat) <- colnames(beta_X_S)
    theta_hat <- map_coefs_to_vector_benchmark(theta_hat, theta_0)
    lambda_used  <- lambda2
  }


  list(
    theta_hat       = theta_hat,
    methyls_nonzero = (abs(theta_hat) > 0),
    n_selected      = sum(abs(theta_hat) > 0),
    lambda_used     = lambda_used,
    fit             = fit,
    alpha_en        = alpha_en,
    lambda_rule     = lambda_rule,
    method          = "IVW"
  )
}

# ============================================================================
# STAGE 2a: SIMPLE IVW NO EGGER
# ============================================================================

stage2_g2g_ivw_sumstat <- function(beta_Y, beta_X_S, se_Y, theta_0) {
  w <- 1 / se_Y^2
  w <- w / mean(w)
  fit_ols   <- lm(beta_Y ~ beta_X_S - 1, weights = w)
  coefs <- coef(fit_ols)
  names(coefs) <- colnames(beta_X_S)
  theta_hat <- map_coefs_to_vector_benchmark(coefs, theta_0)
  list(theta_hat = theta_hat,
       methyls_nonzero = (abs(theta_hat) > 0))
}


# ============================================================================
# CONVENIENCE WRAPPER: full Pipeline B for one simulation replicate
# ============================================================================

#' Run full summary-stat Pipeline B (Stage 1 + Stage 2 IVW + Stage 2 Egger)
#'
#' @param beta_X   p x h marginal SNP-methylation effects
#' @param beta_Y   p vector of marginal SNP-outcome effects
#' @param se_Y     p vector of standard errors for beta_Y
#' @param Sigma_ZZ p x p LD matrix (Z'Z/n scale)
#' @param ld_blocks, cg_blocks  Block structure vectors
#' @param lambda1  Stage 1 penalty
#' @param a1, a2   Stage 1 penalty weights
#' @param lambda2  Stage 2 lasso penalty ("cv" or numeric)
#' @param ...      Additional args passed to optimize_proximal_gradient_sumstat
#'
#' @return List with stage1 result, ivw result, egger result, IV selection info

run_pipeline_B <- function(beta_X, beta_Y, se_Y, Sigma_ZZ,
                           ld_blocks, cg_blocks,
                           lambda1, a1, a2,
                           n = NULL,
                           Sigma_ZX = NULL, theta_0,
                           lambda2 = "cv", ...) {

  # --- Stage 1: IV selection ---
  setup <- setup_optimization_sumstat(
    beta_X   = beta_X,
    Sigma_ZZ = Sigma_ZZ,
    Sigma_ZX = Sigma_ZX,   # Z'X/n passed directly to avoid LD-contamination bug
    ld_blocks = ld_blocks,
    cg_blocks = cg_blocks,
    lambda   = lambda1,
    a1       = a1,
    a2       = a2,
    n        = n
  )

  stage1 <- optimize_proximal_gradient_sumstat(setup, ...)

  iv_info <- extract_iv_selection_sumstat(stage1$B)
  IVs_sel     <- iv_info$IVs_selected
  methyls_sel <- iv_info$methyls_selected

  cat("Pipeline B Stage 1: selected", iv_info$n_IVs, "IVs and",
      iv_info$n_methyls, "methylation sites\n")

  # Guard: need at least 1 IV and 1 methylation site
  # n_methyls == 1 is handled gracefully by stage2_ivw_sumstat (WLS fallback)
  # n_methyls == 0 or n_IVs == 0: cannot proceed at all
  if (iv_info$n_IVs == 0 || iv_info$n_methyls == 0) {
    warning("Pipeline B: no IVs or methylation sites selected. Returning NAs.")
    return(list(stage1=stage1, ivw_en=NULL, ivw_lasso=NULL,
                ivw=NULL, egger=NULL, iv_info=iv_info,
                IVs_selected=IVs_sel, methyls_selected=methyls_sel))
  }

  # Subset summary stats to selected IVs and methylation sites
  beta_Y_S  <- beta_Y[IVs_sel]
  se_Y_S    <- se_Y[IVs_sel]
  # beta_X_S uses marginal betas — required by the MVMR IVW estimating equation:
  #   beta_Y = beta_X * theta + eps
  # Both sides must be marginal (SNP-level) effects. Multivariate (joint)
  # betas would require a different estimating equation and are not compatible
  # with the standard IVW framework used here.
  beta_X_S  <- beta_X[IVs_sel, methyls_sel, drop = FALSE]

  # --- Data-adaptive lambda2 cap ---
  # When mean|beta_Y_S| < 0.20, the calibration grid upper bound can be reached
  # by a single atypical replicate, zeroing out all Stage 2 estimates.
  # Cap lambda2 at mean_bY * 0.15, validated as the MSE minimum plateau
  # (c=0.10-0.20) from a 30-replicate scale factor sweep at signal=0.10.
  # Threshold 0.20 cleanly separates signal=0.10 (mean_bY 0.11-0.25) from
  # signal=0.25+ (mean_bY 0.76+) with no overlap — fully observable in real
  # data without oracle knowledge of the true signal level.
  mean_bY_S <- mean(abs(beta_Y_S))
  if (mean_bY_S < 0.20 && lambda2 > mean_bY_S * 0.15) {
    lambda2_raw <- lambda2
    lambda2     <- mean_bY_S * 0.15
    cat(sprintf("  Pipeline B: adaptive lambda2 cap: %.5f -> %.5f (mean_bY=%.4f)\n",
                lambda2_raw, lambda2, mean_bY_S))
  }

  # --- Stage 2: three variants mirroring Pipeline A Stage 2 ---
  # Pipeline A runs: OLS (unpenalized), lasso, elastic net on instrumented X.
  # Pipeline B mirrors these in summary-stat IVW space:
  #   ivw     : unpenalized IVW (lambda2=0) — analog of OLS/2SLS
  #   ivw_lasso: lasso IVW (alpha=1)        — analog of lasso
  #   ivw_en  : elastic net IVW (alpha=0.5) — analog of elastic net
  # All use the same fixed lambda2 for the penalized variants.

  # 2a: Unpenalized IVW — analog of OLS Stage 2
  # lambda2=0 gives the standard IVW estimator with no variable selection.
  # Equivalent to mr_mvivw but restricted to Pipeline B's selected IV set.
  ivw <- stage2_ivw_sumstat(
    beta_Y = beta_Y_S, beta_X_S = beta_X_S, se_Y = se_Y_S, theta_0 = theta_0,
    lambda2 = 0, alpha_en = 0.5, lambda_rule = "min"
  )

  # 2b: Lasso IVW — analog of lasso Stage 2
  ivw_lasso <- stage2_ivw_sumstat(
    beta_Y = beta_Y_S, beta_X_S = beta_X_S, se_Y = se_Y_S,theta_0 = theta_0,
    lambda2 = lambda2, alpha_en = 1, lambda_rule = "min"
  )

  # 2c: Elastic net IVW — analog of elastic net Stage 2 (PRIMARY)
  # Handles correlated beta_X_S columns better than pure lasso.
  ivw_en <- stage2_ivw_sumstat(
    beta_Y = beta_Y_S, beta_X_S = beta_X_S, se_Y = se_Y_S,theta_0 = theta_0,
    lambda2 = lambda2, alpha_en = 0.5, lambda_rule = "min"
  )

  # --- Stage 2 Egger (unpenalized intercept + EN on theta) ---
  egger_result <- stage2_egger_sumstat(
    beta_Y   = beta_Y_S,
    beta_X_S = beta_X_S,
    se_Y     = se_Y_S,
    lambda2  = lambda2
  )

  list(
    stage1           = stage1,
    iv_info          = iv_info,
    ivw              = ivw,          # unpenalized IVW (OLS analog)
    ivw_lasso        = ivw_lasso,    # lasso IVW
    ivw_en           = ivw_en,       # elastic net IVW (primary)
    egger            = egger_result,
    IVs_selected     = IVs_sel,
    methyls_selected = methyls_sel
  )
}



# =============================================================================
# 6. fit_scaled_group_lasso_withU()
# =============================================================================

fit_scaled_group_lasso_withU <- function(y, X_hat, U, cg_blocks,
                                          lambda_scale = 1.0, sigma_init = NULL,
                                          max_outer = 50, max_inner = 2000,
                                          tol_outer = 1e-4, tol_inner = 1e-6,
                                          warm_start = TRUE, verbose = FALSE) {
  y <- as.numeric(y)
  n <- nrow(X_hat); h <- ncol(X_hat); q <- ncol(U); B <- length(cg_blocks)
  stopifnot(sum(cg_blocks)==h, length(y)==n, nrow(U)==n, lambda_scale>0)

  X_aug         <- cbind(X_hat, U)
  h_aug         <- h + q
  cg_blocks_aug <- c(cg_blocks, rep(1L, q))
  block_end_aug   <- cumsum(cg_blocks_aug)
  block_start_aug <- c(1L, head(block_end_aug, -1L) + 1L)

  omega_base <- rep(0.0, h_aug)
  for (b in seq_len(B)) {
    idx <- block_start_aug[b]:block_end_aug[b]
    omega_base[idx] <- sqrt(cg_blocks[b]/n) + sqrt(2*log(B)/n)
  }

  sigma_hat  <- if (is.null(sigma_init)) sd(y) else sigma_init
  theta_hat  <- rep(0.0, h_aug)
  sigma_path <- numeric(max_outer)

  for (outer_iter in seq_len(max_outer)) {
    inner_res <- solve_group_lasso_proxgrad(
      y=y, X_hat=X_aug, omega=lambda_scale*sigma_hat*omega_base,
      cg_blocks=cg_blocks_aug, max_iter=max_inner, tol=tol_inner,
      theta_init=if(warm_start) theta_hat else NULL, verbose=FALSE
    )
    theta_hat <- inner_res$theta
    sigma_new <- sqrt(sum((y - X_aug %*% theta_hat)^2) / n)
    sigma_path[outer_iter] <- sigma_new
    rel_change <- abs(sigma_new - sigma_hat) / (abs(sigma_hat) + 1e-10)
    sigma_hat  <- sigma_new
    if (rel_change < tol_outer) {
      sigma_path <- sigma_path[seq_len(outer_iter)]; break
    }
  }

  theta_meth <- theta_hat[seq_len(h)]
  theta_conf <- theta_hat[(h+1):h_aug]
  n_active   <- sum(abs(theta_meth) > 1e-10)

  if (verbose)
    cat(sprintf("Final: sigma=%.4f | active meth=%d/%d | conf nz=%d/%d\n",
                sigma_hat, n_active, h, sum(abs(theta_conf)>1e-10), q))

  list(theta_hat=theta_hat, theta_meth=theta_meth, theta_conf=theta_conf,
       sigma_hat=sigma_hat, n_active=n_active, sigma_path=sigma_path)
}


# =============================================================================
# HELPER: Construct score matrix Z_G for a given inference block
#
# Following equations (2.42)-(2.43) of Mitra & Zhang:
#   Gamma_hat = argmin_{Gamma} { (1/2n)||X_G - X_{-G} Gamma||_F^2 + R(Gamma) }
#   Z_G = X_G - X_{-G} %*% Gamma_hat
#   P_G = Z_G (Z_G^T Z_G)^{-1} Z_G^T
#
# KEY SIMPLIFICATION for your setting:
#   Because CpG blocks are approximately uncorrelated with each other,
#   X_{-G} and X_G are approximately orthogonal, so Gamma_hat ≈ 0
#   and Z_G ≈ X_G. This function computes the full version but if
#   between_block_cor is negligible, the penalty will drive Gamma_hat to 0.
#
# PENALTY CHOICE:
#   We use the Frobenius-norm group penalty (equation 2.47) rather than
#   the nuclear norm penalty (equation 2.44), which is appropriate when
#   |G| is not too large (see paper page 16-17). The Frobenius penalty
#   is computable with standard group lasso software.
#
# WITHIN-BLOCK COLLINEARITY:
#   If columns of Z_G are highly collinear, Z_G^T Z_G is near-singular.
#   We handle this via truncated SVD: keep only the top r singular vectors
#   where r is chosen by a condition number threshold.
#   This changes the degrees of freedom from |G| to r.
#
# @param X_G      n x d_b matrix of columns in the inference block
# @param X_neg_G  n x (p - d_b) matrix of all columns NOT in the block
# @param lambda   penalty parameter for the group lasso (equation 2.42)
#                 if NULL, uses a default data-driven choice
# @param cond_threshold  max acceptable condition number for Z_G^T Z_G
#                        before switching to truncated SVD (default: 30)
# @return list with Z_G (score matrix), P_G (projection), k_G (rank/df),
#         Gamma_hat (fitted nuisance regression coefficients),
#         svd_truncated (logical: was truncation applied?)
# =============================================================================
construct_score_matrix <- function(X_G, X_neg_G, lambda = NULL,
                                   cond_threshold = 30) {

  n    <- nrow(X_G)
  d_b  <- ncol(X_G)       # size of inference group |G|
  p_neg <- ncol(X_neg_G)  # number of nuisance variables

  # ------------------------------------------------------------------
  # Step 1: Fit penalized multivariate regression of X_G on X_{-G}
  #         This is equation (2.42) with Frobenius penalty (2.47)
  #         We do this column by column (each column of X_G separately)
  #         since the columns decouple under the Frobenius penalty.
  # ------------------------------------------------------------------
  Gamma_hat <- matrix(0, nrow = p_neg, ncol = d_b)

  if (p_neg > 0) {
    # Default lambda: scaled by n and block size as suggested by paper's
    # omega_{*,k} = sqrt(|G_k|/n) + sqrt(2 log M / n) scaling
    if (is.null(lambda)) {
      # Data-driven default: cross-validated lasso per column
      # (Since blocks are approx uncorrelated, this will select near-zero)
      lambda <- sqrt(log(p_neg) / n)
    }

    # Fit column by column (Frobenius penalty decouples across columns of X_G)
    for (j in 1:d_b) {
      # Each column of X_G regressed on X_{-G} via lasso
      # (Frobenius penalty = sum of squared Frobenius norms = sum of lasso penalties)
      fit_j <- glmnet(x = X_neg_G, y = X_G[, j],
                      alpha = 1,          # lasso
                      lambda = lambda,
                      intercept = FALSE,
                      standardize = FALSE)
      Gamma_hat[, j] <- as.numeric(coef(fit_j)[-1])  # drop intercept
    }
  }

  # ------------------------------------------------------------------
  # Step 2: Compute residual = score matrix Z_G
  #         Equation (2.43): Z_G = X_G - X_{-G} %*% Gamma_hat
  #         Under between-block independence: Z_G ≈ X_G
  # ------------------------------------------------------------------
  if (p_neg > 0) {
    Z_G <- X_G - X_neg_G %*% Gamma_hat
  } else {
    # No nuisance variables (testing the only block) -- Z_G = X_G exactly
    Z_G <- X_G
  }

  # ------------------------------------------------------------------
  # Step 3: Construct projection P_G = Z_G (Z_G^T Z_G)^{-1} Z_G^T
  #         This is equation (2.43)
  #         BUT: check condition number of Z_G^T Z_G first.
  #         High within-block collinearity => near-singular Z_G^T Z_G
  #         => use truncated SVD instead (rank-r projection)
  # ------------------------------------------------------------------

  # SVD of Z_G: Z_G = U D V^T
  # P_G = U_r U_r^T where U_r = first r left singular vectors
  svd_ZG <- svd(Z_G)
  singular_values <- svd_ZG$d

  # Condition number = max singular value / min singular value
  # (use only nonzero singular values)
  nonzero_sv <- singular_values[singular_values > 1e-10]
  cond_num   <- max(nonzero_sv) / min(nonzero_sv)

  svd_truncated <- FALSE

  if (cond_num > cond_threshold) {
    # Apply truncated SVD: keep singular vectors corresponding to
    # singular values above a threshold (relative to max)
    sv_threshold <- singular_values[1] / sqrt(cond_threshold)
    r <- sum(singular_values > sv_threshold)
    r <- max(r, 1)  # keep at least 1 direction
    svd_truncated <- TRUE
    cat(sprintf("  [score matrix] Block collinearity detected (cond = %.1f). ",
                cond_num))
    cat(sprintf("Truncating to r = %d of %d directions.\n", r, d_b))
  } else {
    r <- length(nonzero_sv)
  }

  # Rank-r projection: P_G = U_r %*% t(U_r)
  U_r <- svd_ZG$u[, 1:r, drop = FALSE]
  P_G <- U_r %*% t(U_r)   # n x n projection matrix, rank r

  return(list(
    Z_G           = Z_G,
    P_G           = P_G,
    U_r           = U_r,
    k_G           = r,               # degrees of freedom = rank of P_G
    Gamma_hat     = Gamma_hat,
    singular_values = singular_values,
    cond_num      = cond_num,
    svd_truncated = svd_truncated
  ))
}

# =============================================================================
# stage2_g2g_ivw_sumstat
# =============================================================================

# IVW helper (not in consolidated file)
stage2_g2g_ivw_sumstat <- function(beta_Y, beta_X_S, se_Y, theta_0) {
  w <- 1 / se_Y^2
  w <- w / mean(w)
  fit_ols   <- lm(beta_Y ~ beta_X_S - 1, weights = w)
  coefs <- coef(fit_ols)
  names(coefs) <- colnames(beta_X_S)
  theta_hat <- map_coefs_to_vector_benchmark(coefs, theta_0)
  list(theta_hat = theta_hat,
       methyls_nonzero = (abs(theta_hat) > 0))
}

# =============================================================================
# HELPER: Pipeline A Stage 2
# =============================================================================

compute_g2g_stage2 <- function(y, Z_sel, X_sel, theta_true,
                                cg_blocks_sel, lambda_threshold,
                                lambda_noisy, lambda_signal_val, n) {
  h_sel   <- ncol(X_sel)
  theta_0 <- rep(0, length(theta_true))
  names(theta_0) <- names(theta_true)

  X_hat <- matrix(0, nrow = n, ncol = h_sel)
  colnames(X_hat) <- colnames(X_sel)
  for (j in seq_len(h_sel)) {
    X_hat[, j] <- fitted(lm(X_sel[, j] ~ Z_sel))
  }

  first_stage_r2  <- sapply(seq_len(h_sel),
    function(j) summary(lm(X_hat[, j] ~ as.numeric(y)))$r.squared)
  second_stage_r2 <- summary(lm(y ~ X_hat))$r.squared
  combined_signal <- mean(first_stage_r2) * second_stage_r2
  lam_scale <- if (combined_signal < lambda_threshold) lambda_noisy else lambda_signal_val

  ols_fit   <- lm(y ~ X_hat)
  theta_ols <- map_coefs_to_vector(ols_fit, theta_0, zdata_name = "X_hat")
  mse_ols   <- calculate_theta_diff(theta_ols, theta_true)

  lasso_fit   <- cv.glmnet(X_hat, y, alpha = 1)
  theta_lasso <- map_coefs_to_vector(lasso_fit, theta_0)
  mse_lasso   <- calculate_theta_diff(theta_lasso, theta_true)

  alpha_grid  <- seq(0.1, 0.9, by = 0.1)
  cv_mse_vals <- sapply(alpha_grid, function(a)
    min(cv.glmnet(X_hat, y, alpha = a)$cvm))
  best_alpha  <- alpha_grid[which.min(cv_mse_vals)]
  en_fit      <- cv.glmnet(X_hat, y, alpha = best_alpha)
  theta_en    <- map_coefs_to_vector(en_fit, theta_0)
  mse_en      <- calculate_theta_diff(theta_en, theta_true)

  sgl_fit <- tryCatch({
    fit_scaled_group_lasso(
      y = as.numeric(y), X_hat = X_hat, cg_blocks = cg_blocks_sel,
      lambda_scale = lam_scale, sigma_init = sd(as.numeric(y)),
      max_outer = 50, max_inner = 2000, tol_outer = 1e-4, tol_inner = 1e-6,
      warm_start = TRUE, verbose = FALSE)
  }, error = function(e) NULL)

  if (!is.null(sgl_fit)) {
    theta_sgl <- rep(0, length(theta_true)); names(theta_sgl) <- names(theta_true)
    theta_sgl[colnames(X_hat)] <- sgl_fit$theta_hat
    mse_sgl <- calculate_theta_diff(theta_sgl, theta_true)
  } else {
    theta_sgl <- rep(NA_real_, length(theta_true)); mse_sgl <- NA_real_
  }

  list(mse_ols = mse_ols, mse_lasso = mse_lasso, mse_en = mse_en, mse_sgl = mse_sgl,
       theta_ols = theta_ols, theta_lasso = theta_lasso, theta_en = theta_en,
       theta_sgl = theta_sgl,
       X_hat = X_hat, cg_blocks_sel = cg_blocks_sel, lam_scale = lam_scale)
}


# =============================================================================
# Two arms, both returning the same column layout as mz_group_wald_test(), so the
# existing analysis code can read them interchangeably:
#
#   raw     S_g = X_g' (y - pi_-g)                        <- as in mz_group_wald_logistic
#   decorr  S_g = (X_g - X_-g Gamma)' (y - pi_-g)         <- W-metric residualised
#
# Both use the Fisher information I_g = Sg_mat' W Sg_mat with W = diag(pi(1-pi)),
# so neither needs a scalar sigma^2 and both are heteroskedasticity-correct.
#
# The "decorr" arm is the logistic analogue of construct_score_matrix(): it removes
# the linear part of X_-g from X_g before forming the score, which is what makes the
# chi^2(r_eff) reference valid when theta_hat_-g carries regularisation bias. The
# residualisation is PENALISED (lasso, as in construct_score_matrix) and swept over
# score_lambda_grid -- see the note in .logit_wald_core for why an unpenalised
# Gamma_hat destroys power.
#
#
# theta_init is obtained from the SAME scaled group lasso and the SAME adaptive lambda
# rule that mz_group_wald_test() uses internally, so the arms differ only in the score
# construction, not in the initial estimator.
#
# NOTE ON y: the logistic arms need the raw 0/1 outcome (pi_-g supplies the level).
# The Gaussian arm needs the demeaned outcome (it fits no intercept). Both are passed.
# =============================================================================

.logit_wald_core <- function(y01, X_hat, U, cg_blocks, block_names, null_blocks,
                             theta_init, gamma_init, decorr, score_lambda = 1.0,
                             cond_threshold = 30, fdr_level = 0.05) {
  h <- ncol(X_hat); G <- length(cg_blocks)
  be <- cumsum(cg_blocks); bs <- c(1L, head(be, -1L) + 1L)
  hasU <- !is.null(U) && ncol(U) > 0

  out <- data.frame(block_name = block_names, n_selected = cg_blocks,
                    wald_stat = NA_real_, df = NA_integer_, pval_raw = NA_real_,
                    cond_num = NA_real_, svd_trunc = NA,
                    stringsAsFactors = FALSE)

  for (g in seq_len(G)) {
    if (cg_blocks[g] == 0) next
    ig  <- bs[g]:be[g]
    ing <- setdiff(seq_len(h), ig)
    Xg  <- X_hat[, ig,  drop = FALSE]
    Xng <- X_hat[, ing, drop = FALSE]

    # null-model linear predictor: block g excluded, confounders retained
    eta <- as.numeric(Xng %*% theta_init[ing])
    if (hasU) eta <- eta + as.numeric(U %*% gamma_init)
    pi_ <- pmax(pmin(1 / (1 + exp(-eta)), 1 - 1e-10), 1e-10)
    w   <- pi_ * (1 - pi_)

    Sg_mat <- Xg
    if (decorr) {
      # PENALISED W-metric residualisation of Xg on [Xng, U]: the logistic analogue
      # of construct_score_matrix()'s Z_G = X_G - X_-G Gamma_hat.
      #
      # The penalty is essential, not a refinement. construct_score_matrix uses a
      # lasso at lambda = sqrt(log(p_neg)/n); an UNPENALISED (OLS) Gamma_hat removes
      # everything in span(X_-g) -- with ~21-31 nuisance columns at n=300 that strips
      # the block's real signal along with the bias. Measured with OLS: all three true
      # signal blocks undetected (p = 0.98 / 0.45 / 0.70) while the Gaussian and raw
      # arms each found one, and null W/df pinned at 0.59-0.83 across every rho.
      #
      # Fitting on the sqrt(w)-weighted matrices makes this the W-metric projection;
      # the residual Xg - Anui %*% Gam is then returned on the original scale.
      sw   <- sqrt(w)
      Anui <- if (hasU) cbind(Xng, U) else Xng
      A    <- Anui * sw
      Bm   <- Xg * sw
      lam  <- score_lambda * sqrt(log(max(ncol(Anui), 2)) / nrow(X_hat))
      Gam  <- matrix(0, ncol(Anui), ncol(Xg))
      for (j in seq_len(ncol(Bm))) {
        Gam[, j] <- tryCatch(
          as.numeric(coef(glmnet::glmnet(x = A, y = Bm[, j], alpha = 1, lambda = lam,
                                         intercept = FALSE, standardize = FALSE)))[-1],
          error = function(e) rep(0, ncol(Anui)))
      }
      Sg_mat <- Xg - Anui %*% Gam
    }

    r <- as.numeric(y01) - pi_
    S <- crossprod(Sg_mat, r)
    I <- crossprod(Sg_mat * w, Sg_mat)

    sv <- svd(I)
    if (!is.finite(sv$d[1]) || sv$d[1] <= .Machine$double.eps) next
    keep  <- sv$d > sv$d[1] / cond_threshold
    r_eff <- sum(keep)
    if (r_eff == 0) next
    nz <- sv$d[sv$d > .Machine$double.eps]
    out$cond_num[g]  <- if (length(nz)) max(nz) / min(nz) else NA_real_
    out$svd_trunc[g] <- r_eff < ncol(Sg_mat)

    Iinv <- sv$v[, keep, drop = FALSE] %*%
            diag(1 / sv$d[keep], nrow = r_eff, ncol = r_eff) %*%
            t(sv$u[, keep, drop = FALSE])
    Wst <- as.numeric(t(S) %*% Iinv %*% S)
    out$wald_stat[g] <- Wst
    out$df[g]        <- r_eff
    out$pval_raw[g]  <- pchisq(Wst, df = r_eff, lower.tail = FALSE)
  }

  if (!is.null(null_blocks)) out$null_block <- null_blocks
  padj_bh   <- p.adjust(out$pval_raw, method = "BH")
  padj_bonf <- p.adjust(out$pval_raw, method = "bonferroni")
  out$padj_none         <- out$pval_raw
  out$reject_none       <- as.integer(!is.na(out$pval_raw) & out$pval_raw < fdr_level)
  out$padj_bh           <- padj_bh
  out$reject_bh         <- as.integer(!is.na(padj_bh)   & padj_bh   < fdr_level)
  out$padj_bonferroni   <- padj_bonf
  out$reject_bonferroni <- as.integer(!is.na(padj_bonf) & padj_bonf < fdr_level)
  out
}

# Public entry point. Returns list(raw = <df>, decorr = <df>, raw_sweep = <df>), any
# of which may be NULL.
#
# score_lambda_grid: multipliers on the MZ base penalty sqrt(log(p_neg)/n) used for the
# decorrelating regression. Same style as the existing lambda_noisy_grid sweeps. The
# decorr result is returned stacked with a score_lambda column, one block-set per value.
# Pass NULL to skip the decorr arm entirely.
#
# raw_lambda_grid: SGL lambda_scale values for a lambda_noisy sweep of the RAW arm.
# score_lambda sweeps the score CONSTRUCTION and so applies only to decorr; the raw
# arm's only free knob is the lambda that produced theta_init, which the deployed run
# never varied (the adaptive rule fired "noisy" in every sim of every cell, so every
# stored raw row sits at lambda_noisy = 0.5). Supplying a grid here refits the scaled
# group lasso at each value and re-runs the raw score against it, returning the results
# stacked with a lambda_noisy column.
#
# The grid values BYPASS the adaptive rule -- each is used as lambda_scale directly.
# That matches how the Gaussian sweep in the cell scripts forces its candidates
# (lambda_threshold = Inf, so the "noisy" branch always fires). Include the deployed
# 0.5 in the grid and its rows must reproduce `raw` exactly; that is a free check.
#
# Default NULL: no refits, no cost, existing call sites unchanged.
logit_wald_arms <- function(y01, y_centered, X_hat, U = NULL, cg_blocks,
                            block_names, null_blocks,
                            lambda_threshold, lambda_noisy, lambda_signal,
                            score_lambda_grid = c(0.05, 0.1, 0.25, 0.5, 1.0),
                            raw_lambda_grid = NULL,
                            cond_threshold = 30, fdr_level = 0.05) {
  hasU <- !is.null(U) && ncol(U) > 0
  h <- ncol(X_hat)

  # --- centre X_hat (and U): REQUIRED, not cosmetic ---------------------------
  # The logistic arms fit no intercept, and compute_g2g_stage2 builds X_hat via
  # fitted(lm(X ~ Z)) with an intercept, so X_hat inherits mean(X) -- column means
  # reach 1.10 against a median SD of 0.61. Left uncentred, I_g = X_g' W X_g is
  # dominated by the mean direction, becomes near rank-1, and SVD truncation at
  # cond_threshold guts the df (observed: df 15 -> 4, null p = 6e-07 on a TRUE null).
  # Centring is the FWL equivalent of fitting an intercept and leaves the score's
  # slope content unchanged. y01 is deliberately NOT centred: pi_-g supplies the level.
  X_hat <- sweep(X_hat, 2, colMeans(X_hat))
  if (hasU) U <- sweep(U, 2, colMeans(U))

  # --- same adaptive lambda rule as mz_group_wald_test() ---
  first_r2 <- vapply(seq_len(h), function(j)
    summary(lm(X_hat[, j] ~ as.numeric(y_centered)))$r.squared, numeric(1))
  second_r2 <- if (hasU) summary(lm(y_centered ~ X_hat + U))$r.squared
               else      summary(lm(y_centered ~ X_hat))$r.squared
  combined <- mean(first_r2) * second_r2
  lam <- if (combined < lambda_threshold) lambda_noisy else lambda_signal

  # --- same initial estimator (scaled group lasso) ---
  # Factored out so the raw sweep can refit it at each candidate lambda. The call is
  # byte-identical to what the single-lambda version used; only lambda_scale varies.
  fit_init <- function(lambda_scale) {
    sgl <- tryCatch({
      if (hasU) fit_scaled_group_lasso_withU(y = as.numeric(y_centered), X_hat = X_hat, U = U,
                  cg_blocks = cg_blocks, lambda_scale = lambda_scale,
                  sigma_init = sd(as.numeric(y_centered)), max_outer = 50, max_inner = 2000,
                  tol_outer = 1e-4, tol_inner = 1e-6, warm_start = TRUE, verbose = FALSE)
      else      fit_scaled_group_lasso(y = as.numeric(y_centered), X_hat = X_hat,
                  cg_blocks = cg_blocks, lambda_scale = lambda_scale,
                  sigma_init = sd(as.numeric(y_centered)), max_outer = 50, max_inner = 2000,
                  tol_outer = 1e-4, tol_inner = 1e-6, warm_start = TRUE, verbose = FALSE)
    }, error = function(e) NULL)
    if (is.null(sgl)) return(NULL)
    list(theta = if (hasU) sgl$theta_meth else sgl$theta_hat,
         gamma = if (hasU) sgl$theta_conf else numeric(0))
  }

  go <- function(dec, sl, theta_init, gamma_init) tryCatch(
    .logit_wald_core(y01, X_hat, U, cg_blocks, block_names, null_blocks,
                     theta_init, gamma_init, decorr = dec, score_lambda = sl,
                     cond_threshold = cond_threshold, fdr_level = fdr_level),
    error = function(e) NULL)

  # --- deployed baseline: one SGL fit at the adaptive lambda, both arms off it ------
  init <- fit_init(lam)
  if (is.null(init)) return(list(raw = NULL, decorr = NULL, raw_sweep = NULL))

  raw <- go(FALSE, 1.0, init$theta, init$gamma)  # no residualisation -> score_lambda unused

  decorr <- NULL
  if (!is.null(score_lambda_grid) && length(score_lambda_grid)) {
    dec_list <- lapply(score_lambda_grid, function(sl) {
      d <- go(TRUE, sl, init$theta, init$gamma)
      if (!is.null(d)) d$score_lambda <- sl
      d
    })
    dec_list <- Filter(Negate(is.null), dec_list)
    if (length(dec_list)) decorr <- do.call(rbind, dec_list)
  }

  # --- raw-arm lambda_noisy sweep: one SGL refit per candidate ---------------------
  # The SGL refit is the cost here, not the score; each candidate is a fresh theta_init.
  raw_sweep <- NULL
  if (!is.null(raw_lambda_grid) && length(raw_lambda_grid)) {
    sw_list <- lapply(raw_lambda_grid, function(lv) {
      ini <- fit_init(lv)
      if (is.null(ini)) return(NULL)
      d <- go(FALSE, 1.0, ini$theta, ini$gamma)
      if (!is.null(d)) d$lambda_noisy <- lv
      d
    })
    sw_list <- Filter(Negate(is.null), sw_list)
    if (length(sw_list)) raw_sweep <- do.call(rbind, sw_list)
  }

  list(raw = raw, decorr = decorr, raw_sweep = raw_sweep, lambda_scale_used = lam,
       combined_signal = combined, score_lambda_grid = score_lambda_grid,
       raw_lambda_grid = raw_lambda_grid)
}

