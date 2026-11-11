simulate_SNP_matrix <- function(p, n, LD_var_matrix, maf) {
  # Create latent Gaussian variable matrix Z such that Z_i ~ MVN(0, LD_corr_matrix)
  mu <- numeric(p) # mean 0 vector
  Z <- MASS::mvrnorm(n, mu = mu, Sigma = LD_var_matrix)   
  
  # Find constants based on MAF
  prob_g0 <- dbinom(0, 2, maf)
  prob_g1 <- dbinom(1, 2, maf)
  prob_g2 <- dbinom(2, 2, maf)
  

  #P(X = 0) = P(Z_[ij] < c1j)
  c1 <- qnorm(prob_g0)
  #P(X = 2) = P(Z_[ij] > c2j)
  c2 <- qnorm(prob_g2, lower.tail = FALSE)
  
  # Create SNP matrix where variance = 1
  G <- Z
  G[Z <= c1] <- 0
  G[(Z > c1) & (Z <= c2)] <- 1
  G[Z > c2] <- 2
  
  matrices_list <- list(SNP = G, latent = Z)
  return(matrices_list)
}

simulate_confounders_matrix <- function(n, q) {
  
  # Generate independent standard normal confounders
  U <- matrix(rnorm(n * q, mean = 0, sd = 1), nrow = n, ncol = q)
  
  # Add meaningful column names
  colnames(U) <- paste0("U", 1:q)
  
  return(U)
}

assign_column_names <- function(prefix = "x", m) {
  column_names <- paste0(prefix, 1:m)
  return(column_names)
}

simulate_methyl_sites_blocked_withU <- function(cg_blocks, gene_id,
                                                 SNP_matrix, Coeff_matrix,
                                                 sigma_methyl,
                                                 Conf_matrix, Conf_coeff_matrix,
                                                 rho_within = 0.85, rho_decay = 0.3,
                                                 rho_block = 0.05, rho_gene = 0.01,
                                                 structure = "Simple") {
  if (length(gene_id) != length(cg_blocks))
    stop("gene_id must have same length as cg_blocks")
  if (!(structure %in% c("Simple","AR1")))
    stop("structure must be 'Simple' or 'AR1'")

  h <- sum(cg_blocks); n_blocks <- length(cg_blocks)
  site_block_id <- rep(seq_len(n_blocks), times = cg_blocks)
  site_gene_id  <- rep(gene_id,           times = cg_blocks)
  block_starts  <- c(1, cumsum(cg_blocks[-n_blocks]) + 1)

  if (structure == "Simple") {
    R <- diag(h)
  } else {
    R <- matrix(0, h, h)
    for (i in seq_len(h)) {
      for (j in i:h) {
        if (i == j) {
          R[i,j] <- 1
        } else if (site_block_id[i] == site_block_id[j]) {
          bk <- site_block_id[i]
          d  <- abs((i - block_starts[bk]) - (j - block_starts[bk]))
          R[i,j] <- rho_block + (rho_within - rho_block) * rho_decay^(d-1)
        } else if (site_gene_id[i] == site_gene_id[j]) {
          R[i,j] <- rho_block
        } else {
          R[i,j] <- rho_gene
        }
        R[j,i] <- R[i,j]
      }
    }
    ev <- eigen(R, symmetric=TRUE, only.values=TRUE)$values
    if (any(ev <= 0)) {
      message("Projecting to nearest PD matrix.")
      R <- as.matrix(Matrix::nearPD(R, corr=TRUE)$mat)
    }
  }

  Sigma <- sigma_methyl^2 * R
  X <- t(apply(SNP_matrix %*% Coeff_matrix + Conf_matrix %*% Conf_coeff_matrix,
               1, function(mu) MASS::mvrnorm(1, mu=mu, Sigma=Sigma)))
  colnames(X) <- colnames(Coeff_matrix)

  list(X=X, cg_blocks=cg_blocks, gene_id=gene_id,
       site_block_id=site_block_id, site_gene_id=site_gene_id, R=R, Sigma=Sigma)
}

simulate_fixed_confounder_effects <- function(cg_blocks, gene_id, 
                                                      signal_strength,
                                                      sparsity_within_gene = 0.4,
                                                      seed = 12345) {
  # Set seed for reproducibility
  set.seed(seed)
  
  h <- sum(cg_blocks)
  q_conf <- 3
  
  A <- matrix(0, nrow = q_conf, ncol = h)
  
  # Identify which columns belong to which gene
  site_gene_id <- rep(gene_id, times = cg_blocks)
  gene1_cols <- which(site_gene_id == 1)
  gene2_cols <- which(site_gene_id == 2)
  gene3_cols <- which(site_gene_id == 3)
  
   sparse_effect <- function(cols, strength, sparsity, confounder_id) {
    n_affected <- ceiling(length(cols) * sparsity)
    
    set.seed(seed + confounder_id * 100 + cols[1])
    affected <- sample(cols, n_affected)
    
    effects <- numeric(length(cols))
    indices_in_cols <- which(cols %in% affected)
    effects[indices_in_cols] <- strength
    return(effects)
  }
  
   # U1 → Gene 2, Gene 3 (not Gene 1)
  A[1, gene2_cols] <- sparse_effect(gene2_cols, signal_strength, sparsity_within_gene, confounder_id = 1)
  A[1, gene3_cols] <- sparse_effect(gene3_cols, signal_strength, sparsity_within_gene, confounder_id = 1)
  
  # U2 → Gene 1, Gene 3 (not Gene 2)
  A[2, gene1_cols] <- sparse_effect(gene1_cols, signal_strength, sparsity_within_gene, confounder_id = 2)
  A[2, gene3_cols] <- sparse_effect(gene3_cols, signal_strength, sparsity_within_gene, confounder_id = 2)
  
  # U3 → Gene 1, Gene 2 (not Gene 3)
  A[3, gene1_cols] <- sparse_effect(gene1_cols, signal_strength, sparsity_within_gene, confounder_id = 3)
  A[3, gene2_cols] <- sparse_effect(gene2_cols, signal_strength, sparsity_within_gene, confounder_id = 3)
  
  rownames(A) <- paste0("U", 1:q_conf)
  colnames(A) <- paste0("cg", 1:h)
  
  return(A)
}
simulate_mrg2g_example <- function(n              = 300,
                                   signal         = 0.25,
                                   noise          = 1,
                                   conf_strength  = 0.25,
                                   gamma_strength = 0.25,
                                   cov_structure  = "AR1",
                                   rho_within     = 0.25,
                                   rho_decay      = 0.85,
                                   rho_block      = 0.05,
                                   rho_gene       = 0.03,
                                   maf            = 0.24,
                                   ld_rho_within  = 0.6,
                                   ld_rho_between = 0.1,
                                   standardize    = TRUE,
                                   seed           = NULL) {

  ld_blocks_gene <- c(5, 10, 20, 10, 5, 5, 10, 20, 10, 5)
  ld_blocks      <- rep(ld_blocks_gene, 3)
  p_gene         <- sum(ld_blocks_gene)
  cg_blocks      <- c(15, 15, 10, 10)
  methyl_gene_id <- c(1, 2, 3, 3)
  h              <- sum(cg_blocks)
  p              <- sum(ld_blocks)

  LD_g <- matrix(ld_rho_between, p_gene, p_gene)
  rs <- 1
  for (sz in ld_blocks_gene) {
    LD_g[rs:(rs + sz - 1), rs:(rs + sz - 1)] <- ld_rho_within
    rs <- rs + sz
  }
  diag(LD_g) <- 1

  B1 <- matrix(0, nrow = p_gene, ncol = cg_blocks[1])
  B1[5, 1:8] <- signal;  B1[10, 1:8] <- signal; B1[15, 1:8] <- signal
  B1[20, 1:8] <- signal; B1[25, 1:8] <- signal; B1[30, 1:8] <- signal
  B1[35, 1:8] <- signal; B1[40, 1:8] <- signal; B1[45, 1:8] <- signal
  B1[50, 1:8] <- signal; B1[55, 1:8] <- signal; B1[60, 1:8] <- signal
  B1[65, 1:8] <- signal; B1[70, 1:8] <- signal; B1[75, 1:8] <- signal
  B1[80, 1:8] <- signal; B1[85, 1:8] <- signal; B1[90, 1:8] <- signal
  B1[95, 1:8] <- signal; B1[100, 1:8] <- signal
  rownames(B1) <- assign_column_names("Gene1_rs", p_gene)
  colnames(B1) <- assign_column_names("Gene1_cg", cg_blocks[1])

  B2 <- matrix(0, nrow = p_gene, ncol = cg_blocks[2])
  B2[1:5, 1] <- signal; B2[36:45, 3] <- signal; B2[96:100, 6] <- signal
  B2[c(1, 3, 5), 8] <- signal; B2[c(36, 38, 40, 42, 44), 11] <- signal
  B2[c(96, 98, 100), 1] <- signal; B2[c(37, 39, 41, 43, 45), 3] <- signal
  B2[c(97, 99), 7] <- signal; B2[c(2, 4), 9] <- signal
  B2[c(36, 40, 45), 6] <- signal
  rownames(B2) <- assign_column_names("Gene2_rs", p_gene)
  colnames(B2) <- assign_column_names("Gene2_cg", cg_blocks[2])

  B3 <- matrix(0, nrow = p_gene, ncol = sum(cg_blocks[3:4]))
  B3[1, 1:10] <- signal; B3[1:5, 3] <- signal; B3[c(2, 5), 4] <- signal
  B3[96:100, 3] <- signal; B3[c(97, 99), 5] <- signal; B3[96, 1:10] <- signal
  B3[36:45, 17] <- signal; B3[c(37, 39, 41, 43, 45), 19] <- signal
  B3[c(38, 41), 15] <- signal; B3[36, 11:20] <- signal; B3[44, 11:20] <- signal
  rownames(B3) <- assign_column_names("Gene3_rs", p_gene)
  colnames(B3) <- assign_column_names("Gene3_cg", sum(cg_blocks[3:4]))

  B <- as.matrix(Matrix::bdiag(B1, B2, B3))
  Z_names <- rownames(B) <- c(rownames(B1), rownames(B2), rownames(B3))
  X_names <- colnames(B) <- c(colnames(B1), colnames(B2), colnames(B3))

  B[35, 30] <- signal;  B[95, 25] <- signal;  B[75, 48] <- signal
  B[137, 5] <- signal;  B[198, 2] <- signal;  B[105, 42] <- signal
  B[237, 8] <- signal;  B[297, 4] <- signal;  B[204, 22] <- signal

  theta_true <- c(
    rep(0, 15),
    signal, 0, signal, 0, 0, signal, 0, signal, 0, 0, signal, 0, 0, 0, 0,
    signal, signal, signal, signal, signal, 0, 0, 0, 0, 0,
    0, 0, 0, signal, 0, signal, 0, signal, signal, signal
  )
  names(theta_true) <- X_names

  A <- simulate_confounder_effects_burgess_fixed(
    cg_blocks            = cg_blocks,
    gene_id              = methyl_gene_id,
    signal_strength      = conf_strength,
    sparsity_within_gene = 0.4,
    seed                 = 12345
  )
  colnames(A) <- X_names

  gamma <- rep(gamma_strength, 3)
  names(gamma) <- c("U1", "U2", "U3")

  if (is.null(seed)) set.seed(NULL) else set.seed(seed)

  SimSNP1 <- simulate_SNP_matrix(p_gene, n, LD_g, maf)
  SimSNP2 <- simulate_SNP_matrix(p_gene, n, LD_g, maf)
  SimSNP3 <- simulate_SNP_matrix(p_gene, n, LD_g, maf)
  Z <- cbind(SimSNP1$SNP, SimSNP2$SNP, SimSNP3$SNP)
  colnames(Z) <- Z_names

  U <- simulate_confounders_matrix(n, q = 3)

  SimMethyl <- simulate_methyl_sites_blocked_withU(
    cg_blocks, methyl_gene_id,
    SNP_matrix = Z, Coeff_matrix = B, sigma_methyl = noise,
    Conf_matrix = U, Conf_coeff_matrix = A,
    rho_within = rho_within, rho_decay = rho_decay,
    rho_block = rho_block, rho_gene = rho_gene,
    structure = cov_structure)
  X <- SimMethyl$X
  colnames(X) <- X_names

  epsilon_y <- rnorm(n, 0, noise)
  y <- as.numeric(X %*% theta_true + U %*% gamma) + epsilon_y

  site_block_id <- rep(seq_along(cg_blocks), times = cg_blocks)

  design <- list(
    ld_blocks      = ld_blocks,
    cg_blocks      = cg_blocks,
    methyl_gene_id = methyl_gene_id,
    site_block_id  = site_block_id,
    gene_structure = list(
      list(snp_range = 1:100,   methyl_range = 1:15,
           ld_blocks = ld_blocks_gene, cg_blocks = cg_blocks[1]),
      list(snp_range = 101:200, methyl_range = 16:30,
           ld_blocks = ld_blocks_gene, cg_blocks = cg_blocks[2]),
      list(snp_range = 201:300, methyl_range = 31:50,
           ld_blocks = ld_blocks_gene, cg_blocks = cg_blocks[3:4])
    ),
    block_names    = c("Gene1_NULL", "Gene2_SIG", "Gene3r1_SIG", "Gene3r2_SIG"),
    null_blocks    = c(TRUE, FALSE, FALSE, FALSE)
  )

  truth <- list(
    B               = B,
    A               = A,
    theta           = theta_true,
    gamma           = gamma,
    IVs_true        = apply(B, 1, function(r) l2norm(r) != 0),
    CpG_select_true = apply(B, 2, function(cl) l2norm(cl) != 0),
    CpG_effect_true = theta_true != 0
  )

  raw <- list(Z = Z, X = X, U = U, y = y)

  if (standardize) {
    std <- function(M) {
      M  <- as.matrix(M)
      mu <- colMeans(M)
      s  <- apply(M, 2, sd)
      s[s == 0] <- 1
      out <- sweep(sweep(M, 2, mu, "-"), 2, s, "/")
      attr(out, "center") <- mu
      attr(out, "scale")  <- s
      out
    }
    Z_s <- std(Z); X_s <- std(X); U_s <- std(U); y_s <- std(y)
    sZ <- attr(Z_s, "scale"); sX <- attr(X_s, "scale")
    sU <- attr(U_s, "scale"); sy <- attr(y_s, "scale")

    truth$theta_std <- theta_true * sX / sy
    truth$gamma_std <- gamma * sU / sy
    truth$B_std     <- sweep(sweep(B, 1, sZ, "*"), 2, sX, "/")
    truth$A_std     <- sweep(sweep(A, 1, sU, "*"), 2, sX, "/")

    y_out <- as.numeric(y_s)
    attr(y_out, "center") <- attr(y_s, "center")
    attr(y_out, "scale")  <- sy

    return(list(Z = Z_s, X = X_s, U = U_s, y = y_out,
                raw = raw, truth = truth, design = design,
                standardized = TRUE, seed = seed))
  }

  list(Z = Z, X = X, U = U, y = y,
       raw = raw, truth = truth, design = design,
       standardized = FALSE, seed = seed)
}