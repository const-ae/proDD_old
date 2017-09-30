

detect_differences <- function(X, design, data_description, d0, s0, group_locations, comparison){

  N_genes <- nrow(X)
  N_rep <- ncol(X)
  N_cond <- ncol(design)
  
  # Calculate df for each gene
  tmp <- sapply(1:N_genes, function(i){
    xz <- sum(X[i, ] != 0)
    complete_group_zero <- sum(X[i, ] %*% design == 0)
    c(max(0, xz - (N_cond - complete_group_zero)), max(0, xz))
  })
  df <- tmp[1, ]
  
  # Convert MLE to unbiased n/(n-2)
  tmp_p_vars <- group_locations[, "sdd"]^2 * tmp[2, ]/max(1,df)
  
  # Correct var estimate with global variance information
  var_j_est <- if(is.infinite(d0)){
    rep(s0, N_genes)
  }else{
    (tmp_p_vars  * df + s0 * d0) / (df + d0)
  } 
  
  
  betas <- as.matrix(group_locations[, grepl("beta_\\d+", colnames(group_locations))])
  
  contr_col <- limma::makeContrasts(paste0(comparison[1], "-", comparison[2]), levels=levels(data_description$Condition))[, 1]
  n_contr <- length(which(contr_col != 0))
  col_oi_sel <- rowSums(design[ ,which(contr_col != 0)]) != 0
  
  nnz <- sapply(1:N_genes, function(i)sum(X[i, col_oi_sel] != 0))
  
  mod_p <-  sapply(1:N_genes, function(i){
    nzeros_per_group <- sapply(which(contr_col != 0), function(j)sum(X[i, which(design[, j] != 0)] != 0))
    
    between_var <- ((betas[i, which(contr_col != 0)] ) %*%
                      scalar1(contr_col[which(contr_col != 0)]))^2 * mean(nzeros_per_group)
    within_var <- var_j_est[i] 
    r <- between_var / within_var 
    names(r) <- NULL
    (1 - pf(r, df1=2 - 1, df2=df[i] + d0))
  })
  
  data.frame(df, var_post=var_j_est, var_prior=tmp_p_vars, p_value=mod_p)
}