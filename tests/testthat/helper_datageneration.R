


generate_zero_inflated_data_with_effect <- function(N_genes=1000, N_rep=20, lfc_sd=1, perc_changed=0.1, 
                                                    mu0=9, mu0_sd=0.4, location=8, scale=-0.3, nu0=40, sigma0=0.6){
  real_mus <- rnorm(n=N_genes, mu0, mu0_sd)
  changed_label <- runif(N_genes) < perc_changed
  log_fold_change <- ifelse(changed_label, rnorm(N_genes, mean=0, sd=lfc_sd),  rnorm(N_genes, mean=0, sd=0))
  real_mus1 <- real_mus + log_fold_change/2
  real_mus2 <- real_mus - log_fold_change/2
  var_j <- geoR::rinvchisq(n=N_genes, nu0, sigma0)
  
  X <- t(sapply(1:N_genes, function(i)ifelse(runif(N_rep) < p_sigmoid(real_mus1[i], location, scale), 0, rnorm(n=N_rep, real_mus1[i], sd=sqrt(var_j[i])))))
  Y <- t(sapply(1:N_genes, function(i)ifelse(runif(N_rep) < p_sigmoid(real_mus2[i], location, scale), 0, rnorm(n=N_rep, real_mus2[i], sd=sqrt(var_j[i])))))
  list(X=X, Y=Y, changed_label=changed_label)
}
