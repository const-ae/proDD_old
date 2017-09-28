


#' For each protein estimate latent intensity per condition
#'
#' MLE for the group mean of each condition in each protein, taking into account
#' the likelihood of observing a zero (quantified by the sigmoid).
#' The function can also deal with all zero groups.
#' @param X a numerical matrix that contains the log intensity values. Missing
#'   data has to be encoded as 0
#' @param design a design matrix that has as many columns as conditions
#'   and as many rows as samples (i.e. columns) in X. If a sample belongs
#'   to a condition it has value 1, if not the matrix contains a 0.
#' @param d0 estimated degrees of freedom of the variance moderation prior. Can be found
#'   using \code{estimate_variance_moderation}
#' @param s0 estimated variance of the variance moderation prior. Can be found
#'   using \code{estimate_variance_moderation}
#' @param location_est inflection point of the sigmoid. Can be calculated using
#'   \code{estimate_sigmoid}
#' @param scale_est steepness of the sigmoid. Usually negative and can be calculate
#'   using \code{estimate_sigmoid} 
#' 
#' @export
estimate_group_means <- function(X, design, d0, s0, location_est, scale_est) {
  
  
  N_genes <- nrow(X)
  N_rep <- ncol(X)
  N_cond <- ncol(design)
  
  
  
  cll <- function(par, x, design, location, scale){
    N_cond <- ncol(design)
    betas <- sapply(1:(N_cond-1), function(i) par[[paste0("beta_", i)]])
    betas <- c(betas, - sum(betas))
    mus <- par[["mu"]] + betas
    if(par[["sds"]] <= 0 | par[["sdd"]] <= 0){
      - Inf
    }else{
      sum(ifelse(x == 0,
                 log(p_sigmoid(par[["mu"]], location, scale)),
                 log(1-p_sigmoid(par[["mu"]], location, scale)) +
                   dnorm(x - par[["mu"]], sd=par[["sds"]], log = TRUE))) +
        sum(sapply(1:N_cond, function(i)
          sum(ifelse(x[which(design[, i] == 1)] == 0,
                     log(p_sigmoid(mus[i], location, scale)),
                     log(1-p_sigmoid(mus[i], location, scale)) +
                       dnorm(x[which(design[, i] == 1)] - mus[i], sd=par[["sdd"]], log = TRUE)))))
    }
  }
  
  
  
  
  result <- do.call(rbind, lapply(1:nrow(X), function(i){
    
    x <- X[i, ]
    
    if(sum(x) > 0){
      mu_start <- ifelse(all(x == 0), 0, mean(x[x != 0]))
      beta_start <- 0
      sds_start <- ifelse(sum(x != 0) <= 1, 1, sd(x[x!=0])+1)
      sdd_start <- sds_start
      max_ll <- maxLik::maxLik(cll, 
                               start=c("mu"=mu_start, "sds"=sds_start, "sdd"=sdd_start,
                                       sapply(1:(N_cond-1), function(j)structure(beta_start, names=paste0("beta_",j)))),
                               method="BFGS", x=x, design=design, location=location_est, scale=scale_est)
      if(max_ll$code != 0){
        warning(paste0("Parameter fit did not succeed (code=", max_ll$code, "). Message: ", max_ll$message, "\n"))
      }
      mu_err <- maxLik:::stdEr.maxLik(max_ll)[["mu"]]^2 
      sdd_err=maxLik:::stdEr.maxLik(max_ll)[["sdd"]]^2
      
      # In case there are multiple groups with all zero
      # there are multiple ways to distribute the beta
      # of which one is randomly chosen. To fix this I will 
      # here set the beta of all random groups to a common
      # value - there mean
      betas <- c(max_ll$estimate[grepl("beta_", names(max_ll$estimate))],
                 structure(- sum(max_ll$estimate[grepl("beta_", names(max_ll$estimate))]), names=paste0("beta_", N_cond)))
      all_zero_groups <- sapply(1:N_cond, function(j)sum(x[which(design[, j] == 1)]) == 0)
      betas[all_zero_groups] <- mean(betas[all_zero_groups])
      
      # return list with important variables
      c(mu=max_ll$estimate[["mu"]], mu_err=mu_err, betas,
        sds=max_ll$estimate[["sds"]], sdd=max_ll$estimate[["sdd"]], sdd_err=sdd_err, max_ll=max_ll$maximum)
    }else{
      # return same variables just that they are NA
      c(mu=NA, mu_err=NA, structure(rep(NA, N_cond), names=paste0("beta", 1:N_cond)),
        sds=NA, sdd=NA, sdd_err=NA, max_ll=NA)
    }
    
  }))
  
  result
  
}