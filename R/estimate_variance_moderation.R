
#' Estimate variance moderation prior
#'
#' Function to estimate the parameters of an inverse chi-squared distribution
#' that is assummed as the prior of the variance for each protein.
#' The method was initially developed in Smyth, 2004 in the context of microarrays
#' to increase the power to detected differentially expressed genes. The
#' estimated inverse chi-squared distribution can be used to moderate the
#' variance estimate for each individual protein.
#' @param X a numerical matrix that contains the log intensity values. Missing
#'   data can either be encoded as 0 or NA.
#' @param design a design matrix that has as many columns as conditions
#'   and as many rows as samples (i.e. columns) in X. If a sample belongs
#'   to a condition it has value 1, if not the matrix contains a 0.
#' 
#' @export
estimate_variance_moderation <- function(X, design){
  '
  Identical to the result that limma produces
  '
  
  sum_across_cond <-  X %*% design
  mask <- X != 0 | is.na(X)
  nz <- mask %*% design
  mean_cond <- sum_across_cond / nz
  
  tmp <- sapply(1:nrow(X), function(i){
    x <- X[i, ]
    ss <- sum(sapply(1:ncol(design), function(j)sum(((x - mean_cond[i, j]) * design[, j] * mask[i, ])^2)), na.rm=TRUE)
    sum_nz <- sum(mask[i, ])
    n_cond <- ncol(design)
    if(sum_nz - n_cond > 0){
      c(ss/(sum_nz-n_cond), sum_nz-n_cond)
    }else{
      c(NA, NA)
    }
  })
  
  
  sg <- tmp[1, ! is.na(tmp[1, ])]
  dg <- tmp[2, ! is.na(tmp[2, ])]
  
  n <- length(dg)
  
  zg <- log(sg)
  eg <- zg - digamma(dg/2) + log(dg/2)
  
  emean <- mean(eg, na.rm=TRUE)
  evar <- sum((eg - emean)^2, na.rm=TRUE)/(n - 1)
  evar <- evar - mean(trigamma(dg/2), na.rm=TRUE)
  
  if(evar > 0){
    d0 <- 2 * limma:::trigammaInverse(evar)
    s0 <- exp(emean + digamma(d0/2) - log(d0/2))
  }else{
    d0 <- Inf
    s0 <- exp(emean)
  }
  
  return(list(nu_est=d0, sigma2_est=s0))
}
