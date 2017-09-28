

#' Estimate dependency between missing values and intensity
#' 
#' Function that estimates the increasing likelihood to observe a missing
#' value if the latent protein intensity is lower. The decrease is parameterized
#' as a sigmoid and its location and scale are estimated. The function constructs
#' regularized mean estimates for each protein in each condition. Those estimates
#' are used to find the best sigmoid function that explains the drop outs (i.e. zeros).
#' Internally the function uses MCMC (implemented with Stan), it can thus take
#' a considerable time to calculate the results (1000x6 matrix took ~ 2h).
#' 
#' @param X a numerical matrix that contains the log intensity values. Missing
#'   data has to be encoded as 0
#' @param data_description a dataframe with a `Condition` factor column which 
#'   assign each column of \code{X} to a condition
#' @param d0 estimated degrees of freedom of the variance moderation prior. Can be found
#'   using \code{estimate_variance_moderation}
#' @param s0 estimated variance of the variance moderation prior. Can be found
#'   using \code{estimate_variance_moderation}
#' @param n_iter the number of MCMC iterations
#' @param ... additional parameters that are passed on to \code{rstan::sampling}
estimate_sigmoid <- function(X, data_description, d0, s0, n_iter=2000, ...){
  if(d0 > 200)
    d0 <- 200
  
  design <- as.numeric(data_description$Condition)
  N_cond <- length(unique(design))
  sigmoid_model <- get_stan_model()
  fit <- rstan::sampling(sigmoid_model, 
                         data=list(Y=X,N_genes=nrow(X), N_samp=ncol(X), 
                                   N_cond=N_cond, design=design, f_nu=d0, f_sigma=s0), 
                         iter=n_iter, ...)
  
  list(location_est=mean(rstan::extract(fit, par='location')$location),
       scale_est=1/mean(rstan::extract(fit, par='scale_inv')$scale_inv),
       fit=fit)
}