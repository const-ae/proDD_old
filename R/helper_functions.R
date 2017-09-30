
#' @useDynLib proDD, .registration = TRUE 
NULL

#'
#'
#'Function that will be transformed with a function
#'  that is very close to log10 for values between 1-infinity, but has 
#'  the feature that pseudoLog10(0) = 0 and pseudoLog10(- x) = - pseudoLog10(x)
#' @param x a numeric vector 
#'
#' @export
pseudoLog10 <- function(x) { asinh(x/2)/log(10) }

p_sigmoid <- function(mu_i, location, scale) plogis((mu_i - location) / scale)

get_stan_model <- function(){

  stanmodels$estimate_sigmoid
  
}


scalar1 <- function(x) {x / sqrt(sum(x^2))}