
library(proDD)

context("Sigmoid Estimator")

set.seed(1234)

test_that("trivial test succeed", {
  expect_equal(7, 56/8)
  expect_equal(4, 2 * 2)
})

test_that("retrieving stan model works", {
  model <- get_stan_model()
  # Should fail if not a real model
  model@mk_cppmodule(model)
})




test_that("retrieving stan model works", {
  model <- get_stan_model()
  # Should fail if not a real model
  model@mk_cppmodule(model)
})


test_that("sigmoid_estimation works", {
  sigmoid_est_model <- get_stan_model()
  
  data <- generate_zero_inflated_data_with_effect(N_genes=10, N_rep=3, perc_changed = 0, mu0=8.2, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  

  design <- as.numeric(data_description$Condition)
  N_cond <- length(unique(design))
  fit <- suppressWarnings(rstan::sampling(sigmoid_est_model, 
                         data=list(Y=X, N_genes=nrow(X), N_samp=ncol(X), N_cond=N_cond, design=design, f_nu=5, f_sigma=0.4), 
                         iter=2000, control=list(max_treedepth=10), chains=1, seed=1234,
                         show_messages=FALSE))
  location_est=mean(rstan::extract(fit, par='location')$location)
  scale_est=1/mean(rstan::extract(fit, par='scale_inv')$scale_inv)
  expect_equal(location_est, 6.899262763)
  expect_equal(scale_est, -2.058711656)
})



test_that("estimate_sigmoid works", {

  data <- generate_zero_inflated_data_with_effect(N_genes=10, N_rep=3, perc_changed = 0, mu0=8.2, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)

  sig_est <- suppressWarnings(estimate_sigmoid(X, data_description, d0 = 5, s0=0.4, chains=1, seed=1234))
  expect_equal(sig_est$location_est, 6.250291, tolerance=1e-6)
  expect_equal(sig_est$scale_est, -1.950324, tolerance=1e-6)
})


test_that("sigmoid_estimation works with bad variance hyperparameters", {
  sigmoid_est_model <- get_stan_model()
  
  data <- generate_zero_inflated_data_with_effect(N_genes=50, N_rep=3, perc_changed = 0, mu0=8.2,
                                                  nu0=5, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  
  sig_est <- suppressWarnings(estimate_sigmoid(X, data_description, d0 = Inf, s0=0.4, chains=1, seed=1234))
  expect_equal(sig_est$location_est, 7.010875, tolerance=1e-6)
  expect_equal(sig_est$scale_est, -0.9440992, tolerance=1e-6)
})






