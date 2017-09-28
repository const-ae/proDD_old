library(proDD)

context("Variance Moderation Estimator")

set.seed(1234)


test_that("estimate_variance_moderation works", {

  data <- generate_zero_inflated_data_with_effect(N_genes=1000, N_rep=3, perc_changed = 0, mu0=8.2, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  vm_est <- estimate_variance_moderation(X, design)
  expect_equal(vm_est$nu_est, 6.936226, tolerance=1e-6)
  expect_equal(vm_est$sigma2_est, 0.4204029, tolerance=1e-6)
})



test_that("estimate_variance_moderation can return infinity", {
  
  data <- generate_zero_inflated_data_with_effect(N_genes=100, N_rep=3, perc_changed = 0, mu0=8.2, nu0=50, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  vm_est <- estimate_variance_moderation(X, design)
  expect_equal(vm_est$nu_est, Inf, tolerance=1e-6)
  expect_equal(vm_est$sigma2_est, 0.574127, tolerance=1e-6)
})


test_that("estimate_variance_moderation is accurate if large change", {
  
  data <- generate_zero_inflated_data_with_effect(N_genes=1000, N_rep=6, perc_changed = 0.8, lfc_sd=20, mu0=8.2, nu0=5, sigma0=0.4)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:6), paste0("B_", 1:6))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 6), rep("B", 6))), 
                                  Replicate=c(1:6, 1:6))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  vm_est <- estimate_variance_moderation(X, design)
  expect_equal(vm_est$nu_est,  5.709778, tolerance=1e-6)
  expect_equal(vm_est$sigma2_est, 0.5365091, tolerance=1e-6)
})
