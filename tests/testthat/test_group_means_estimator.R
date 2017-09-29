
library(proDD)

context("Group Mean Estimator")

set.seed(1234)


test_that("estimate_group means works", {
  
  data <- generate_zero_inflated_data_with_effect(N_genes=100, N_rep=3, perc_changed = 0, mu0=8.2, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  group_locations <- estimate_group_means(X, design, 5, 0.4, 8, -0.3)
  expect_is(group_locations, "data.frame")
  expect_named(group_locations, c("mu", "mu_err", "beta_1", "beta_2", "sds", "sdd", "sdd_err", "max_ll"))
})

















