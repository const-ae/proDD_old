
library(proDD)

context("Difference Detection")


test_that("estimate_group means works", {
  
  data <- generate_zero_inflated_data_with_effect(N_genes=100, N_rep=3, perc_changed = 0, mu0=8.7, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  group_locations <- estimate_group_means(X, design, 5, 0.4, 8, -0.3)
  result <- detect_differences(X, design, data_description, d0=5, s0=0.4, group_locations=group_locations, comparison=c("A", "B"))
  
  plot(ppoints(length(result$p_value)), sort(result$p_value), log="xy"); abline(0,1)
})


test_that("complete workflow works", {
  
  data <- generate_zero_inflated_data_with_effect(N_genes=100, N_rep=3, perc_changed = 0, mu0=9, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  vm_est <- estimate_variance_moderation(X, design)
  sig_est <- estimate_sigmoid(X, data_description, vm_est$nu_est, vm_est$sigma2_est, chains=1)
  group_locations <- estimate_group_means(X, design, vm_est$nu_est, vm_est$sigma2_est, sig_est$location_est, sig_est$scale_est)
  result <- detect_differences(X, design, data_description, d0=vm_est$nu_est, s0=vm_est$sigma2_est,
                               group_locations=group_locations, comparison=c("A", "B"))
  print(list( vm_est$nu_est, vm_est$sigma2_est, sig_est$location_est, sig_est$scale_est))
  plot(ppoints(length(result$p_value)), sort(result$p_value), log="xy"); abline(0,1)
})
