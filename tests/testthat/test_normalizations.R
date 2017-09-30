
library(proDD)

context("Normalization")


test_that("ma plot works", {
  data <- generate_zero_inflated_data_with_effect(N_genes=1000, N_rep=3, perc_changed = 0, mu0=8.2, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  data$Y <- data$Y * 1.2
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  
  ma_plot(X, design, reference_cond = 2)
  
})



test_that("median normalization works", {
  data <- generate_zero_inflated_data_with_effect(N_genes=1000, N_rep=3, perc_changed = 0, mu0=8.2, nu0=5, sigma0=0.4, location=8, scale=-0.3)
  data$Y <- data$Y * 1.2
  X <- cbind(data$X, data$Y)
  colnames(X) <- c(paste0("A_", 1:3), paste0("B_", 1:3))
  X <- X[rowSums(X) != 0, ]
  data_description <-  data.frame(Condition=as.factor(c(rep("A", 3), rep("B", 3))), 
                                  Replicate=c(1:3, 1:3))
  data_description$Sample <- paste0(data_description$Condition, data_description$Replicate)
  
  design <-  model.matrix(Sample ~ Condition - 1, data_description)
  
  X <- median_normalization(X, design, reference_cond = 2)
  ma_plot(X, design, reference_cond = 2)
  
})