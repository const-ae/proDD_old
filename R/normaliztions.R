

#' Median Normalization
#' 
#' @param X a numerical matrix that contains the log intensity values. Missing
#'   data has to be encoded as 0
#' @param design a design matrix that has as many columns as conditions
#'   and as many rows as samples (i.e. columns) in X. If a sample belongs
#'   to a condition it has value 1, if not the matrix contains a 0.
#' @param reference_cond the index of the design column with the condition
#'   that is used as the reference.
#'   
#' @export
median_normalization <- function(X, design,  reference_cond=1){
  ref_rowwise_mean_wo_zero <- apply(X[, which(design[ , reference_cond] != 0)], 1, function(row)mean(row[row != 0]))
  
  correction_fct <- sapply((1:ncol(design))[-reference_cond], function(idx){
    i <-  which(design[, idx] != 0)
    rowwise_mean_wo_zero <- apply(X[, i], 1, function(row)mean(row[row != 0]))
    median(rowwise_mean_wo_zero / ref_rowwise_mean_wo_zero, na.rm = TRUE)
  })

  for(idx in (1:ncol(design))[-reference_cond]){
    i <-  which(design[, idx] != 0)
    X[, i] <- X[, i] / correction_fct[idx]
  } 
  X
}

#' MA Plot
#' 
#' @param X a numerical matrix that contains the log intensity values. Missing
#'   data has to be encoded as 0
#' @param design a design matrix that has as many columns as conditions
#'   and as many rows as samples (i.e. columns) in X. If a sample belongs
#'   to a condition it has value 1, if not the matrix contains a 0.
#' @param reference_cond the index of the design column with the condition
#'   that is used as the reference.
#'   
#' @export
ma_plot <- function(X, design, reference_cond=1){
  ref_rowwise_mean_wo_zero <- apply(X[, which(design[ , reference_cond] != 0)], 1, function(row)mean(row[row != 0]))
  # Let's make new MA plots
  par(mfrow=c(ceiling(sqrt(ncol(design)-1)),ceiling((ncol(design)-1)/ceiling(sqrt(ncol(design)-1)))),
      oma = c(5,4,0,0) + 0.1,
      mar = c(0,0,1,1) + 0.1)
  for(idx in (1:ncol(design))[-reference_cond]){
    i <-  which(design[, idx] != 0)
    tmp_a <- apply(X[, i], 1, function(row)mean(row[row != 0]))
    plot(1/2 * (log2(tmp_a) + log2(ref_rowwise_mean_wo_zero)),log2(tmp_a) - log2(ref_rowwise_mean_wo_zero), xlab="A", ylab="M")
    abline(h=median(log2(tmp_a) - log2(ref_rowwise_mean_wo_zero), na.rm=TRUE), col="blue")
    abline(h=0, col="red")
  }
  title(xlab = "A",
        ylab = "M",
        outer = TRUE, line = 3)

}
