rm(list=ls())
cat("\n**************************************************************************************************************************\n**************************************************************************************************************************\n\n")
library(Matrix)
# library(Rmagic)
# ------------------------------------- Read in / Set parameters --------------------------------------------------
# for GSE136831 data - 18 fabric，18 COPD patients, 28 health control, total 312k cells
args = commandArgs(trailingOnly=TRUE)
Subclass_CellType <- as.character(args[1]) # the cell sub type
# Subclass_CellType = c("Macrophage", "Macrophage_Alveolar", "cMonocyte", "ncMonocyte")
Subclass_CellType
work_dir <- as.character(args[2])
work_dir <- paste0(work_dir, "/")
work_dir
groupType <- c("COPD", "Control")
# groupType
# ------------------------------------- part 6 --------------------------------------------------
cat("\n============================================================================================================\n\n")

# calculate bayesian correlation for all gene combination
BaCo <- function(X){
  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE) 
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}

for (curr_groupType in groupType) {
    cat("RDS read in as: ", paste0(work_dir, Subclass_CellType, "_", curr_groupType, ".transposed_MAGIC_matrix_processed.rds"))
    data_MAGIC_processed_part_byGroup <- readRDS(paste0(work_dir, Subclass_CellType, "_", curr_groupType, ".transposed_MAGIC_matrix_processed.rds"))
    cat("\n\nDim info of [", Subclass_CellType, "]-[", curr_groupType, "] transposed_MAGIC_matrix_processed rds:\n")
    print(dim(data_MAGIC_processed_part_byGroup))
    cat("\n\n")

    # BaCo works with genes as row, so have to transpose.
    data_MAGIC_processed_correlationMTX <- BaCo(t(data_MAGIC_processed_part_byGroup))
    cat("Dim info of transposed _MAGIC_processed rds bayesian corr matrix:\n")
    print(dim(data_MAGIC_processed_correlationMTX))
    cat("\n\n")
    cat("\nHead of current result to be saved into RDS:\n")
    print(head(data_MAGIC_processed_correlationMTX, c(6L, 6L)))
    cat("\n\n")

    cat("Current memory usage:\n")
    system('free -h')
    gc(TRUE)
    system('free -h')
    saveRDS(data_MAGIC_processed_correlationMTX, file = paste0(Subclass_CellType, "_", curr_groupType, "_MAGIC_processed_correlationMTX.rds"))
    cat("RDS saved as", paste0(Subclass_CellType, "_", curr_groupType, "_MAGIC_processed_correlationMTX.rds\n"))
}