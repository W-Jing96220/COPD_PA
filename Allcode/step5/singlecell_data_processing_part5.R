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

# ------------------------------------- part 7 --------------------------------------------------
cat("\n============================================================================================================\n\n")

curr_groupType = groupType[1] # COPD
cat("RDS read in as: ", paste0(work_dir, Subclass_CellType, "_", curr_groupType, "_MAGIC_processed_correlationMTX.rds"))
data_MAGIC_processed_correlationMTX_COPD <- readRDS(paste0(work_dir, Subclass_CellType, "_", curr_groupType, "_MAGIC_processed_correlationMTX.rds"))
cat("\n\nDim info of [", Subclass_CellType, "]-[", curr_groupType, "] correlationMTX rds:\n")
print(dim(data_MAGIC_processed_correlationMTX_COPD))
cat("\n\n")

curr_groupType = groupType[2] # Control
cat("RDS read in as: ", paste0(work_dir, Subclass_CellType, "_", curr_groupType, "_MAGIC_processed_correlationMTX.rds"))
data_MAGIC_processed_correlationMTX_Control <- readRDS(paste0(work_dir, Subclass_CellType, "_", curr_groupType, "_MAGIC_processed_correlationMTX.rds"))
cat("\n\nDim info of [", Subclass_CellType, "]-[", curr_groupType, "] correlationMTX rds:\n")
print(dim(data_MAGIC_processed_correlationMTX_Control))
cat("\n\n")

# sanity check, to see if COPD and Control have same gene set (which it should)
Control_COPD_common_geneSet = intersect(colnames(data_MAGIC_processed_correlationMTX_COPD), colnames(data_MAGIC_processed_correlationMTX_Control))
cat("Control_COPD_common_geneSet length (should equal col num from dim info): ", length(Control_COPD_common_geneSet), "\n\n")

# Transform corr mtx to pairwise list
cat("Transform COPD corr matrix to pairwise list...\n")
data_MAGIC_processed_correlationMTX_COPD[lower.tri(data_MAGIC_processed_correlationMTX_COPD, diag=TRUE)] = NA # put NA to lower triangle including diag
cat("Dim info of COPD bayesian corr matrix (make lower triangle including diag all NA):\n")
print(dim(data_MAGIC_processed_correlationMTX_COPD))
data_MAGIC_processed_correlationMTX_COPD_toList <- as.data.frame(as.table(data_MAGIC_processed_correlationMTX_COPD)) # as a dataframe
cat("Dim info of COPD pairwise list (direct transform from matrix, so still with NA):\n")
print(dim(data_MAGIC_processed_correlationMTX_COPD_toList))
data_MAGIC_processed_correlationMTX_COPD_toList <- na.omit(data_MAGIC_processed_correlationMTX_COPD_toList) # remove NA
cat("Dim info of COPD pairwise list (NA removed, length should be n(n-1)/2):\n")
print(dim(data_MAGIC_processed_correlationMTX_COPD_toList)) # length should be n(n-1)/2
cat("Take absolute value for all the correlation...\n")
data_MAGIC_processed_correlationMTX_COPD_toList$Freq = abs(data_MAGIC_processed_correlationMTX_COPD_toList$Freq)
cat("\n\n")

cat("Transform Control corr matrix to pairwise list...\n")
data_MAGIC_processed_correlationMTX_Control[lower.tri(data_MAGIC_processed_correlationMTX_Control, diag=TRUE)] = NA # put NA to lower triangle including diag
cat("Dim info of Control bayesian corr matrix (make lower triangle including diag all NA):\n")
print(dim(data_MAGIC_processed_correlationMTX_Control))
data_MAGIC_processed_correlationMTX_Control_toList <- as.data.frame(as.table(data_MAGIC_processed_correlationMTX_Control)) # as a dataframe
cat("Dim info of Control pairwise list (direct transform from matrix, so still with NA):\n")
print(dim(data_MAGIC_processed_correlationMTX_Control_toList))
data_MAGIC_processed_correlationMTX_Control_toList <- na.omit(data_MAGIC_processed_correlationMTX_Control_toList) # remove NA
cat("Dim info of Control pairwise list (NA removed, length should be n(n-1)/2):\n")
print(dim(data_MAGIC_processed_correlationMTX_Control_toList)) # length should be n(n-1)/2
cat("Take absolute value for all the correlation...\n")
data_MAGIC_processed_correlationMTX_Control_toList$Freq = abs(data_MAGIC_processed_correlationMTX_Control_toList$Freq)
cat("\n\n")

cat("Current memory usage:\n")
gc(TRUE)
system('free -h')
cat("\n\n")
   
saveRDS(data_MAGIC_processed_correlationMTX_COPD_toList, file = paste0(Subclass_CellType, "_", groupType[1], "_MAGIC_processed_correlationMTX_toList.rds"))
cat(groupType[1], " RDS saved as", paste0(Subclass_CellType, "_", groupType[1], "_MAGIC_processed_correlationMTX_toList.rds\n"))
saveRDS(data_MAGIC_processed_correlationMTX_Control_toList, file = paste0(Subclass_CellType, "_", groupType[2], "_MAGIC_processed_correlationMTX_toList.rds"))
cat(groupType[2], " RDS saved as", paste0(Subclass_CellType, "_", groupType[2], "_MAGIC_processed_correlationMTX_toList.rds\n"))

cat("COPD and Control two _MAGIC_processed_correlationMTX_toList rds saved!\n\n")

# ------------------------------------- part 8 --------------------------------------------------
cat("Start to save COPD and Control two _MAGIC_processed_correlationMTX_toList rds into txt...\n\n")

# data_MAGIC_processed_correlationMTX_COPD_toList <- readRDS(paste0(Subclass_CellType, "_", groupType[1], "_MAGIC_processed_correlationMTX_toList.rds"))
write.table(data_MAGIC_processed_correlationMTX_COPD_toList, paste0(Subclass_CellType, "_", groupType[1], "_MAGIC_processed_correlationMTX_toList.txt"), append = FALSE, sep = "\t", quote = FALSE, dec = ".", row.names = FALSE, col.names = TRUE)
cat(groupType[1], " txt saved as", paste0(Subclass_CellType, "_", groupType[1], "_MAGIC_processed_correlationMTX_toList.txt\n"))

# data_MAGIC_processed_correlationMTX_Control_toList <- readRDS(paste0(Subclass_CellType, "_", groupType[2], "_MAGIC_processed_correlationMTX_toList.rds"))
write.table(data_MAGIC_processed_correlationMTX_Control_toList, paste0(Subclass_CellType, "_", groupType[2], "_MAGIC_processed_correlationMTX_toList.txt"), append = FALSE, sep = "\t", quote = FALSE, dec = ".", row.names = FALSE, col.names = TRUE)
cat(groupType[2], " txt saved as", paste0(Subclass_CellType, "_", groupType[2], "_MAGIC_processed_correlationMTX_toList.txt\n"))

cat("COPD and Control two _MAGIC_processed_correlationMTX_toList txt saved!\n\n")