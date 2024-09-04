rm(list=ls())
cat("\n**************************************************************************************************************************\n**************************************************************************************************************************\n\n")
######################



#######################


library(reticulate)
library(Matrix)
library(Rmagic)
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

# ------------------------------------- part 2 --------------------------------------------------
cat("\n============================================================================================================\n\n")
cat("RDS read in as: ", paste0(work_dir, Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed.rds\n"))
data_part_byGroup_tmp <- readRDS(paste0(work_dir, Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed.rds"))

# data preparation per Ref: http://htmlpreview.github.io/?https://github.com/KrishnaswamyLab/MAGIC/blob/master/Rmagic/inst/examples/bonemarrow_tutorial.html
# remove lowly expressed genes and cells with small library size
# keep genes expressed in at least 10 cells
keep_cols_gene <- colSums(data_part_byGroup_tmp > 0) > 10
cat("Number of genes that expressed in at least 10 cells:\n")
sum(keep_cols_gene)
cat("\n\n")

cat("Partition data based on number of genes that expressed in at least 10 cells:\n")
data_part_byGroup_tmp <- data_part_byGroup_tmp[, keep_cols_gene]
cat("\nShape of current partitioned data:\n")
dim(data_part_byGroup_tmp)
cat("\nHead of current partitioned data:\n")
head(data_part_byGroup_tmp)

# keep cells with at least 1000 UMIs
keep_rows_cell <- rowSums(data_part_byGroup_tmp) > 1000
sum(keep_rows_cell)

cat("\nFurther partition data based on number of cells that have at least 1000 UMIs:\n")
data_part_byGroup_tmp <- data_part_byGroup_tmp[keep_rows_cell, ]
cat("\nShape of current partitioned data:\n")
dim(data_part_byGroup_tmp)
cat("\nHead of current partitioned data:\n")
head(data_part_byGroup_tmp)

# Normalizing data
cat("\nNormalizing data:\n")
data_part_byGroup_tmp <- library.size.normalize(data_part_byGroup_tmp)
data_part_byGroup_tmp <- sqrt(data_part_byGroup_tmp)
cat("\nShape of current data:\n")
dim(data_part_byGroup_tmp)
cat("\nHead of current partitioned data:\n")
head(data_part_byGroup_tmp)

saveRDS(data_part_byGroup_tmp, file = paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_ready4MAGIC.rds"))
cat("Intermediate [", Subclass_CellType, "] transposed_ready4MAGIC rds saved as: ", paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_ready4MAGIC.rds"), "\n\n")

# ------------------------------------- part 3 --------------------------------------------------
# data_part_byGroup_tmp <- readRDS(paste0(Subclass_CellType, "_", Reduce(paste0, groupType), "_ready4MAGIC.rds"))
cat("\n============================================================================================================\n\n")
# get mem usage
cat("Current memory usage:\n")
gc(TRUE)
cat("\n\n")

# dropout correction - impute zero
# run MAGIC
data_part_byGroup_tmp_MAGIC <- magic(data_part_byGroup_tmp, genes='all_genes')
cat("Current memory usage:\n")
gc(TRUE)
cat("\n\n")

# Print a MAGIC object
Rmagic:::print.magic(data_part_byGroup_tmp_MAGIC)

cat("\nShape of current partitioned data:\n")
dim(data_part_byGroup_tmp_MAGIC$result)
cat("\nHead of current partitioned data:\n")
data_part_byGroup_tmp_MAGIC$result[1:6, 1:34]

saveRDS(data_part_byGroup_tmp_MAGIC, file = paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC.rds"))
cat("The [", Subclass_CellType, "] transposed_MAGIC rds is saved as: ", paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC.rds"), "\n\n")

cat("Current memory usage:\n")
gc(TRUE)
cat("\n\n")

# ------------------------------------- part 4 --------------------------------------------------
cat("\n============================================================================================================\n\n")

# Convert a MAGIC object to a matrix
cat("Start to convert a MAGIC object to matrix...\n")
data_part_byGroup_tmp_MAGIC_res <- Rmagic:::as.matrix.magic(data_part_byGroup_tmp_MAGIC)
cat("Current memory usage:\n")
gc(TRUE)
cat("\n\n")

cat("\nShape of current data:\n")
dim(data_part_byGroup_tmp_MAGIC_res)
cat("\nHead of current data:\n")
data_part_byGroup_tmp_MAGIC_res[1:6, 1:34]

saveRDS(data_part_byGroup_tmp_MAGIC_res, file = paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_ori.rds"))
cat("The [", Subclass_CellType, "] transposed_MAGIC_matrix_ori rds (matrix with negative values) is saved as: ", paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_ori.rds"), "\n\n")

# Replace any resulting negative expression with 0, Ref: doi: 10.1093/nargab/lqaa002
cat("Start to replace any resulting negative expression from MAGIC process with 0...\n")
replaceNegative <- function(x) (abs(x)+x)/2
data_part_byGroup_tmp_MAGIC_processed <- replaceNegative(data_part_byGroup_tmp_MAGIC_res)

cat("\nShape of current data:\n")
dim(data_part_byGroup_tmp_MAGIC_processed)
cat("\nHead of current data:\n")
data_part_byGroup_tmp_MAGIC_processed[1:6, 1:34]

saveRDS(data_part_byGroup_tmp_MAGIC_processed, file = paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_processed.rds"))
cat("Processed [", Subclass_CellType, "] transposed_MAGIC_matrix_processed rds savedas: ", paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_processed.rds"), "\n\n")

