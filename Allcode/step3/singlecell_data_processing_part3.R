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
# ------------------------------------- part 5 --------------------------------------------------
cat("\n============================================================================================================\n\n")
cat("RDS read in as: ", paste0(work_dir, Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_processed.rds"))
data_part_byGroup_tmp_MAGIC_processed <- readRDS(paste0(work_dir, Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_processed.rds"))
cat("Dim info of [", Subclass_CellType, "] transposed_MAGIC_matrix_processed rds:\n")
dim(data_part_byGroup_tmp_MAGIC_processed)
cat("\n\n")

cat("Read in Meta...\n")
######
meta <- read.delim
###### need change here
diseaseType <- meta$Disease_Identity
unique(diseaseType) # "Control" "IPF" "COPD"
cat("\n\n")

for (curr_groupType in groupType) {
    cat("\n-------------------------------------------------------------------------------------------------------------------\n")
    # select meta based on group type
    cat("\nTake subset of meta, only get [", curr_groupType, "]:\n")
    select_meta_by_groupType <- meta[meta$Disease_Identity == curr_groupType, ]
    cat("\nSubset meta's [", curr_groupType, "] shape:\n")
    print(dim(select_meta_by_groupType))
    cat("\nSubset meta's [", curr_groupType, "] head:\n")
    print(head(select_meta_by_groupType))
    cat("\n\n")

    # further select meta based on Manuscript_Identity
    cat("\nTake subset of select_meta_by_groupType, only get [", Subclass_CellType, "]:\n")
    select_meta_by_groupType_then_by_Manuscript_Identity <- select_meta_by_groupType[select_meta_by_groupType$Manuscript_Identity == Subclass_CellType, ]
    cat("\nSub-subset meta's [", curr_groupType, "]-[", Subclass_CellType, "] shape:\n")
    print(dim(select_meta_by_groupType_then_by_Manuscript_Identity))
    cat("\nSub-subset meta's [", curr_groupType, "]-[", Subclass_CellType, "] head:\n")
    print(head(select_meta_by_groupType_then_by_Manuscript_Identity))
    cat("\n\n")

    # extract the CellBarcode_Identity from the sub-subset of meta
    cat("Here is HEAD of the CellBarcode_Identity column of sub-subset of meta:\n")
    select_cell <- select_meta_by_groupType_then_by_Manuscript_Identity$CellBarcode_Identity
    print(head(select_cell))
    cat("\n\n")

    cat("Partition data matrix using CellBarcode_Identity column of sub-subset of meta:\n")
    data_MAGIC_processed_part_byGroup <- data_part_byGroup_tmp_MAGIC_processed[rownames(data_part_byGroup_tmp_MAGIC_processed) %in% select_cell, ]
    cat("\nShape of current result to be saved into RDS:\n")
    print(dim(data_MAGIC_processed_part_byGroup))
    cat("\nHead of current result to be saved into RDS:\n")
    print(head(data_MAGIC_processed_part_byGroup, c(6L, 34L)))
    cat("\n\n")

    saveRDS(data_MAGIC_processed_part_byGroup, file = paste0(Subclass_CellType, "_", curr_groupType, ".transposed_MAGIC_matrix_processed.rds"))
    cat("RDS saved as: ", paste0(Subclass_CellType, "_", curr_groupType, ".transposed_MAGIC_matrix_processed.rds"))
}




