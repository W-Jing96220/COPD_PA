rm(list=ls())
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

# ------------------------------------- Read in full data --------------------------------------------------
# load full data
data <- readMM(paste0(work_dir, "GSE136831_RawCounts_Sparse.mtx")) # [1]  45947 312928
cell <- read.delim(paste0(work_dir,"GSE136831_AllCells.cellBarcodes.txt"), header = FALSE) # [1] 312928      1
gene <- read.delim(paste0(work_dir,"GSE136831_AllCells.GeneIDs.txt")) # [1] 45947     2
meta <- read.delim(paste0(work_dir,"GSE136831_AllCells.Samples.CellType.MetadataTable.txt")) # [1] 312928      9
# assign col- and row- name to data matrix
rownames(data) <- gene$HGNC_EnsemblAlt_GeneID
colnames(data) <- cell[,1]
cat("\n\ndata matrix head:\n")
head(data)
cat("\n\n")

# ------------------------------------- Data exploration/extraction --------------------------------------------------
diseaseType <- meta$Disease_Identity
unique(diseaseType) # "Control" "IPF" "COPD"
cat("\n\n")

# select meta based on Manuscript_Identity
cat("\nTake subset of meta, only get [", Subclass_CellType, "]:\n")
select_meta_by_Manuscript_Identity <- meta[meta$Manuscript_Identity == Subclass_CellType, ]
cat("\nSubset [", Subclass_CellType, "] shape:\n")
dim(select_meta_by_Manuscript_Identity)
cat("\nSubset [", Subclass_CellType, "] head:\n")
head(select_meta_by_Manuscript_Identity)
cat("\n\n")

# select meta based on group type
cat("Take Sub-subset of meta, only get [", Subclass_CellType, "]-[", groupType, "]:\n")
select_meta_by_Manuscript_Identity_then_by_groupType <- select_meta_by_Manuscript_Identity[select_meta_by_Manuscript_Identity$Disease_Identity %in% groupType, ]
cat("\nSub-subset [", Subclass_CellType, "]-[", groupType, "] shape:\n")
dim(select_meta_by_Manuscript_Identity_then_by_groupType)
cat("\nSub-subset [", Subclass_CellType, "]-[", groupType, "] head:\n")
head(select_meta_by_Manuscript_Identity_then_by_groupType)
cat("\n\n")

# display count of CellType_Category in the selected results
cat("Here are the CellType_Category column summary for [", Subclass_CellType, "]:\n")
table(select_meta_by_Manuscript_Identity$CellType_Category)
cat("\n\n")

# display count of Manuscript_Identity in the selected results
cat("Here are the Manuscript_Identity column summary for [", Subclass_CellType, "]:\n")
table(select_meta_by_Manuscript_Identity$Manuscript_Identity)
cat("\n\n")


# check if column Subclass_Cell_Identity is identical to column Manuscript_Identity
cat("Sub-subset column Subclass_Cell_Identity == Manuscript_Identity? :\n")
all(select_meta_by_Manuscript_Identity_then_by_groupType$Subclass_Cell_Identity == select_meta_by_Manuscript_Identity_then_by_groupType$Manuscript_Identity)

# display count of Manuscript_Identity in the selected results
cat("Here are the Manuscript_Identity column summary for subset [", Subclass_CellType, "]'s sub-subset (", groupType, "):\n")
table(select_meta_by_Manuscript_Identity_then_by_groupType$Manuscript_Identity)
cat("\n\n")

# display count of Subclass_Cell_Identity in the selected results
cat("Here are the Subclass_Cell_Identity column summary for subset [", Subclass_CellType, "]'s sub-subset (", groupType, "):\n")
table(select_meta_by_Manuscript_Identity_then_by_groupType$Subclass_Cell_Identity)
cat("\n\n")

# extract the CellBarcode_Identity from the sub-subset of meta
cat("Here is HEAD of the CellBarcode_Identity column of sub-subset of meta:\n")
select_cell <- select_meta_by_Manuscript_Identity_then_by_groupType$CellBarcode_Identity
head(select_cell)
cat("\n\n")

cat("Partition data matrix using CellBarcode_Identity column of sub-subset of meta:\n")
data_part_byGroup_tmp <- data[ ,select_cell]
cat("\nShape of current result to be saved into RDS:\n")
dim(data_part_byGroup_tmp)
cat("\nHead of current result to be saved into RDS:\n")
head(data_part_byGroup_tmp)
cat("\n\n")

saveRDS(data_part_byGroup_tmp, file = paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".rds"))
cat("RDS saved as: ", paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".rds\n"))

# transpose matrix, move gene to col and cell to row, required by Rmagic
cat("\nTranspose the data matrix:\n")
data_part_byGroup_tmp <- Matrix::t(data_part_byGroup_tmp)
cat("\nHead of transposed result to be saved into RDS:\n")
head(data_part_byGroup_tmp)
cat("\n\n")

saveRDS(data_part_byGroup_tmp, file = paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed.rds"))

cat("RDS saved as: ", paste0(Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed.rds\n"))