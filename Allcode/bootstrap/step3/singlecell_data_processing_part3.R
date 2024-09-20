rm(list=ls())
cat("\n**************************************************************************************************************************\n**************************************************************************************************************************\n\n")
library(Matrix)
# library(Rmagic)
# ------------------------------------- Read in / Set parameters --------------------------------------------------
# for GSE136831 data - 18 fabric, 18 COPD patients, 28 health control, total 312k cells
args = commandArgs(trailingOnly=TRUE)
Subclass_CellType <- as.character(args[1]) # the cell sub type
# Subclass_CellType = c("Macrophage", "Macrophage_Alveolar", "cMonocyte", "ncMonocyte")
Subclass_CellType
work_dir <- as.character(args[2]) 
work_dir <- paste0(work_dir, "/")
work_dir
groupType <- c("COPD", "Control")
# groupType

# ------------------------------------- Define the splitting function --------------------------------------------------
split_data_with_identified_item <- function(data, identified_index, split_ratio = 0.8, n_splits = 10) {
  # List to store each split's train data
  split_list <- list()
  
  # Number of rows in the dataset
  total_rows <- nrow(data)
  
  # Calculate split size (subtract 1 to account for the identified item)
  split_size <- round(split_ratio * total_rows) - 1
  
  # Ensure split_size is non-negative
  split_size <- max(split_size, 0)

  # Loop to generate multiple splits
  for (i in 1:n_splits) {
    set.seed(i)  # Set a different seed for each split
    
    # Ensure the identified item is included
    remaining_indices <- setdiff(1:total_rows, identified_index)  # Exclude the identified item
    
    # Sample indices while ensuring we do not exceed the available indices
    if (split_size > length(remaining_indices)) {
      stop("Split size exceeds the number of available indices.")
    }
    
    sampled_indices <- sample(remaining_indices, split_size, replace = FALSE)
    
    # Combine the identified item with the sampled data
    final_indices <- c(identified_index, sampled_indices)
    
    # Create the split
    train_data <- data[final_indices, ]
    
    # Save the current split to the list
    split_list[[i]] <- train_data
  }
  
  # Return the list of splits
  return(split_list)
}

# ------------------------------------- part 5 --------------------------------------------------
cat("\n============================================================================================================\n\n")
cat("RDS read in as: ", paste0(work_dir, Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_processed.rds"))
data_part_byGroup_tmp_MAGIC_processed <- readRDS(paste0(work_dir, Subclass_CellType, "_", Reduce(paste0, groupType), ".transposed_MAGIC_matrix_processed.rds"))
cat("Dim info of [", Subclass_CellType, "] transposed_MAGIC_matrix_processed rds:\n")
dim(data_part_byGroup_tmp_MAGIC_processed)
cat("\n\n")

cat("Read in Meta...\n")
meta <- read.delim("/lustre06/project/6040537/jingw/GSE136831_AllCells.Samples.CellType.MetadataTable.txt")
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

    # Split the data into 10 subsamples while keeping the first item
    identified_index <- 1  # Change this if needed to keep a different item
    splits <- split_data_with_identified_item(data_MAGIC_processed_part_byGroup, identified_index, split_ratio = 0.8, n_splits = 10)

    # Save each split as an RDS file
    for (i in 1:length(splits)) {
        saveRDS(splits[[i]], file = paste0(Subclass_CellType, "_", curr_groupType, "_split_", i, ".rds"))
        cat("RDS saved as: ", paste0(Subclass_CellType, "_", curr_groupType, "_split_", i, ".rds"), "\n")
    }
}
