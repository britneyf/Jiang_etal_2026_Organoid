# Load libraries
library("anndata")
library("DESeq2")
library("reticulate")

# Use the path to your Python executable in the virtual environment
use_python("/lila/home/forsythb/.virtualenvs/r-reticulate/bin/")

# Look at where Python is located
py_config()

# Import scanpy
sc <- import("scanpy")

# Read h5ad file
adata <- sc$read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/postprocess_adata/postprocess_adata.020524/adata.combined.postprocess.h5ad')

# Extract count matrix from the raw layer
count_matrix <- adata$raw$X
count_matrix <- t(count_matrix)
count_matrix <- round(count_matrix)

# Define the column data
coldata <- adata$obs[,c("Tumor_Site","Culture_Media", "ZFP_Expression")]

# Convert column data to factors
coldata$Tumor_Site <- factor(coldata$Tumor_Site)
coldata$Culture_Media <- factor(coldata$Culture_Media)
coldata$ZFP_Expression <- factor(coldata$ZFP_Expression)

# Design a DESeqDataSet
dds_1 <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ Tumor_Site + Culture_Media + ZFP_Expression
)

# Add a pseudo-count of 1 to all counts
count_matrix <- counts(dds_1) + 1

# Specify the row names as the gene names
rownames(count_matrix) <- adata$raw$var_names

# Set the reference levels
coldata$Tumor_Site<-relevel(coldata$Tumor_Site,ref="Primary")
coldata$Culture_Media<-relevel(coldata$Culture_Media,ref="BASE")
coldata$ZFP_Expression<-relevel(coldata$ZFP_Expression,ref="CTRL")

# Design a DESeqDataSet with modified count matrix
dds_1 <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ Tumor_Site + Culture_Media + ZFP_Expression 
)

# Perform DESeq analysis
dds_1 <- DESeq(dds_1)

# Extract results
result_names_1 <- resultsNames(dds_1)
result_names_1

# Specify the contrast and extract results
res_1 <- results(dds_1)
res_1

# Specify the date
date <- '020524'

# Specify the directory path
directory_path <- paste("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/deseq/deseq.", date, sep="")

# Check if the directory exists, and create it if not
if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}

# Specify the file path where you want to save the results
result_file <- paste(directory_path, "/sc_zfpexp_results.csv", sep="")

# Write the results table to a CSV file
write.csv(as.data.frame(res_1), file = result_file)

# Save the DESeqDataSet object
saveRDS(dds_1, file = paste(directory_path, "/sc_zfpexp_dds.rds", sep=""))

# Save the DESeq results
saveRDS(results(dds_1), file = paste(directory_path, "/sc_zfpexp_results.rds", sep=""))
