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

# Read in the adata for primary tumor site
adata_primary = sc$read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/Organoid/input/degs/primary_tumor.h5ad')

# Get the counts matrix
counts = as.matrix(adata_primary$X)

# Normalize
ms <- rowSums(counts)
norm_df <- counts/ms * median(ms)

# Log transform
counts <- log(norm_df + 1)

# Add pseudocount
counts <- counts + 1

# Round since error that not all values are integers
counts <- round(counts)

# Transpose
counts <- t(counts)

# Assign row and column names to the counts matrix
rownames(counts) <- adata_primary$var_names
colnames(counts) <- adata_primary$obs_names

# Define the column data
coldata <- adata_primary$obs[,c("Culture_Media", "ZFP_Expression")]

# Convert column data to factors
coldata$Culture_Media <- factor(coldata$Culture_Media)
coldata$ZFP_Expression <- factor(coldata$ZFP_Expression)

# Set the reference levels
coldata$Culture_Media<-relevel(coldata$Culture_Media,ref="BASE")
coldata$ZFP_Expression<-relevel(coldata$ZFP_Expression,ref="CTRL")

# Design a DESeqDataSet with modified count matrix
dds_1 <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ Culture_Media + ZFP_Expression 
)

# Filter the DDS object
keep <- rowSums(counts(dds_1)) >= 10
dds_1 <- dds_1[keep,]

# Perform DESeq analysis
dds_1 <- estimateSizeFactors(dds_1)
dds_1 <- estimateDispersionsGeneEst(dds_1)
dispersions(dds_1) <- mcols(dds_1)$dispGeneEst
dds_1 <- nbinomWaldTest(dds_1)

# Extract results
result_names_1 <- resultsNames(dds_1)
result_names_1

# Specify the contrast and extract results
res_1 <- results(dds_1)
res_1

# # Specify the date
# date <- '020224'

# Specify the directory path
directory_path <- paste("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/deseq/")

# Check if the directory exists, and create it if not
if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}

# Specify the file path where you want to save the results
result_file <- paste(directory_path, "/primary_zfpexp_results.csv", sep="")

# Write the results table to a CSV file
write.csv(as.data.frame(res_1), file = result_file)

# Save the DESeqDataSet object
saveRDS(dds_1, file = paste(directory_path, "/primary_zfpexp_dds.rds", sep=""))

# Save the DESeq results
saveRDS(results(dds_1), file = paste(directory_path, "/primary_zfpexp_results.rds", sep=""))