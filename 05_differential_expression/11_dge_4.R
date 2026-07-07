# Set the CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# If not already downloaded in environment, then download DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")

# Install the R anndata package
install.packages("anndata")

# Run this only if you do not already have an installation of miniconda
#reticulate::install_miniconda()

# Install the Python anndata package
anndata::install_anndata()

# Load libraries
library(DESeq2)
library(anndata)

# Read h5ad file
adata <- read_h5ad("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/metacells_noscvi/adata.combined.postprocessing.h5ad")

# Extract count matrix from the raw layer
count_matrix <- adata$raw$X
count_matrix <- t(count_matrix)
count_matrix <- round(count_matrix)

# Define the column data
coldata <- adata$obs[,c("Tumor_Type","Culture_Media", "ZFP_Expression", "Batch", "Patient")]

# Convert column data to factors
coldata$Tumor_Site <- factor(coldata$Tumor_Type)
coldata$Culture_Media <- factor(coldata$Culture_Media)
coldata$ZFP_Expression <- factor(coldata$ZFP_Expression)
coldata$Batch <- factor(coldata$Batch)
coldata$Patient <- factor(coldata$Patient)

# Design a DESeqDataSet
dds_1 <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ Tumor_Type + Culture_Media + Patient + ZFP_Expression
)

# Add a pseudo-count of 1 to all counts
count_matrix <- counts(dds_1) + 1

# Specify the row names as the gene names
rownames(count_matrix) <- adata$raw$var_names

# Design a DESeqDataSet with modified count matrix
dds_1 <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ Tumor_Site + Culture_Media + Patient + ZFP_Expression:Tumor_Site 
)

# Perform DESeq analysis
dds_1 <- DESeq(dds_1)

# Extract results
result_names_1 <- resultsNames(dds_1)
result_names_1

# Specify the contrast and extract results
res_1 <- results(dds_1)
res_1

# Specify the file path where you want to save the results
result_file <- "/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/deseq/zfpexprint_zfpkd_vs_ctrl_results.csv"

# Write the results table to a CSV file
write.csv(as.data.frame(res_1), file = result_file)

# Save the DESeqDataSet object
saveRDS(dds_1, file = "/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/deseq/zfpexp_inttumor_zfpkd_vs_ctrl_dds.rds")

# Save the DESeq results
saveRDS(results(dds_1), file = "/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/deseq/zfpexp_inttumor_zfpkd_vs_ctrl_results.rds")