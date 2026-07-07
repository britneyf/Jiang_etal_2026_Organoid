# Load libraries
library("anndata")
library("DESeq2")
library("reticulate")
# Load the limma library
library("edgeR")

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
# count_matrix <- adata$X
# count_matrix <- as(count_matrix, "matrix")
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

# Specify the count matrix and column data
count_matrix <- counts(dds_1)  # Use the DESeqDataSet you previously created
coldata <- colData(dds_1)

# Set up the design matrix
design_matrix <- model.matrix(~ Tumor_Site + Culture_Media + ZFP_Expression, data = coldata)

# Create a contrast matrix
contrast_matrix <- makeContrasts(
  ZFP_ExpressionshZFP - ZFP_ExpressionCTRL,
  levels = colnames(design_matrix)
)

# Fit the linear model
fit <- lmFit(count_matrix, design_matrix)

# Apply empirical Bayes smoothing
fit <- eBayes(fit)

# Perform the contrast analysis
contrast_results <- contrasts.fit(fit, contrast = contrast_matrix)
contrast_results <- decideTests(contrast_results)

# Extract and save results
result_table <- topTable(contrast_results, coef = 1)  # coef = 1 corresponds to the specified contrast
date <- '020524'
directory_path <- paste("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/limma/limma.", date, sep = "")

if (!dir.exists(directory_path)) {
  dir.create(directory_path, recursive = TRUE)
}

result_file <- paste(directory_path, "/sc_zfpexp_limma_results.csv", sep = "")
write.csv(result_table, file = result_file)

# Save the limma object
saveRDS(contrast_results, file = paste(directory_path, "/sc_zfpexp_limma_results.rds", sep = ""))
