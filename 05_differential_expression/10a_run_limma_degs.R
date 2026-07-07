# Remove all objects from the current workspace (R memory).
rm(list = ls())

# Load libraries
library("reticulate")
library(edgeR)
library(anndata)

# Use the path to your Python executable in the virtual environment
use_python("/lila/home/forsythb/.virtualenvs/r-reticulate/bin/")

# Look at where Python is located
py_config()

# Import scanpy
sc <- import("scanpy")

# Read in the adata
input_file <- commandArgs(trailingOnly = TRUE)[1]
basename <- commandArgs(trailingOnly = TRUE)[2]
adata <- sc$read_h5ad(input_file)

# Convert AnnData to dataframe
adata_df <- adata$to_df()

# Transpose
df <- t(adata_df)

# Define groups for contrast
group <- factor(paste(adata$obs$ZFP_Expression, adata$obs$Tumor_Site, sep='_'))

# Define the contrast matrix
contrast_matrix <- makeContrasts(
    g.zfp = 'ZFP_KD_Primary - ZFP_KD_Metastatic',
    g.ctrl = 'CTRL_Primary - CTRL_Metastatic',
    levels = levels(group))

# Create the design matrix
design <- model.matrix(~0 + group, data=df)

# Limma fit model
v <- voom(df, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrast_matrix)
efit <- eBayes(vfit)

# Look at the results
de.summary <- summary(decideTests(efit, p.value=0.01, lfc=0))
g_zfp <- topTable(efit, coef="g.zfp", n=Inf, sort.by="p")
g_ctrl <- topTable(efit, coef="g.ctrl", n=Inf, sort.by="p")

# Define output file paths
output_file_g_zfp <- paste0("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/limma/g_zfp_results_", basename, ".csv")
output_file_g_ctrl <- paste0("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/limma/g_ctrl_results_", basename, ".csv")

# Save the results
write.csv(g_zfp, file=output_file_g_zfp, row.names=TRUE)
write.csv(g_ctrl, file=output_file_g_ctrl, row.names=TRUE)
