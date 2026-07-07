# By: Britney Forsyth

# Load libraries
suppressMessages(library(MAST))
suppressMessages(library(lme4))
suppressMessages(library(data.table))
suppressMessages(require("Matrix"))
library("reticulate")

# Use the path to your Python executable in the virtual environment
use_python("/lila/home/forsythb/.virtualenvs/r-reticulate/bin/")

# Look at where Python is located
py_config()

# Read the leiden value from command-line argument
args <- commandArgs(trailingOnly = TRUE)
leiden_value <- as.numeric(args[1])

# Import scanpy
sc <- import("scanpy")

# Read in the adata for primary tumor site
adata <- sc$read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/postprocess_adata/adata.combined.postprocess.leiden.h5ad')

# Get the counts matrix
counts <- as.matrix(adata$X)

# Normalize
ms <- rowSums(counts)
norm_df <- counts / ms * median(ms)

# Log transform
counts <- log(norm_df + 1)

# Assign row names and column names
colnames(counts) <- adata$var_names$to_list()
rownames(counts) <- adata$obs_names$to_list()

# Clusters
# Read in the clusters from adata$obs
obs <- "leiden_res_.25"  # Observation column name
obs_df <- as.data.frame(adata$obs)
cell_clusters <- as.character(obs_df[, obs])
names(cell_clusters) <- rownames(obs_df)
cell_clusters <- cell_clusters[rownames(counts)]
# obs <- "leiden_res_0.1"  # Observation column name
# obs_df <- as.data.frame(adata$obs)
# cell_clusters <- as.character(obs_df[, obs])
# names(cell_clusters) <- rownames(obs_df)
# cell_clusters <- cell_clusters[rownames(counts)]

# Function to perform pairwise DE analysis
# Function to perform pairwise DE analysis
pairwise_de <- function(df1, df2, title) {
    # Combine data frames
    o_df <- rbind(df1, df2)
    condition <- rep(c(1, 2), times = c(nrow(df1), nrow(df2)))

    # Prepare data for MAST
    wellKey <- rownames(o_df)
    cdata <- data.frame(wellKey = wellKey, condition = condition)
    fdata <- data.frame(primerid = colnames(o_df))

    # Create single-cell assay object
    sca <- FromMatrix(t(o_df), cdata, fdata)
    #sca <- FromMatrix(o_df, cdata, fdata)
    cdr2 <- colSums(assay(sca) > 0)
    colData(sca)$cngeneson <- scale(cdr2)
    colData(sca)$cond <- factor(unlist(as.list(condition)))
    colData(sca)$cond <- relevel(colData(sca)$cond, 2)

    # Fit model
    zlmCond <- zlm(~ cond + cngeneson, sca)
    summaryDt <- summary(zlmCond, doLRT = 'cond1')$datatable

    # Extract significant results
    fcHurdle <- merge(summaryDt[contrast == 'cond1' & component == 'H', .(primerid, `Pr(>Chisq)`)],
                      summaryDt[contrast == 'cond1' & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)],
                      by = 'primerid')
    fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
    #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
    setorder(fcHurdle, fdr)

    # Export data
    write.csv(as.data.frame(fcHurdle), file.path(out_dir, sprintf("%s.csv", title)), quote = FALSE)
}

# Loop over each cluster and perform pairwise DE analysis
out_dir <- '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/mast_25/'
dir.create(out_dir, showWarnings = FALSE)

df1 <- counts[cell_clusters == leiden_value, ]
df2 <- counts[cell_clusters != leiden_value, ]
pairwise_de(df1, df2, sprintf("%s_%s", obs, leiden_value))

#'''OLD CODE'''

# # Load libraries
# suppressMessages(library(MAST))
# suppressMessages(library(lme4))
# suppressMessages(library(data.table))
# suppressMessages(require("Matrix"))
# library("reticulate")

# # Use the path to your Python executable in the virtual environment
# use_python("/lila/home/forsythb/.virtualenvs/r-reticulate/bin/")

# # Look at where Python is located
# py_config()

# # Import scanpy
# sc <- import("scanpy")

# # Get the leiden value from command line argument
# #leiden_value <- as.numeric(commandArgs(trailingOnly = TRUE))

# # Read in the adata for primary tumor site
# adata_subset = sc$read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/postprocess_adata/adata.combined.postprocess.leiden.h5ad')

# # Subset adata for the specified leiden value
# #adata_subset <- adata[adata$obs$leiden_res_0.1 == leiden_value, , drop = FALSE]
# #adata_subset <- subset(adata, subset = adata$obs$leiden_res_0.1 == leiden_value)

# # Get the counts matrix
# counts = as.matrix(adata_subset$X)

# # Normalize
# ms <- rowSums(counts)
# norm_df <- counts/ms * median(ms)

# # Log transform
# counts <- log(norm_df + 1)

# # Transpose
# counts <- t(counts)

# # Assign row and column names to thte counts matrix
# rownames(counts) <- adata_subset$var_names$to_list()
# colnames(counts) <- adata_subset$obs_names$to_list()

# # Make a single cell assay object
# sca <- FromMatrix(counts, adata_subset$obs, adata_subset$var)

# # Prepare the SCA data
# cdr2 <- colSums(assay(sca)>0)
# colData(sca)$cngeneson <- scale(cdr2)

# # Set covariates 
# zfp_expression <- factor(colData(sca)$ZFP_Expression)
# zfp_expression <- relevel(zfp_expression, 'CTRL')
# # culture_media <- factor(colData(sca)$Culture_Media)
# # culture_media <- relevel(culture_media, 'BASE')
# # tumor_site <- factor(colData(sca)$Tumor_Site)
# # tumor_site <- relevel(tumor_site, 'Primary')
# cluster <- factor(colData(sca)$leiden_res_0.1)

# # Create a design matrix
# zlmCond <- zlm(~ zfp_expression + cluster, sca)

# # Fit the model
# summaryCond <- summary(zlmCond, doLRT='zfp_expressionZFP_KD')
# summaryCond <- summaryCond$datatable

# fcHurdle <- merge(summaryCond[contrast=='zfp_expressionZFP_KD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
#                   summaryCond[contrast=='zfp_expressionZFP_KD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

# fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
# fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>log2(1.5)], as.data.table(mcols(sca)), by='primerid')
# setorder(fcHurdleSig, fdr)

# # Export data
# write.csv(as.data.frame(fcHurdleSig), 
#           sprintf("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/mast/mast_leiden_v3.csv"), 
#           quote = FALSE)

# # Export data
# #write.csv(as.data.frame(fcHurdle), sprintf("%s/%s.csv", '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/mast/mast.022124/', 'primary_1_zfpexp'), quote=FALSE)
