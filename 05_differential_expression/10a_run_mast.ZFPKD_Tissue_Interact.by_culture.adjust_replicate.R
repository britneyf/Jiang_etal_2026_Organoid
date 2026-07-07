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

# Import scanpy
sc <- import("scanpy")

# Read the sample name from command-line argument
args <- commandArgs(trailingOnly = TRUE)
sample = args[1]

# Read in the adata for the sample
in_dir='/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/stratify/'
h5ad_file=sprintf('%s/%s.h5ad', in_dir, sample)
adata <- sc$read_h5ad(h5ad_file)

counts <- as.matrix(adata$raw$X)
colnames(counts) <- adata$raw$var_names$to_list()
rownames(counts) <- adata$raw$obs_names$to_list()

# # Normalized data frame
# in_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/stratify/'
# mtx_file = sprintf('%s/matrix/%s.mtx', in_dir, sample)
# counts = as.matrix(mtx_file)

# bc_file = sprintf('%s/barcodes/%s_bc.csv', in_dir, sample)
# bc = read.csv(bc_file,header=F,col.names = 'Cell', stringsAsFactors = F)

# g_file = sprintf('%s/genes/%s_genes.csv', in_dir, sample)
# g = read.csv(g_file,header=F,col.names = 'Gene', stringsAsFactors = F)
# rownames(counts) = bc$Cell
# colnames(counts) = g$Gene

# Normalize
ms <- rowSums(counts)
norm_df <- counts / ms * median(ms)
norm_df = norm_df[,!grepl("^MT-|^MTMR|^MTND|^RPS|^RPL|^MRP|^FAU$|UBA52", colnames(norm_df))]

# Pseudocount and log-transform
counts <- log(counts + 1)

# Clusters
# Read in the clusters from adata$obs

# obs_file = sprintf('%s/obs_df/%s_obs.csv', in_dir, sample)
# obs_df = read.table(obs_file, sep =',', row.names=1, header=T)
obs_df <- as.data.frame(adata$obs)
obs <- "ZFP_Expression"  # Observation column name
cell_clusters <- as.character(obs_df[, obs])
names(cell_clusters) <- rownames(obs_df)
cell_clusters <- cell_clusters[rownames(counts)]
# obs <- "leiden_res_0.1"  # Observation column name
# obs_df <- as.data.frame(adata$obs)
# cell_clusters <- as.character(obs_df[, obs])
# names(cell_clusters) <- rownames(obs_df)
# cell_clusters <- cell_clusters[rownames(counts)]

# Function to perform pairwise DE analysis
pairwise_de <- function(df1, df2) {
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
    setorder(fcHurdle, fdr)

    # Export data
    write.csv(as.data.frame(fcHurdle), file.path(out_dir, sprintf("mast.%s.csv", sample)), quote = FALSE)
}
# # Function to perform pairwise DE analysis
# pairwise_de <- function(df1, df2) {
#     # Combine data frames
#     o_df <- rbind(df1, df2)
#     condition <- rep(c(1, 2), times = c(nrow(df1), nrow(df2)))

#     # Prepare data for MAST
#     wellKey <- rownames(o_df)
#     cdata <- data.frame(wellKey = wellKey, condition = condition)
#     fdata <- data.frame(primerid = colnames(o_df))

#     # Create single-cell assay object
#     sca <- FromMatrix(t(o_df), cdata, fdata)
#     #sca <- FromMatrix(o_df, cdata, fdata)
#     cdr2 <- colSums(assay(sca) > 0)
#     colData(sca)$cngeneson <- scale(cdr2)
#     colData(sca)$cond <- factor(unlist(as.list(condition)))
#     colData(sca)$cond <- relevel(colData(sca)$cond, 2)

#     # Fit model
#     zlmCond <- zlm(~ cond + cngeneson, sca)
#     summaryDt <- summary(zlmCond, doLRT = 'cond1')$datatable

#     # Extract significant results
#     fcHurdle <- merge(summaryDt[contrast == 'cond1' & component == 'H', .(primerid, `Pr(>Chisq)`)],
#                       summaryDt[contrast == 'cond1' & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)],
#                       by = 'primerid')
#     fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
#     #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
#     setorder(fcHurdle, fdr)

#     # Export data
#     write.csv(as.data.frame(fcHurdle), file.path(out_dir, sprintf("mast.ZFPKD_vs_CTRL.%s.csv", sample)), quote = FALSE)
# }

# Loop over each cluster and perform pairwise DE analysis
out_dir <- '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/mast/mast_nonormalization/'
dir.create(out_dir, showWarnings = FALSE)

df1 <- counts[cell_clusters == 'ZFP_KD', ]
df2 <- counts[cell_clusters == 'CTRL', ]
pairwise_de(df1, df2)