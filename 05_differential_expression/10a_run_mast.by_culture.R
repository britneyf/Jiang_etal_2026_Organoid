# Load libraries
suppressMessages(library(MAST))
suppressMessages(library(lme4))
suppressMessages(library(data.table))
suppressMessages(require("Matrix"))
library("reticulate")
library(stringr)

# Use the path to your Python executable in the virtual environment
use_python("/lila/home/forsythb/.virtualenvs/r-reticulate/bin/")

# Look at where Python is located
py_config()

# Import scanpy
sc <- import("scanpy")

# Read the culture name from command-line argument
args <- commandArgs(trailingOnly = TRUE)
culture = args[1]

# Read in the adata for the sample
in_dir='/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/stratify_by_culture/'
h5ad_file=sprintf('%s/%s.h5ad', in_dir, culture)
adata <- sc$read_h5ad(h5ad_file)

# Set up the counts matrix
counts <- as.matrix(adata$raw$X)
colnames(counts) <- adata$raw$var_names$to_list()
rownames(counts) <- adata$raw$obs_names$to_list()

# Normalize
ms <- rowSums(counts)
norm_df <- counts / ms * median(ms)
norm_df = norm_df[,!grepl("^MT-|^MTMR|^MTND|^RPS|^RPL|^MRP|^FAU$|UBA52", colnames(norm_df))]

# Pseudocount and log-transform
counts <- log(counts + 1)

# Get the observations
obs_df <- as.data.frame(adata$obs)
condition <- factor(obs_df$ZFP_Expression, levels = c('CTRL','ZFP_KD'))
Tumor_Site <- factor(obs_df$Tumor_Site, levels = c('Primary','Metastatic'))
o_df = counts

# Prepare data for MAST
wellKey <- rownames(o_df)
cdata <- data.frame(wellKey = wellKey, condition = condition)
fdata <- data.frame(primerid = colnames(o_df))

# Create single-cell assay object
sca <- FromMatrix(t(o_df), cdata, fdata)
#sca <- FromMatrix(o_df, cdata, fdata)
cdr2 <- colSums(assay(sca) > 0)
colData(sca)$cngeneson <- scale(cdr2)
colData(sca)$cond <- condition
colData(sca)$Tumor_Site <- Tumor_Site

# Make the output directory
out_dir <- '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/mast/mast.by_culture.nonorm/'
dir.create(out_dir, showWarnings = FALSE)

# Fit model
zlmCond <- zlm(~ cngeneson + condition + Tumor_Site + condition:Tumor_Site, sca)
summaryDt <- summary(zlmCond, doLRT = 'conditionZFP_KD:Tumor_SiteMetastatic')$datatable

#Export data
write.csv(as.data.frame(summaryDt), file.path(out_dir, sprintf("mast.by_culture.%s.full_model.csv", culture)), quote = FALSE)