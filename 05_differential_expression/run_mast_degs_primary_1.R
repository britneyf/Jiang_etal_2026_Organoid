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

# Read in the adata for primary tumor site
adata_primary = sc$read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/postprocess_adata/adata.primary.hvg.h5ad')

# Get the counts matrix
counts = as.matrix(adata_primary$X)

# Normalize
ms <- rowSums(counts)
norm_df <- counts/ms * median(ms)

# Log transform
counts <- log(norm_df + 1)

# Transpose
counts <- t(counts)

# Assign row and column names to thte counts matrix
rownames(counts) <- adata_primary$var_names$to_list()
colnames(counts) <- adata_primary$obs_names$to_list()

# Make a single cell assay object
sca <- FromMatrix(counts, adata_primary$obs, adata_primary$var)

# Prepare the SCA data
cdr2 <- colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)

# Set covariates 
zfp_expression <- factor(colData(sca)$ZFP_Expression)
zfp_expression <- relevel(zfp_expression, 'CTRL')
culture_media <- factor(colData(sca)$Culture_Media)
culture_media <- relevel(culture_media, 'BASE')

# Create a design matrix
zlmCond <- zlm(~ zfp_expression + culture_media, sca)

# Fit the model
summaryCond <- summary(zlmCond, doLRT='zfp_expressionZFP_KD')
summaryCond <- summaryCond$datatable

fcHurdle <- merge(summaryCond[contrast=='zfp_expressionZFP_KD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryCond[contrast=='zfp_expressionZFP_KD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>log2(1.5)], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)

# Export data
write.csv(as.data.frame(fcHurdleSig), sprintf("%s/%s.csv", '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/mast/', 'primary_1_zfpexp'), quote=FALSE)

# Export data
#write.csv(as.data.frame(fcHurdle), sprintf("%s/%s.csv", '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/mast/mast.022124/', 'primary_1_zfpexp'), quote=FALSE)
