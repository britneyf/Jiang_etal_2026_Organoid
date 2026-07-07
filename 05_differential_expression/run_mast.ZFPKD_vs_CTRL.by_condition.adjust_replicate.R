# By: Britney Forsyth

# Load libraries
suppressMessages(library(MAST))
suppressMessages(library(data.table))
suppressMessages(require("Matrix"))

# Read the leiden value from command-line argument
args <- commandArgs(trailingOnly = TRUE)
sample = args[1]

# Normalized data frame
in_dir = '/home/chanj3/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/stratified'
mtx_file = sprintf('%s/raw_counts/%s/counts.%s.mtx', in_dir, sample, sample)
counts = as.matrix(readMM(mtx_file))

bc_file = sprintf('%s/raw_counts/%s/counts.%s.barcodes.csv', in_dir, sample, sample)
bc = read.csv(bc_file,header=F,col.names = 'Cell', stringsAsFactors = F)

g_file = sprintf('%s/raw_counts/%s/counts.%s.genes.csv', in_dir, sample, sample)
g = read.csv(g_file,header=F,col.names = 'Gene', stringsAsFactors = F)
rownames(counts) = bc$Cell
colnames(counts) = g$Gene

# Normalize
ms <- rowSums(counts)
norm_df <- counts / ms * median(ms)

norm_df = norm_df[,!grepl("^MT-|^MTMR|^MTND|^RPS|^RPL|^MRP|^FAU$|UBA52", colnames(norm_df))]

# Log transform
counts <- log(norm_df + 1)

# Clusters
# Read in the clusters from adata$obs

obs_file = sprintf('%s/obs/obs.%s.csv', in_dir, sample)
obs_df = read.table(obs_file, sep =',', row.names=1, header=T)

obs <- "ZFP_Expression"  # Observation column name
cell_clusters <- as.character(obs_df[, obs])
names(cell_clusters) <- rownames(obs_df)
cell_clusters <- cell_clusters[rownames(counts)]

# obs <- "leiden_res_0.1"  # Observation column name
# obs_df <- as.data.frame(adata$obs)
# cell_clusters <- as.character(obs_df[, obs])
# names(cell_clusters) <- rownames(obs_df)
# cell_clusters <- cell_clusters[rownames(counts)]



out_dir <- '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/stratified/mast.ZFPKD_vs_CTRL.adjust_replicate'
dir.create(out_dir, showWarnings = FALSE)

df1 <- counts[cell_clusters == 'ZFP_KD', ]
df2 <- counts[cell_clusters == 'CTRL', ]

    replicate = factor(as.numeric(grepl('ZFPKD_2',obs_df$Sample)))
    replicate = relevel(replicate, 2)
    
    replicate = c(replicate[cell_clusters == 'ZFP_KD'], replicate[cell_clusters == 'CTRL'])

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
    colData(sca)$replicate <- replicate

    # Fit model
    zlmCond <- zlm(~ cond + replicate + cngeneson, sca)
    summaryDt <- summary(zlmCond, doLRT = 'cond1')$datatable

    # Extract significant results
    fcHurdle <- merge(summaryDt[contrast == 'cond1' & component == 'H', .(primerid, `Pr(>Chisq)`)],
                      summaryDt[contrast == 'cond1' & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)],
                      by = 'primerid')
    fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
    #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
    setorder(fcHurdle, fdr)

    # Export data
    write.csv(as.data.frame(fcHurdle), file.path(out_dir, sprintf("mast.ZFPKD_vs_CTRL.%s.adjust_replicate.042524.csv", sample)), quote = FALSE)




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
