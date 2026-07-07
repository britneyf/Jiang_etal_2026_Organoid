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

mtx_files = list.files(sprintf('%s/raw_counts/*%s/counts.*.mtx', in_dir, sample))

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

condition <- factor(obs_df$ZFP_Expression, levels = c('CTRL','ZFP_KD'))


out_dir <- '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/stratified/mast.ZFPKD_vs_CTRL.adjust_replicate.042524'
dir.create(out_dir, showWarnings = FALSE)

o_df = counts

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
#    fcHurdle <- merge(summaryDt[contrast == 'cond1' & component == 'H', .(primerid, `Pr(>Chisq)`)],
#                      summaryDt[contrast == 'cond1' & component == 'logFC', .(primerid, coef, ci.hi, ci.lo)],
#                      by = 'primerid')
#    fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
#    #fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
#    setorder(fcHurdle, fdr)

    # Export data
#    write.csv(as.data.frame(fcHurdle), file.path(out_dir, sprintf("mast.ZFPKD_vs_CTRL.Metastatic_vs_Primary.%s.adjust_replicate.csv", sample)), quote = FALSE)


write.csv(as.data.frame(fcHurdle), file.path(out_dir, sprintf("mast.ZFPKD_vs_CTRL.Metastatic_vs_Primary.%s.full_model.csv", sample)), quote = FALSE)
