suppressMessages(library(MAST))
suppressMessages(library(data.table))
suppressMessages(require("Matrix"))

in_dir = '/home/chanj3/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/stratified/'

args <- commandArgs(trailingOnly = TRUE)
culture = args[1]

sprintf('%s/raw_counts/', in_dir)

mtx_files = list.files(sprintf('%s/raw_counts/', in_dir), recursive=T, pattern=sprintf('counts.*%s.mtx', culture))


counts_list = list()
obs_list = list()


for (mtx_file in mtx_files) {
    sample = dirname(mtx_file)
    
    mtx_file = paste0(in_dir, 'raw_counts/', mtx_file)
    
    counts = as.matrix(readMM(mtx_file)) 
    bc_file = sprintf('%s/raw_counts/%s/counts.%s.barcodes.csv', in_dir, sample, sample)
    bc = read.csv(bc_file,header=F,col.names = 'Cell', stringsAsFactors = F)
    
    g_file = sprintf('%s/raw_counts/%s/counts.%s.genes.csv', in_dir, sample, sample)
    g = read.csv(g_file,header=F,col.names = 'Gene', stringsAsFactors = F)
    rownames(counts) = bc$Cell
    colnames(counts) = g$Gene
    
    counts_list[[sample]] = counts
    
    obs_file = sprintf('%s/obs/obs.%s.csv', in_dir, sample)
    obs_list[[sample]] = read.table(obs_file, sep =',', row.names=1, header=T)
}

counts = do.call(rbind, counts_list)

obs_df = do.call(rbind, obs_list)

ind = obs_df$Culture_Media == culture

counts = counts[ind,]
obs_df = obs_df[ind,]

# Normalize
ms <- rowSums(counts)
norm_df <- counts / ms * median(ms)

norm_df = norm_df[,!grepl("^MT-|^MTMR|^MTND|^RPS|^RPL|^MRP|^FAU$|UBA52", colnames(norm_df))]

# Log transform
counts <- log(norm_df + 1)

# Clusters
# Read in the clusters from adata$obs

condition <- factor(obs_df$ZFP_Expression, levels = c('CTRL','ZFP_KD'))

Tumor_Site <- factor(obs_df$Tumor_Site, levels = c('Primary','Metastatic'))

replicate = grepl('ZFPKD_2', obs_df$Sample)
replicate[grepl('146_M.*ZFPKD_2', obs_df$Sample)] = replicate[grepl('146_M.*ZFPKD_2', obs_df$Sample)] + 1
replicate = factor(replicate)

out_dir <- '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/stratified/mast.ZFPKD_Tissue_Interact.by_culture.adjust_replicate'
dir.create(out_dir, showWarnings = FALSE)

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
colData(sca)$replicate <- replicate


# Fit model
zlmCond <- zlm(~ cngeneson + replicate + condition + Tumor_Site + condition:Tumor_Site, sca)

summaryDt <- summary(zlmCond, doLRT = 'conditionZFP_KD:Tumor_SiteMetastatic')$datatable

#Export data
write.csv(as.data.frame(summaryDt), file.path(out_dir, sprintf("mast.ZFPKD_Tumor_Interact.%s.full_model.csv", culture)), quote = FALSE)