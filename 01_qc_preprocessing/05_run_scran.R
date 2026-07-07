# AUTHOR: 
# Britney T. Forsyth

# DESCRIPTION:
# R script to run SCRAN. Performs normalization of adata using counts matrix, gene, and barcodes file from concatenation.
# Uses a deconvolution approach to partition cells into pools and performs normalization across cells in each pool. 

# Remove all objects from the current workspace (R memory).
rm(list = ls())

# Set a CRAN mirror
#options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required libraries
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }

# Install additional packages
#packages <- c("bitops", "boot", "class", "cluster", "codetools", 
#              "GenomeInfoDb", "KernSmooth", "Matrix", "nlme", "nnet", 
#              "RCurl", "rpart", "spatial", "SummarizedExperiment", 
#              "survival")

#install.packages(packages)

# BiocManager::install("scran")
# BiocManager::install("scRNAseq")
# BiocManager::install("scater")

# Load libraries
suppressMessages(library(scran))
suppressMessages(library(Matrix))
suppressMessages(library(scater))
suppressMessages(library(scRNAseq))

# Print warnings
#warnings()

args <- commandArgs(trailingOnly = TRUE)

# Load DGC matrix file and metadata
mtx_file <- args[1]
g_file <- args[2]
bc_file <- args[3]
#ofile <- args[4]

# Transpose matrix for python - R connection
counts <- t(as.matrix(readMM(mtx_file)))

# Read cell names
bc <- read.csv(bc_file, header=F, col.names = 'Cell', stringsAsFactors = F)

# Read gene names
g <- read.csv(g_file,header=F,col.names = 'Gene', stringsAsFactors = F)

print("Data loaded")

ind = rowSums(counts)!=0
counts = counts[ind,]
g = g[ind,]

# Create single cell experiment
sce <- SingleCellExperiment(assays = list(counts = counts))

# Store the gene names in this object
rownames(sce) <- g
rowData(sce) <- "Gene"

# Store the gene names in this object
colnames(sce) <- bc$Cell

print("SCE created")

# Calculate stats
qcstats <- perCellQCMetrics(sce)
qcfilter <- quickPerCellQC(qcstats)
sce <- sce[,!qcfilter$discard]

print("QC filter ")
#summary(qcfilter$discard)

# Compute clusters
clusters <- quickCluster(sce)

# Compute factors
sce <- computeSumFactors(sce, clusters=clusters)

print("Size factors")
summary(sizeFactors(sce))

# Normalize
sce <- logNormCounts(sce, log=F)
print(sce)

# Save this gene matrix to a tsv file
mat <- t(normcounts(sce))

# Barcodes file
barcode_file <- file("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/scran/scran.barcodes.csv")
writeLines(rownames(mat), barcode_file)
close(barcode_file)

# Genes file
gene_file <- file("/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/scran/scran.genes.csv")
writeLines(colnames(mat), gene_file)
close(gene_file)

# Convert to sparse matrix and write to file
sce_norm <- Matrix(mat, sparse = TRUE)
writeMM(sce_norm, "/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/scran/scran.matrix.mtx")

# End script