# Import packages
import os
import pandas as pd
import re
import numpy as np
import glob
from pathlib import Path
from scipy import sparse
from copy import deepcopy
import csv
import itertools
import warnings
import hotspot
import sys
import matplotlib.pyplot as plt
from scipy.io import mmread
import scanpy as sc

# Set the input and output directories
date='082824'
idir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/postprocess_adata/'
hs_dir = f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/hotspot/{date}/'
os.makedirs(hs_dir, exist_ok=True)

#################
'''
LOAD DATA
'''
ofile = idir +  'adata.combined.postprocess.h5ad'
adata = sc.read_h5ad(ofile)

# Subset the data for 'Metastatic' or 'Primary Tumor'
#adata_subset = adata[adata.obs['Tumor_Site'].isin(['Metastatic', 'Primary'])].copy()
#adata_subset = adata[adata.obs['Tumor_Site'].isin(['Metastatic'])].copy()
adata_subset = adata[adata.obs['Tumor_Site'].isin(['Primary'])].copy()

# Filter genes
sc.pp.filter_genes(adata_subset, min_cells=1)

# Load cell type gene information
#cell_types_df = pd.read_excel("/data/chanjlab/CRC_ZFP36L2.092023/ref/GenesSets_MasterSheet.xlsx", header=None)
cell_types_df = pd.read_excel("/data/chanjlab/CRC_ZFP36L2.092023/ref/GenesSets_MasterSheet.xlsx", header=None, sheet_name="HotSpot Modules_Grouped")

# Create a gene dictionary
genes_dict = {}
for index, row in cell_types_df.iterrows():
    cell_type = row.iloc[0]  # Extract the cell type
    genes = row.iloc[2:].dropna().tolist()  # Drop the first two elements from the row and remove NaN values
    genes_dict[cell_type] = genes

# Create a list of genes
genes_list = []
for index, row in cell_types_df.iterrows():
    genes = row.iloc[2:].dropna().tolist()
    genes_list.extend(genes)
genes_list = [x for x in genes_list if x in adata_subset.var_names]

# Additional markers and target genes
NE_markers = ['CHGA','CHGB','SYP','ENO2','NCAM1','INSM1']

biogrid_fn = '/data/chanjlab/CRC_ZFP36L2.092023/ref/BIOGRID.ZFP36L2_targets.txt'
ZFP36L2_targets = pd.read_csv(biogrid_fn, sep='\t').loc[:, 'Official Symbol Interactor B'].tolist()

ppi_fn = '/data/chanjlab/CRC_ZFP36L2.092023/ref/IP_MS_up_146P_146Li_partner_protein.txt'
ZFP36L2_targets2 = pd.read_csv(ppi_fn, index_col=0, sep='\t').index.tolist()

rna_fn = '/data/chanjlab/CRC_ZFP36L2.092023/ref/ZFP36L2_RNA_targets.txt'
ZFP36L2_targets3 = pd.read_csv(rna_fn, index_col=0, sep='\t').index.tolist()

df = pd.read_csv('/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/output_12.04.23_df/knnDREMI.Tumor_vs_TAISC_pathways.csv', sep=',', index_col=0)
knndremi_genes = df.index[df.score > 1.5].tolist()

# Combine all genes
genes_to_add = list(set(genes_list + NE_markers + ZFP36L2_targets + ZFP36L2_targets2 + ZFP36L2_targets3 + knndremi_genes))

# Identify highly variable genes
sc.pp.highly_variable_genes(
    adata_subset,
    n_top_genes=3000,
    subset=False,
    flavor="seurat_v3",
    #layer='raw'
)
# adata_subset.obs['Batch'] = adata_subset.obs['Batch'].astype('category')

# sc.pp.highly_variable_genes(
#     adata_subset,
#     n_top_genes=3000,
#     subset=False,
#     flavor="seurat",
#     batch_key="Batch"
# )

# Mark the genes to keep as highly variable
adata_subset.var.loc[adata_subset.var.index.intersection(genes_to_add), 'highly_variable'] = True

# Exclude bad genes
bad_genes = adata_subset.var_names.str.contains("^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")
adata_subset.var.loc[bad_genes, 'highly_variable'] = False

# Subset the data to keep only highly variable genes
adata_subset = adata_subset[:, adata_subset.var.highly_variable].copy()

# Compute neighbors
sc.pp.neighbors(adata_subset, n_neighbors=50, n_pcs=adata_subset.obsm['X_pca'].shape[1])

# Create a counts layer from raw data
adata_subset.layers['counts'] = adata_subset.raw[adata_subset.obs_names, adata_subset.var_names].X

# Initialize Hotspot
hs = hotspot.Hotspot(
    adata_subset,
    layer_key="counts",
    model='danb',
    distances_obsp_key="distances",
    umi_counts_obs_key="total_counts"
)

# Create KNN graph
hs.create_knn_graph(weighted_graph=False, n_neighbors=15)

# Compute autocorrelations
hs_results = hs.compute_autocorrelations(jobs=30)
hs_results.to_csv(hs_dir + 'hotspot.organoid.primary.grouped.informative_genes.tsv', sep='\t')

# Select genes with significant FDR
hs_genes = hs_results.loc[hs_results.FDR < 0.05].index

# Compute local correlations
local_correlations = hs.compute_local_correlations(hs_genes, jobs=30)
local_correlations.to_csv(hs_dir + 'hotspot.organoid.primary.grouped.local_correlations.tsv', sep='\t')

# Create modules
modules = hs.create_modules(min_gene_threshold=50, core_only=True, fdr_threshold=0.05)
modules.to_csv(hs_dir + 'hotspot.organoid.primary.grouped.modules.tsv', sep='\t')

# Calculate module scores
module_scores = hs.calculate_module_scores()
module_scores.to_csv(hs_dir + 'hotspot.organoid.primary.grouped.module_scores.tsv', sep='\t')

# Plot local correlations
#hs.plot_local_correlations()
#plt.savefig(hs_dir + 'hotspot.organoid.grouped.local_correlation_plot.png')

#########################
