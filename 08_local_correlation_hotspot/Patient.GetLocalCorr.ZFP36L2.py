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

idir = '/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/outcombined.10.06.2023/updated_analyses_adata/'
idir = '/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/adatas/'
hs_dir = '/home/chanj3/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/output/hotspot/hvg_1000_within_batch.MT_excluded.orig_pca.042524/'
os.makedirs(hs_dir, exist_ok=True)

#################
'''
LOAD DATA
'''

from scipy.io import mmread

import scanpy as sc

ofile = idir +  'Epithelial.h5ad'

adata = sc.read_h5ad(ofile)
bc = pd.read_csv('/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/outcombined.10.06.2023/updated_analyses_adata/obs.epithelial.csv', sep =',', index_col=0).index

adata = adata[bc,:].copy()

sc.pp.filter_genes(adata, min_cells=1)

cell_types_df = pd.read_excel("/data/chanjlab/CRC_ZFP36L2.092023/ref/cell_type_kg.xlsx", header=None)



# Create a gene dictionary
genes_dict = {}

# Iterate over the generator object and extract genes from each row
for index, row in cell_types_df.iterrows():
    cell_type = row.iloc[0]  # Extract the cell type
    # Drop the first two elements from the row and remove NaN values
    genes = row.iloc[2:].dropna().tolist()
    genes_dict[cell_type] = genes
    
# Define an empty list to store genes
genes_list = []

# Iterate over the generator object and extract genes from each row
for index, row in cell_types_df.iterrows():
    # Drop the first two elements from the row and remove NaN values
    genes = row.iloc[2:].dropna().tolist()
    # Extend the genes list with the extracted genes
    genes_list.extend(genes)

#print(genes_list)
genes_list = [x for x in genes_list if x in adata.var_names]
    
NE_markers = ['CHGA','CHGB','SYP','ENO2','NCAM1','INSM1']

biogrid_fn = '~/chanjlab/CRC_ZFP36L2.092023/ref/BIOGRID.ZFP36L2_targets.txt'
ZFP36L2_targets = pd.read_csv(biogrid_fn, sep = '\t').loc[:,'Official Symbol Interactor B'].to_list()

ppi_fn = '/data/chanjlab/CRC_ZFP36L2.092023/ref/IP_MS_up_146P_146Li_partner_protein.txt'
ZFP36L2_targets2 = pd.read_csv(ppi_fn, index_col=0, sep = '\t').index.to_list()

rna_fn = '/data/chanjlab/CRC_ZFP36L2.092023/ref/ZFP36L2_RNA_targets.txt'
ZFP36L2_targets3 = pd.read_csv(rna_fn, index_col=0, sep = '\t').index.to_list()

df = pd.read_csv('/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/output_12.04.23_df/knnDREMI.Tumor_vs_TAISC_pathways.csv', sep = ',', index_col=0)
knndremi_genes = df.index[df.score > 1.5].to_list()

# Combine genes_list and NE_markers
genes_to_add = list(set(genes_list + NE_markers + ZFP36L2_targets + ZFP36L2_targets2 + ZFP36L2_targets3 + knndremi_genes))

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1000,
    subset=False,
    flavor="seurat_v3",
    layer = 'raw'
)

adata.var.loc[adata.var.index.intersection(genes_to_add),'highly_variable'] = True

bad_genes = adata.var_names.str.contains(
    "^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")

adata.var.loc[bad_genes,'highly_variable'] = False

adata = adata[:,adata.var.highly_variable]

adata = adata.copy()

sc.pp.neighbors(adata, n_neighbors=50, n_pcs=adata.obsm['X_pca'].shape[1])

adata.layers['counts'] = adata.raw[adata.obs_names, adata.var_names].X

hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
#    latent_obsm_key="X_pca",
    distances_obsp_key="distances",
    umi_counts_obs_key="total_counts"
)


hs.create_knn_graph(weighted_graph=False, n_neighbors = 15)

hs_results = hs.compute_autocorrelations(jobs=30)
hs_results.to_csv(hs_dir + 'hotspot.informative_genes.042524.tsv',sep='\t')

hs_genes = hs_results.loc[hs_results.FDR < 0.05].index # Select genes

local_correlations = hs.compute_local_correlations(hs_genes, jobs=30) # jobs for parallelization

local_correlations.to_csv(hs_dir + 'hotspot.local_correlations.042524.tsv',sep='\t')

modules = hs.create_modules(min_gene_threshold=50, core_only=True, fdr_threshold=0.05)

modules.to_csv(hs_dir + 'hotspot.modules.042524.tsv',sep='\t')

module_scores = hs.calculate_module_scores()

module_scores.to_csv(hs_dir + 'hotspot.module_scores.042524.tsv',sep='\t')

hs.plot_local_correlations()

plt.savefig(hs_dir + 'hotspot.local_correlation_plot.042524.png')
