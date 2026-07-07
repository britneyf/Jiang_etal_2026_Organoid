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

idir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/postprocess_adata/'
out_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/local_corr.ZFP36L2/'
os.makedirs(out_dir, exist_ok=True)

from scipy.io import mmread

import scanpy as sc

ofile = idir +  'adata.combined.postprocess.leiden.h5ad'
adata = sc.read_h5ad(ofile)

sc.pp.neighbors(adata, n_neighbors=30)
knn = adata.obsp['distances']

from scipy.stats import zscore

cell_types_df = pd.read_excel("/data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/cell_type_kg.xlsx", header=None)

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

df = pd.read_csv('/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/output_12.04.23_df/knnDREMI.Tumor_vs_TAISC_pathways.csv', sep = ',', index_col=0)
knndremi_genes = df.index[df.score > 1.5].to_list()

# Combine genes_list and NE_markers
genes_to_add = list(set(genes_list + NE_markers + ZFP36L2_targets + knndremi_genes))


sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    subset=False,
    flavor="seurat",
    batch_key = 'Batch'
)

adata.var.loc[adata.var.index.intersection(genes_to_add),'highly_variable'] = True

bad_genes = adata.var_names.str.contains(
    "^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")

adata.var.loc[bad_genes,'highly_variable'] = False

mat = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs.index).loc[:,adata.var.highly_variable].apply(zscore)

g = mat.columns.difference(['ZFP36L2'])

rho_mat = pd.DataFrame(0, index=adata.obs.index, columns = g)

# Define a function to calculate entropy for a given array of distances
def calculate_entropy(distances):
    return entropy(distances, base=2)

# Create an empty DataFrame to store entropy values
entropy_df = pd.DataFrame(index=adata.obs_names, columns=np.arange(30))

# Iterate over each cell's kNN neighbors
for i in range(knn_distances.shape[0]):
    # Get the indices of nonzero distances
    idx = np.nonzero(knn_distances[i])[0]
    # Extract distances for the kNN neighbors
    neighbor_distances = knn_distances[i, idx]
    # Calculate entropy for the kNN neighbor distances
    neighbor_entropy = calculate_entropy(neighbor_distances)
    # Store the entropy values in the DataFrame
    entropy_df.loc[adata.obs_names[i], idx] = neighbor_entropy

rho_mat.astype(np.float16).to_csv(out_dir + 'local_corr.ZFP36L2.csv', sep = ',')
