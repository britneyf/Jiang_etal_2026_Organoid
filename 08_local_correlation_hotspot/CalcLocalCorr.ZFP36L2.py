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
import sys
import matplotlib.pyplot as plt
from scipy.io import mmread
import scanpy as sc
from scipy.stats import zscore

# Set the input and output directories
idir = '/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/adatas/'
out_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/calclocalcorr/080524/'
os.makedirs(out_dir, exist_ok=True)

# Choose the adata file
ofile = idir +  'Epithelial.h5ad' # This is the patient data
adata = sc.read_h5ad(ofile)

sc.pp.neighbors(adata, n_neighbors=20)
knn = adata.obsp['distances']

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

biogrid_fn = '/data/chanjlab/CRC_ZFP36L2.092023/ref/BIOGRID.ZFP36L2_targets.txt'
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

# adata.obs['Batch'] = adata.obs['Batch'].astype('category')
# sc.pp.highly_variable_genes(
#     adata,
#     n_top_genes=3000,
#     subset=False,
#     flavor="seurat",
#     batch_key = 'Batch'
# )

adata.var.loc[adata.var.index.intersection(genes_to_add),'highly_variable'] = True

bad_genes = adata.var_names.str.contains(
    "^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")

adata.var.loc[bad_genes,'highly_variable'] = False

mat = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs.index).loc[:,adata.var.highly_variable].apply(zscore)

#mat = pd.DataFrame(adata.X.todense(), columns=adata.var_names, index=adata.obs.index).loc[:,adata.var.highly_variable].apply(zscore)

g = mat.columns.difference(['ZFP36L2'])

rho_mat = pd.DataFrame(0, index=adata.obs.index, columns = g)

for i in np.arange(knn.shape[0]):
    idx = knn[i,:].nonzero()[1]
    mat0 = mat.iloc[[i] + idx.tolist(),:]
    rho_mat.iloc[i,:] = mat0.loc[:,mat0.columns.difference(['ZFP36L2'])].mul(mat0.ZFP36L2, axis=0).sum(axis=0)

rho_mat.astype(np.float16).to_csv(out_dir + 'local_corr.ZFP36L2.seurat.patient.csv', sep = ',')
