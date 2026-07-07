import os
import scvi
import scanpy as sc
import sys
import numpy as np
import pandas as pd
import pickle
import re

# Suppress the deprecation warning
scvi.settings.dl_pin_memory_gpu_training = False

n_latent = int(sys.argv[1])

ref_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/scripts/organoid_analysis_pipeline_scripts/'
comb_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/concatenated/concatenated.020624/'
out_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/scvi/scvi.021224/scvi.latent_%d/' % n_latent
os.makedirs(out_dir, exist_ok=True)

sample_df = pd.read_csv(ref_dir + 'organoid_sample_description.txt', sep='\t',index_col=0)

adata_fn = comb_dir + 'adata.combined.h5ad'
adata = sc.read_h5ad(adata_fn)

# Extract the matrix (CSR matrix)
mtx_fn = adata.X

# Extract genes and barcodes
g_fn = pd.Index(adata.var_names)
bc_fn = pd.Index(adata.obs_names)

# Use adata.X directly as the count matrix
mtx_fn = adata.X.A

mtx_fn = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/scran/scran.020624/' +  'scran.matrix.mtx'
g_fn = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/scran/scran.020624/' + 'scran.genes.csv'
bc_fn = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/scran/scran.020624/' + 'scran.barcodes.csv'

fn = g_fn
with open(fn, 'r') as f:
    g = [i.strip() for i in f.readlines()]

fn = bc_fn
with open(fn, 'r') as f:
    bc = [i.strip() for i in f.readlines()]

from scipy.io import mmread
norm_df = mmread(mtx_fn)
norm_df = norm_df.tocsr().toarray()
norm_df = pd.DataFrame(norm_df, index=bc, columns=g)

samples = [re.sub('_[ACGT]+-1$','',i) for i in bc]
samples = pd.Series(samples, index = bc)

bc_intersect = adata.obs.index.intersection(norm_df.index)
g_intersect = adata.var_names.intersection(norm_df.columns)
adata = adata[bc_intersect, g_intersect]
adata.layers["counts"] = adata.X.copy()
norm_df = norm_df.loc[adata.obs.index, adata.var_names]

adata.X = norm_df.values
adata.layers['without_log'] = adata.X

#adata.X = np.log2(adata.X+1)
adata.raw = adata  # Save the raw data before log transformation
sc.pp.log1p(adata, base=2)

sample_df2 = sample_df.loc[samples,:]
sample_df2.index = adata.obs.index

adata.obs = pd.concat([adata.obs.drop('Sample', axis=1), sample_df2], axis=1)

#adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated(keep='first')]

adata.obs['Sample'] = adata.obs['Sample'].astype('category')


#####

adata.obs['Sample'] = adata.obs.Batch.astype(str).astype('category')

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000,
    subset=True,
    layer="counts",
    flavor="seurat_v3",  # Change to "seurat"
    batch_key="Sample"
)

scvi.model.SCVI.setup_anndata(
    adata,
    #layer="counts",
    batch_key = 'Sample',
    categorical_covariate_keys = ['Sample'],
    #continuous_covariate_keys=["mito_frac", "RBP_frac"]
    continuous_covariate_keys=["mito_frac"]
)


model = scvi.model.SCVI(adata, n_latent = n_latent)

model

model.train()

model.save(out_dir + "my_model/")

latent = model.get_latent_representation()

adata.obsm["X_scVI"] = latent

denoised = model.get_normalized_expression(adata, library_size=1e4)

adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)
        
adata.write_h5ad(out_dir + 'adata.scvi.latent_%d.h5ad' % n_latent)

