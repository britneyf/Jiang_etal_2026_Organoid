'''
AUTHOR: 
Britney T. Forsyth

DESCRIPTION:
Python script to process adata. Performs PCA on normalized and log-transformed data. Calculate nearest neighbors and perform clustering. Compute UMAP embedding based on the PCA results. Compute diffusion map and perform differential expression analysis.
'''
# Import packages
import os
import matplotlib.pyplot as plt
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
import scanpy as sc

# Set directories
mode='combined'
out_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/postprocess_adata/'
os.makedirs(out_dir,exist_ok=True)

sample_df = pd.read_csv('/data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/' + 'organoid_sample_description.txt', sep='\t',index_col=0)

# These are the output files from scran
mtx_fn = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/scran/' +  'scran.matrix.mtx'
g_fn = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/scran/' + 'scran.genes.csv'
bc_fn = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/scran/' + 'scran.barcodes.csv'

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

ifile = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/concatenated/' + 'adata.combined.h5ad'

adata = sc.read_h5ad(ifile)
adata.raw = adata

bc_intersect = adata.obs.index.intersection(norm_df.index)
g_intersect = adata.var_names.intersection(norm_df.columns)
adata = adata[bc_intersect, g_intersect]

norm_df = norm_df.loc[adata.obs.index, adata.var_names]

adata.X = norm_df.values
adata.layers['without_log'] = adata.X

#adata.X = np.log2(adata.X+1)
sc.pp.log1p(adata, base = 2)

# Filter sample_df2 to exclude samples containing '*_shZFP36L2_4'
exclude_pattern = '_shZFP36L2_4'
filtered_samples = [sample for sample in samples.index if exclude_pattern not in sample]
sample_df2 = sample_df.loc[samples[filtered_samples]]

# Assign the modified sample information to AnnData
sample_df2.index = adata.obs.index

adata.obs = pd.concat([adata.obs.drop('Sample', axis=1), sample_df2], axis=1)

#################
    
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA

import numpy.matlib
def kneepoint(vec):
    curve =  [1-x for x in vec]
    nPoints = len(curve)
    allCoord = np.vstack((range(nPoints), curve)).T
    np.array([range(nPoints), curve])
    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(vecFromFirst * numpy.matlib.repmat(lineVecNorm, nPoints, 1), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
    idxOfBestPoint = np.argmax(distToLine)
    return idxOfBestPoint

def RunPCA(cts, var_threshold, n_components=300):
    pca = PCA(n_components=n_components, svd_solver='randomized')
    pca.fit(cts)
    num_components = 0
    num_components = max(num_components,kneepoint(np.cumsum(pca.explained_variance_ratio_)))
    num_components = max(num_components,np.where(np.cumsum(pca.explained_variance_ratio_) > var_threshold)[0][0])
    var_explained = np.cumsum(pca.explained_variance_ratio_)[num_components]
    print('# Components = %d' % (num_components+1))
    print('Variance explained = %f' % var_explained)
    return pca, num_components, var_explained

######################
'''
PCA
'''
# norm_df = pd.DataFrame(adata.X, index=adata.obs_names, columns = adata.var_names)

# bad_genes = norm_df.columns.str.contains(
#     "^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")
# norm_df = norm_df.loc[:,~bad_genes]

# Remove bad genes
bad_genes = adata.var_names.str.contains("^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")
good_genes = ~bad_genes
norm_df = pd.DataFrame(adata.X.tocsr().toarray()[:, good_genes], index=adata.obs_names, columns=adata.var_names[good_genes])
adata = adata[:, good_genes]

print('Performing PCA')
n_components=500
pca = PCA(n_components=n_components, svd_solver='randomized')
pca.fit(norm_df)

#By Kneepoint
num_components = 0
num_components = max(num_components,kneepoint(np.cumsum(pca.explained_variance_ratio_)))
print('# Components = %d' % (num_components+1))

var_explained = np.cumsum(pca.explained_variance_ratio_)[num_components]
print('Variance explained = %f' % var_explained)

pca = PCA(n_components=num_components, svd_solver='randomized')
pca_merge = pd.DataFrame(pca.fit_transform(norm_df.values),
                index=norm_df.index)
adata.obsm['X_pca'] = pca_merge.loc[adata.obs_names,:].values
adata.uns['num_components'] = num_components
adata.uns['var_explained'] = var_explained

######################
'''
NEAREST NEIGHBORS
'''
print('Performing nearest neighbors')
n_neighbors=30
min_dist = 0.3
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=pca_merge.shape[1])

######################
'''
CLUSTERING
'''
print('Phenograph Clustering')
import phenograph
clusters_merge, _, _ = phenograph.cluster(pca_merge, k=30)
clusters_merge = pd.Series(clusters_merge, pca_merge.index)

adata.obs['phenograph'] = clusters_merge.loc[adata.obs_names].astype('str').astype('category')

######################
'''
LEIDEN CLUSTERING
'''
print('Leiden Clustering')
sc.tl.leiden(adata, resolution = 1.8)

######################
'''
UMAP
'''
print('Performing UMAP')
sc.tl.paga(adata, groups = 'phenograph')
sc.pl.paga(adata, plot=False)

sc.tl.umap(adata, init_pos='paga', min_dist=min_dist)


######################
'''
DIFFUSION MAP
'''
print('Performing Diffusion Map')
sc.tl.diffmap(adata, n_comps = 20)

######################
'''
DEG
'''
stats = {}
# Log-transform the count data
#sc.pp.log1p(adata, base=2)

# Perform DEG analysis
print('Performing DEG')
sc.tl.rank_genes_groups(adata, groupby='phenograph', method='wilcoxon')
# stats[group_name, 'logfoldchanges'] = np.log2(
#     pd.concat(stats[group_name, 'scale1'], axis=1)
# )

######################

ofile = out_dir + 'adata.combined.postprocess.h5ad'
adata.write_h5ad(ofile)

adata.obs.to_csv(out_dir + 'obs.combined.postprocess.csv', sep ='\t')

######################
'''
IMPUTE
'''

'''
import scanpy.api as sca
sca.pp.magic(adata, name_list='all_genes', k=30, t=3)
imp_merge = pd.DataFrame(adata.X, index=adata.obs_names, columns = adata.var_names)

store = pd.HDFStore(out_dir + 'imp_df.042023.h5')
store['imp_merge'] = imp_merge.astype(np.float32)
store.close()
'''
