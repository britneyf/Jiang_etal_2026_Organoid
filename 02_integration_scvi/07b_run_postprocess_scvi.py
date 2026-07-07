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

n_latent = int(sys.argv[1])

mode='combined'
out_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/postprocess_scvi_test/scvi.latent_%d/' % n_latent
os.makedirs(out_dir,exist_ok=True)

comb_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/scvi/'

ifile = comb_dir + 'scvi.latent_%d/adata.scvi.latent_%d.h5ad' % (n_latent, n_latent)
adata = sc.read_h5ad(ifile)

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
norm_df = pd.DataFrame(adata.X, index=adata.obs_names, columns = adata.var_names)

bad_genes = norm_df.columns.str.contains(
    "^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")
norm_df = norm_df.loc[:,~bad_genes]

# # Remove bad genes
# bad_genes = adata.var_names.str.contains("^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^MRP|^FAU$|UBA52|MALAT")
# good_genes = ~bad_genes
# norm_df = pd.DataFrame(adata.X.tocsr().toarray()[:, good_genes], index=adata.obs_names, columns=adata.var_names[good_genes])
# adata = adata[:, good_genes]

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
print('Performing DEG')
sc.tl.rank_genes_groups(adata, groupby='phenograph', method='wilcoxon')

######################

ofile = out_dir + 'adata.scvi.latent_%d.h5ad' % n_latent
adata.write_h5ad(ofile)

adata.obs.to_csv(out_dir + 'obs.scvi.latent_%d.042023.csv' % n_latent, sep ='\t')

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
