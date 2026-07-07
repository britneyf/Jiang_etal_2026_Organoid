import sys
import os
import pandas as pd
import re
import numpy as np
import glob
from pathlib import Path
from scipy import sparse
from copy import deepcopy
import seaborn as sns
import csv
import itertools
import warnings
import scanpy as sc
import pickle as pkl
import SEACells
import gzip

sample = sys.argv[1]
n_SEACells = int(sys.argv[2])
target_mc = int(sys.argv[3])

mode = 'combined'
ind_dir = f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/input/metacells/{sample}/'
comb_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/postprocess_adata/postprocess_adata.013124/'
out_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/metacells/metacells.020124/'
os.makedirs(out_dir, exist_ok=True)

ind_file = ind_dir + f'adata.{sample}.h5ad'
ad = sc.read_h5ad(ind_file)

# # Read the sample information from the CSV file
# aliases = pd.read_csv('/data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/Organoid_Sample_Description.txt', sep='\t', header=0)

# # Iterate through each sample in the DataFrame
# for index, row in aliases.iterrows():
#     sample = row['Sample']
    
#     ind_dir = f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/input/metacells/{sample}/'

#     ind_file = ind_dir + f'adata.{sample}.h5ad'
#     ad = sc.read_h5ad(ind_file)

# mode='combined'
# comb_dir = '/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.combined.ambient_corrected.042023/'
# ind_dir = '/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.individual.042023/%s/' % sample
# out_dir = '/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.combined.ambient_corrected.042023/'
# os.makedirs(out_dir,exist_ok=True)

# aliases = pd.read_csv('/home/chanj3/scripts/HTA.multiome_plasticity.combined.042023/HTA.multiome_plasticity.manifest.txt', sep='\t',header=None,index_col=0, names=['AWS'])

# ind_file = ind_dir + 'adata.%s.042023.h5ad' % sample
# ad = sc.read_h5ad(ind_file)

#################
# Directly use the 'raw' layer from ad
#ad.X = ad.X.A if sparse.issparse(ad.X) else ad.X

# if ad.raw is not None:
#     # 'raw' layer exists, proceed with the copy
#     ad_mc = sc.AnnData(ad.raw.X.copy())
# else:
#     # 'raw' layer does not exist, handle accordingly
#     print("Warning: 'raw' layer does not exist in the original AnnData object.")
#     ad_mc = sc.AnnData(ad.X.copy())  # Assuming 'X' should be used in this case
    
ad_mc = sc.AnnData(ad.raw.X.copy())   
ad_mc.obs_names = ad.raw.obs_names
ad_mc.obs = ad.obs
ad_mc.var_names = ad.raw.var_names
ad_mc.obsm = ad.obsm
ad_mc.uns = ad.uns

## User defined parameters

## Core parameters 

build_kernel_on = 'X_pca' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

model = SEACells.core.SEACells(ad_mc, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)

model.construct_kernel_matrix()
M = model.kernel_matrix

# Initialize archetypes
model.initialize_archetypes()

model.fit(min_iter=10, max_iter=200)

SEACell_ad = SEACells.core.summarize_by_SEACell(ad_mc, SEACells_label='SEACell', summarize_layer='X')

#odir = f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/metacells_noscvi/{sample}/SEACells/{target_mc}/'
out_dir = f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/metacells/metacells.020124/{sample}/{target_mc}/'
os.makedirs(out_dir, exist_ok=True)
ofile = out_dir + f'metacells_adata.{sample}.{target_mc}.h5ad'
SEACell_ad.write_h5ad(ofile)

# odir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/metacells/{sample}/SEACells/$target_mc/
# #odir=ind_dir + 'SEACells/%d/' % n_SEACells
# os.makedirs(odir,exist_ok=True)
# ofile = odir + 'SEACell_adata.%s.%d.h5ad' % (sample, n_SEACells)
# SEACell_ad.write_h5ad(ofile)

model.save_assignments(out_dir)

model.get_hard_assignments().to_csv(out_dir + 'mc_assignments.csv')

'''
mc_results = {}
mc_results['hard_assignments'] = model.get_hard_assignments()
mc_results['A_'] = model.A_

fp=gzip.open(ind_dir + 'SEACells/SEACell_model.%s.%d.p.gz' % (sample, n_SEACells),'wb')
pkl.dump(mc_results, fp)
fp.close()
'''
