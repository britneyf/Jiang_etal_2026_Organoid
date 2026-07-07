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

mode='combined'
comb_dir = '/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.combined.ambient_corrected.042023/'
ind_dir = '/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.individual.042023/%s/' % sample
out_dir = '/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.combined.ambient_corrected.042023/'
os.makedirs(out_dir,exist_ok=True)

aliases = pd.read_csv('/home/chanj3/scripts/HTA.multiome_plasticity.combined.042023/HTA.multiome_plasticity.manifest.txt', sep='\t',header=None,index_col=0, names=['AWS'])

ind_file = ind_dir + 'adata.%s.042023.h5ad' % sample
ad = sc.read_h5ad(ind_file)

#################

## User defined parameters

## Core parameters 

build_kernel_on = 'X_pca' # key in ad.obsm to use for computing metacells
                          # This would be replaced by 'X_svd' for ATAC data

## Additional parameters
n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

model = SEACells.core.SEACells(ad, 
                  build_kernel_on=build_kernel_on, 
                  n_SEACells=n_SEACells, 
                  n_waypoint_eigs=n_waypoint_eigs,
                  convergence_epsilon = 1e-5)

model.construct_kernel_matrix()
M = model.kernel_matrix

# Initialize archetypes
model.initialize_archetypes()

model.fit(min_iter=10, max_iter=200)

SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')

odir=ind_dir + 'SEACells/%d/' % n_SEACells
os.makedirs(odir,exist_ok=True)
ofile = odir + 'SEACell_adata.%s.%d.h5ad' % (sample, n_SEACells)
SEACell_ad.write_h5ad(ofile)

model.save_assignments(odir)

model.get_hard_assignments().to_csv(odir + 'mc_assignments.csv')

'''
mc_results = {}
mc_results['hard_assignments'] = model.get_hard_assignments()
mc_results['A_'] = model.A_

fp=gzip.open(ind_dir + 'SEACells/SEACell_model.%s.%d.p.gz' % (sample, n_SEACells),'wb')
pkl.dump(mc_results, fp)
fp.close()
'''
