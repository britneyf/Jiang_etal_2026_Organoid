'''
DESCRIPTION:
Bsub script to impute data.
'''

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
import scanpy as sc

# Read arguments and adata
ifile=sys.argv[1]
basename=sys.argv[2]
adata = sc.read_h5ad(ifile)


#sc.pp.filter_genes(adata, min_cells=1)
#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)

#################

'''
IMPUTE
'''
import scanpy.external as sce
adata2 = adata.copy()

adata_magic = sce.pp.magic(adata2, name_list='all_genes', copy = True, knn=5, t=3, n_pca=adata.obsm['X_pca'].shape[1])

adata.layers['imputed'] = adata_magic.X.astype(np.float16)

output_path = f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/impute/{basename}.impute.h5ad'
adata.write_h5ad(output_path)

