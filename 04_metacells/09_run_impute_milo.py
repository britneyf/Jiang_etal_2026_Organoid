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

ifile=sys.argv[1]

import scanpy as sc

adata = sc.read_h5ad(ifile)

#sc.pp.filter_genes(adata, min_cells=1)

#################

'''
IMPUTE
'''
import scanpy.external as sce
adata2 = adata.copy()

adata_magic = sce.pp.magic(adata2, name_list='all_genes', copy = True, knn=5, t=3, n_pca=adata.obsm['X_pca'].shape[1])

adata.layers['imputed'] = adata_magic.X.astype(np.float16)

adata.write_h5ad(ifile)


'''
store = pd.HDFStore(fn) #out_dir + 'imp_df.%s.090219.h5' % mode)
store['imp_merge'] = imp_merge.astype(np.float16)
store.close()
'''
