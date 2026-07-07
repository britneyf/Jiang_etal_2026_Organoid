# By Britney T Forsyth
# Description: Script to test out different numbers of archetypes with seacells

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
import matplotlib
import matplotlib.pyplot as plt
from scipy.io import mmread
import scanpy as sc
from scipy.stats import zscore
from anndata import concat
import seaborn as sns
from pathlib import Path
import matplotlib.patheffects as PathEffects
from adjustText import adjust_text
from scipy.stats import zscore
import gseapy as gp
from gseapy.plot import dotplot
from py_pcha import PCHA
from sklearn.decomposition import PCA
import umap
from sklearn.neighbors import NearestNeighbors
import SEACells
from tqdm.notebook import tqdm
from tqdm import tqdm

# Functions
# Greedy initialization only
def initialize_archetypes_greedy_only(self):
    """Initialize B matrix which defines cells as SEACells using only greedy selection."""
    k = self.k

    # Use greedy initialization only
    from_greedy = self.k
    greedy_ix = self._get_greedy_centers()
    if self.verbose:
        print(f"Selecting {from_greedy} cells from greedy initialization.")

    all_ix = np.hstack([greedy_ix])

    unique_ix, ind = np.unique(all_ix, return_index=True)
    all_ix = unique_ix[np.argsort(ind)][:k]
    self.archetypes = all_ix

# Refine SEACells class to include the new method
SEACells.core.SEACells.initialize_archetypes_greedy_only = initialize_archetypes_greedy_only

# Load the adata_g
adata_g = sc.read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/Patients.HTAN/output/local_corr.ZFP36L2/adata.local_corr.gene_modules.h5ad')

# Test out different numbers of archetypes
for n_SEACells in range(9, 21):
    print(f"Testing with n_SEACells = {n_SEACells}")

    # Initialize SEACells model
    model = SEACells.core.SEACells(adata_g, 
                      build_kernel_on='X_pca', 
                      n_SEACells=n_SEACells,
                      convergence_epsilon=1e-5)

    # Construct the kernel matrix
    model.construct_kernel_matrix()
    M = model.kernel_matrix

    # Initialize archetypes using the greedy initialization method
    model.initialize_archetypes_greedy_only()

    # Fit the model
    model.fit(min_iter=100, max_iter=5000)

    # Save the kernel matrix and archetypes to files
    np.save(f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/archetypal_analysis/081224/kernel_matrix_n{n_SEACells}.npy', M)
    np.save(f'/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/archetypal_analysis/081224/archetypes_n{n_SEACells}.npy', model.archetypes)
