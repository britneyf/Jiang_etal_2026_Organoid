# Import packages
import scanpy as sc
import doubletdetection
import time
import numpy as np
import sys
import os
import re
import pandas as pd
import csv
import joblib
import seaborn as sns
from scipy.sparse import issparse
from cellbender.remove_background.downstream import anndata_from_h5

# Load the data
sample_file = sys.argv[1]  # First argument is sample barcodes file
batch_file = sys.argv[2]   # Second argument is the .combined.h5 batch file
odir = sys.argv[3]         # Third argument is the output directory
sample_name = sys.argv[4]  # Fourth argument is the sample name

# Check if the output directory exists, and create it if it doesn't
if not os.path.exists(odir):
    os.makedirs(odir)

# Set output directory
sc.settings.figdir = odir

# Read sample data
adata_sample = sc.read_10x_mtx(sample_file)

# Extract sample barcodes
adata_sample_barcodes = set(adata_sample.obs_names)

# Read batch data
adata_batch = anndata_from_h5(batch_file)

# Extract batch barcodes
adata_batch_barcodes = set(adata_batch.obs_names)

# Find the common barcodes between sample and batch
common_barcodes = adata_sample_barcodes.intersection(adata_batch_barcodes)

# Filter data based on common barcodes
adata_sample_filtered = adata_sample[adata_sample.obs_names.isin(common_barcodes)]
adata_batch_filtered = adata_batch[adata_batch.obs_names.isin(common_barcodes)]

adata_batch_filtered = adata_batch_filtered[adata_batch_filtered.obs['cell_probability'] > 0.9]
adata_batch_filtered.var_names_make_unique()
sc.pp.filter_genes(adata_batch_filtered, min_cells=1)

# Run doublet detection
clf = doubletdetection.BoostClassifier(
    n_iters=10,
    clustering_algorithm="louvain",
    standard_scaling=True,
    pseudocount=0.1,
    n_jobs=-1,
)

# Calculate doublets
doublets = clf.fit(adata_batch_filtered.X).predict(p_thresh=1e-7, voter_thresh=0.8)
doublet_score = clf.doublet_score()

adata_batch_filtered.obs["doublet"] = doublets
adata_batch_filtered.obs["doublet_score"] = doublet_score

# Save convergence plot with sample name
#f = doubletdetection.plot.convergence(clf, save = odir + sample_name + '.convergence_test.pdf', show=True, p_thresh=1e-7, voter_thresh=0.8)

sc.pp.normalize_total(adata_batch_filtered)
sc.pp.log1p(adata_batch_filtered)
sc.pp.highly_variable_genes(adata_batch_filtered)
sc.tl.pca(adata_batch_filtered)
sc.pp.neighbors(adata_batch_filtered)
sc.tl.umap(adata_batch_filtered)

# Save UMAP plot with sample name
# umap_plot_file = f"{odir}/{sample_name}_umap_doublets.png"
#sc.pl.umap(adata_batch_filtered, color=["doublet", "doublet_score"], save = odir + sample_name + '.umap_doublets.png')

# # Save violin plot with sample name
# violin_plot_file = f"{odir}/{sample_name}_violin_doublets.png"
#sc.pl.violin(adata_batch_filtered, "doublet_score", save = odir + sample_name + '.violin_doublets.png')

# # Save threshold plot with sample name
# threshold_plot_file = f"{odir}/{sample_name}_threshold_test.pdf"
#f3 = doubletdetection.plot.threshold(clf, save = odir + sample_name + '.threshold_test.pdf', show=True, p_step=6)

# # Save doublets file with sample name
doublets_file = f"{odir}{sample_name}.doublets.txt"
adata_batch_filtered.obs.loc[:, ['doublet', 'doublet_score']].to_csv(doublets_file)
