# Import packages
import numpy as np
import pandas as pd
import scanpy as sc
import re
import os
import sys
from scipy.sparse import csr_matrix
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import matplotlib

# Read in adata
adata=sc.read_h5ad('/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/postprocess_adata/adata.combined.postprocess.imputed.h5ad')

# Create a dictionary to store the stratified datasets
stratified_data = {}

# Get unique values for each condition
tumor_sites = adata.obs['Tumor_Site'].unique()
culture_media = adata.obs['Culture_Media'].unique()
zfp_expressions = adata.obs['ZFP_Expression'].unique()

# Create the output directory
output_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/stratify_with_zfpkd/'
os.makedirs(output_dir, exist_ok=True)

# Iterate over each combination of conditions
for site in tumor_sites:
    for media in culture_media:
        for expression in zfp_expressions:
            # Generate boolean masks for each condition
            mask_site = adata.obs['Tumor_Site'] == site
            mask_media = adata.obs['Culture_Media'] == media
            mask_expression = adata.obs['ZFP_Expression'] == expression
            
            # Combine the masks using logical AND to get the final mask
            final_mask = mask_site & mask_media & mask_expression
            
            # Apply the mask to subset the original AnnData object
            subset_adata = adata[final_mask, :].copy()
            
            # Recompute PCA for the subset
            sc.pp.pca(subset_adata)
            
            # Add the subset to the dictionary with a unique key
            key = f'{site}_{media}_{expression}'
            stratified_data[key] = subset_adata
            
            # Save the subset to a file
            output_file = os.path.join(output_dir, f'{key}.h5ad')
            subset_adata.write(output_file)
