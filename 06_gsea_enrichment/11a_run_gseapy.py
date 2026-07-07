# Import packages
import gseapy as gp
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
import matplotlib.pyplot as plt
from scipy.io import mmread
import scanpy as sc
import seaborn as sns

# Set the directories
ref_dir = '/data/chanjlab/CRC_ZFP36L2.092023/ref/'
gmt_file = 'pathways_for_gsea.curated.small.042624.gmt'
mast_dir = '/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/mast/mast_nonormalization/'

# Take in the MAST results to create a ranked gene list
mast_file = sys.argv[1]
deg_df=pd.read_csv(mast_dir + mast_file,  sep=',', index_col=1)

# Check if 'Unnamed: 0' column exists and drop it if it does
if 'Unnamed: 0' in deg_df.columns:
    deg_df = deg_df.drop(columns=['Unnamed: 0'])

# Create a ranked genes list
ind = deg_df.loc[:,'Pr(>Chisq)'] == 0
min_pval = deg_df.loc[~ind,'Pr(>Chisq)'].min()
deg_df.loc[ind,'Pr(>Chisq)'] = min_pval

rnk = -np.log2(deg_df.loc[:,'Pr(>Chisq)']) * np.sign(deg_df.coef)

# Run pre-ranked GSEA
pre_res = gp.prerank(rnk=rnk,
                     gene_sets=ref_dir + gmt_file,
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )

# Create the output directory if it doesn't exist
odir='/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/gseapy/'
if not os.path.exists(odir):
    os.makedirs(odir)
    
# Remove the '.csv' extension from mast_file
if mast_file.endswith('.csv'):
    mast_file = mast_file[:-4]

# Construct the output file name
output_file = os.path.join(odir, mast_file + '.gsea.csv')

# Save the DataFrame to CSV
pre_res.res2d.to_csv(output_file)