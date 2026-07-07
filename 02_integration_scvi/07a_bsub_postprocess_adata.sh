#!/bin/bash
'''
AUTHOR: 
Britney T. Forsyth

DESCRIPTION:
Bsub script to process adata. Performs PCA on normalized and log-transformed data. Calculate nearest neighbors and perform clustering. Compute UMAP embedding based on the PCA results. Compute diffusion map and perform differential expression analysis.
'''

# Set up project directory and log directory
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_noreplicates/postprocess_adata/
mkdir -p "$log_dir"
cd $log_dir

# Define bsub parameters
bsub -J postprocess_adata -n 15 -R span[ptile=15] -R rusage[mem=10] -W 3:00 -o postprocess_adata.%J.stdout -eo postprocess_adata.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate doubletdetection; python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/07a_run_postprocess_adata.py"
