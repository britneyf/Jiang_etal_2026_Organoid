#!/bin/bash
# By Britney T Forsyth
# Description: Bsub script to test out different numbers of archetypes

project="CRC_ZFP36L2.092023/Organoid"
date="081224"

# Log directory
log_dir="/data/chanjlab/$project/logs_noreplicates/archetypal_analysis/$date"
mkdir -p $log_dir

# Output directory
odir="/data/chanjlab/$project/output_noreplicates/archetypal_analysis/$date"
mkdir -p $odir

# Define bsub parameters
bsub -J archetypal_analysis -n 15 -R span[ptile=15] -R rusage[mem=10] -W 3:00 -o archetypal_analysis.%J.stdout -eo archetypal_analysis.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/14_run_archetypal_analysis.py"