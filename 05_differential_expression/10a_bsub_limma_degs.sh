#!/bin/bash
mode=organoid.samples.110223
project=CRC_ZFP36L2.092023/Organoid

# Create log directory if it doesn't exist
log_dir=/data/chanjlab/$project/logs.042224/limma/
mkdir -p "$log_dir"

# Create output directory if it doesn't exist
odir=/data/chanjlab/$project/output.042224/limma/
mkdir -p "$odir"

# Loop through files in input directory
for input_file in /data/chanjlab/$project/output.042224/limma/input/*.h5ad; do
    # Extract file name without extension
    filename=$(basename -- "$input_file")
    filename_noext="${filename%.*}"
    
    # Run the limma_deg script for each input file
    bsub -gpu "num=1" -q gpuqueue -J "limma.$filename_noext" \
        -n 10 -R span[ptile=10] -R rusage[mem=15] -W 3:00 \
        -o "$log_dir/limma.$filename_noext.%J.stdout" \
        -eo "$log_dir/limma.$filename_noext.%J.stderr" \
        "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/10a_run_limma_degs.R '$input_file' '$filename_noext'"
done
