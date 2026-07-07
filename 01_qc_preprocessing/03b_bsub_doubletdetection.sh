# By Britney T Forsyth
# Description: Bsub script to run doublet detection for one sample for testing

project="CRC_ZFP36L2.092023/Organoid"

# Set the base directory for samples
base_dir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/cellbender/"

# Get a list of sample directories from the sample_matrix_and_cellbender_paths.txt file
sample_info_file="/data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/samplematrix_cellbenderpaths_for_doubletdetection.txt"
# Log directory
log_dir="/data/chanjlab/$project/logs_new/doubletdetection"
mkdir -p $log_dir

# Set a fixed output directory for doublet detection
odir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/doubletdetection/"
mkdir -p $odir

# Loop over all samples from the file
while IFS=$'\t' read -r batch_name sample_name sample_matrix_path cellbender_path; do

    # Run doublet detection for the one sample
    bsub -J doublet_detection_$sample_name -n 15 \
         -R span[ptile=15] -R rusage[mem=15] -W 5:00 \
         -o $log_dir/doublet_detection_$sample_name.%J.stdout -eo $log_dir/doublet_detection_$sample_name.%J.stderr \
         "source ~/anaconda3/etc/profile.d/conda.sh; conda activate doubletdetection; \
         python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/03b_run_doubletdetection.py \
         $sample_matrix_path $cellbender_path $odir $sample_name"
done < "$sample_info_file"

