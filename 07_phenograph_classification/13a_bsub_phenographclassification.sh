#!/bin/bash

# By Britney T Forsyth
# Description: Bsub script to run phenograph classification with variable MNN and KNN.

project="CRC_ZFP36L2.092023/Organoid"
date="091624"

# Directory setup
log_dir="/data/chanjlab/$project/logs_noreplicates/harmony_phenograph/$date/cell_state"
mkdir -p $log_dir
odir="/data/chanjlab/$project/output_noreplicates/harmony_phenograph/$date/cell_state"
mkdir -p $odir

# Change to the directory with .h5ad files
cd /data/chanjlab/$project/output_noreplicates/stratify/
unique_samples=($(ls *.h5ad | sed -E 's/\.h5ad$//'))

# Parameter arrays
knn_values=(10 15 20)
mnn_values=(20 30 40)

# Job submission loop
for knn in "${knn_values[@]}"; do
    for mnn in "${mnn_values[@]}"; do
        for sample in "${unique_samples[@]}"; do
            bsub -J "classification_${sample}_knn${knn}_mnn${mnn}" -n 15 \
                 -R "span[ptile=15]" -R "rusage[mem=15]" -W 5:00 \
                 -o "$log_dir/classification_${sample}_knn${knn}_mnn${mnn}.%J.stdout" \
                 -eo "$log_dir/classification_${sample}_knn${knn}_mnn${mnn}.%J.stderr" \
                 "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/13a_run_phenographclassification.py $sample --knn $knn --mnn $mnn"
        done
    done
done

# #!/bin/bash
# # By Britney T Forsyth
# # Description: Bsub script to run phenograph classification. 

# project="CRC_ZFP36L2.092023/Organoid"
# date="091324"

# # Log directory
# log_dir="/data/chanjlab/$project/logs_noreplicates/harmony_phenograph/$date/cell_state_alt1"
# mkdir -p $log_dir

# # Output directory
# odir="/data/chanjlab/$project/output_noreplicates/harmony_phenograph/$date/cell_state_alt1"
# mkdir -p $odir

# # Activate conda environment
# # source ~/anaconda3/etc/profile.d/conda.sh
# # conda activate scrna_clone

# # Change directory to the location of .h5ad files
# cd /data/chanjlab/$project/output_noreplicates/stratify/

# # Extract sample names from the .h5ad files
# unique_samples=($(ls *.h5ad | sed -E 's/\.h5ad$//'))

# # Submit jobs for each sample
# for sample in "${unique_samples[@]}"; do
#     bsub -J "classification_$sample" -n 15 \
#          -R "span[ptile=15]" -R "rusage[mem=15]" -W 5:00 \
#          -o "$log_dir/classification_$sample.%J.stdout" -eo "$log_dir/classification_$sample.%J.stderr" \
#          "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/13a_run_phenographclassification.py $sample"
# done
