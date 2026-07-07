# By Britney T Forsyth
# Description: Bsub script to run Cell Bender. Requires batch_description.txt file with batch name, path to batch matrix file, and number of expected cells per batch. 

#!/bin/bash
project=CRC_ZFP36L2.092023/Organoid

# Create log directory if it doesn't exist
log_dir=/data/chanjlab/$project/logs_new/cellbender/
mkdir -p "$log_dir"

while read batch_name h5_file expected_cells; do
    bsub -gpu "num=1" -q gpuqueue -J cellbender.$batch_name \
        -n 10 -R span[ptile=10] -R rusage[mem=15] -W 3:00 \
        -o /data/chanjlab/$project/logs_new/cellbender/cellbender.$batch_name.%J.stdout \
        -eo /data/chanjlab/$project/logs_new/cellbender/cellbender.$batch_name.%J.stderr \
        "sh /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/02a_run_cellbender_fororganoid.sh $batch_name $h5_file $expected_cells"

done < /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/batch_description.txt