#!/bin/bash
mode=organoid.samples.110223
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs/cellbender

# Create log directory if it doesn't exist
mkdir -p "$log_dir"

while read batch_name h5_file expected_cells; do
    bsub -gpu "num=1" -q gpuqueue -J cellbender.$batch_name \
        -n 10 -R span[ptile=10] -R rusage[mem=15] -W 3:00 \
        -o /data/chanjlab/$project/logs/cellbender/cellbender.$batch_name.%J.stdout \
        -eo /data/chanjlab/$project/logs/cellbender/cellbender.$batch_name.%J.stderr \
        "sh /data/chanjlab/forsythb/britney_run_cellbender_copy.sh $batch_name $h5_file $expected_cells"

done < /data/chanjlab/CRC_ZFP36L2.092023/Organoid/batch_description.txt