#!/bin/bash
mode=organoid.samples.110223
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_new/postprocess_scvi

# Create log directory if it doesn't exist
mkdir -p "$log_dir"

cd $log_dir

for n_latent in 75 100 125 150; do
    bsub -J postprocess_scvi.$n_latent -n 15 -R span[ptile=15] -R rusage[mem=10] -W 3:00 -o postprocess_scvi.$n_latent.%J.stdout -eo postprocess_scvi.$n_latent.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate doubletdetection; python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/7b_run_postprocess_scvi.py $n_latent"
done

cd -
