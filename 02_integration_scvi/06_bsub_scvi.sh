#!/bin/bash
date=021224
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_new/scvi/scvi.$date/

# Create log directory if it doesn't exist
mkdir -p $log_dir

cd $log_dir

for n_latent in 75 100 125 150; do
    out_dir=/data/chanjlab/$project/output_new/scvi/scvi.$date/scvi.latent_$n_latent
    mkdir -p $out_dir

    #rm -rf $out_dir/my_model
    bsub -J scvi.$n_latent -q gpuqueue -gpu "num=1" -n 10 -R span[ptile=10] -R rusage[mem=12] -W 6:00 -o scvi.$n_latent.%J.stdout -eo scvi.$n_latent.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scvi; /home/forsythb/anaconda3/envs/scvi/bin/python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/6_run_scvi.py $n_latent"

done

cd -
