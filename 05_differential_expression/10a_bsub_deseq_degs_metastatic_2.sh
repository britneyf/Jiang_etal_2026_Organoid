date=022124
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_new/deseq/
odir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/deseq/
mkdir -p $odir $log_dir

cd $log_dir

bsub -J deseq_m2 -n 10 -R span[ptile=10] -R rusage[mem=20] -W 5:00 -o deseq_m2.%J.stdout -eo deseq_m2.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/10a_run_deseq_degs_metastatic_2.R"