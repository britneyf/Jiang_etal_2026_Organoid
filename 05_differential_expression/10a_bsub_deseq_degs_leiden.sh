# By: Britney Forsyth

project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_new/deseq/
odir=/data/chanjlab/$project/output_new/deseq/
mkdir -p $odir $log_dir


bsub -J deseq_leiden -n 10 -R span[ptile=10] -R rusage[mem=25] -W 10:00 -o deseq_leiden.%J.stdout -eo deseq_leiden.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/10a_run_deseq_degs_leiden.R"