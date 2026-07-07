# By: Britney Forsyth

project=CRC_ZFP36L2.092023/Organoid

# Make log directory and output directory
log_dir=/data/chanjlab/$project/logs_new/mast.ZFPKD_vs_CTRL_per_condition.adjust_replicate.042524/
out_dir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/stratified/mast.ZFPKD_vs_CTRL/
mkdir -p $out_dir $log_dir

for sample in `ls /home/chanj3/chanjlab/CRC_ZFP36L2.092023/Organoid/output.042224/stratified/raw_counts`; do
    bsub -sla llSC4 -R A100 -q cpuqueue -J "mast_$sample" -n 15 \
         -R "span[ptile=10]" -R "rusage[mem=20]" -W 24:00 \
         -o "$log_dir/mast.$sample.%J.stdout" -eo "$log_dir/mast.$sample.%J.stderr" \
         "source ~/anaconda3/etc/profile.d/conda.sh; conda activate multiome_env; ~/anaconda3/envs/multiome_env/bin/Rscript --vanilla /home/chanj3/chanjlab/CRC_ZFP36L2.092023/Organoid/scripts/organoid_analysis_pipeline_scripts/run_mast.ZFPKD_vs_CTRL.by_condition.adjust_replicate.R $sample"
done


# # Loop through unique leiden values
# bsub -J mast_leiden -n 15 -R span[ptile=10] -R rusage[mem=30] -W 10:00 -o mast_leiden.%J.stdout -eo mast_leiden.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/run_mast_by_leiden_cluster.R"
