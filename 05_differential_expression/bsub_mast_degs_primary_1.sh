
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_new/mast/
odir=/data/chanjlab/$project/output_new/mast/
mkdir -p $odir $log_dir

cd $log_dir

bsub -J mast_p1 -n 10 -R span[ptile=10] -R rusage[mem=20] -W 5:00 -o mast.%J.stdout -eo mast.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/run_mast_degs_primary_1.R"