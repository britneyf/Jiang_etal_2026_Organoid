'''
AUTHOR: 
Britney T. Forsyth

DESCRIPTION:
Bsub script to hard stratify the full adata into unique tumor and culture media conditions. Recompute PCA, clutering, and UMAPs for stratified data. 
'''
# Define project, log, and output directories. 
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_noreplicates/stratify_with_zfpkd/
out_dir=/data/chanjlab/$project/output_noreplicates/stratify_with_zfpkd/
mkdir -p $out_dir $log_dir
cd $log_dir

# Define bsub parameters
bsub -J stratify -n 10 \
     -R span[ptile=10] -R rusage[mem=20] -W 5:00 \
     -o stratify.%J.stdout -eo stratify.%J.stderr \
     "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna_clone;
     /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/13_run_stratify_and_reprocess.py"