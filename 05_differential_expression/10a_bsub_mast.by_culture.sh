# By: Britney Forsyth

# This is the project directory
project=CRC_ZFP36L2.092023/Organoid

# Make log directory and output directory
log_dir=/data/chanjlab/$project/logs_noreplicates/mast/mast.by_culture.nonorm/
out_dir=/data/chanjlab/$project/output_noreplicates/mast/mast.by_culture.nonorm/
mkdir -p $out_dir $log_dir

# Directory where .h5ad files are located
data_dir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/stratify_by_culture/"

# Loop over each .h5ad file in the specified directory
for file_path in "${data_dir}"/*.h5ad; do
    # Extract base filename without extension
    culture=$(basename "${file_path}" .h5ad)

    # Submit job using bsub
    bsub -sla llSC4 -R A100 -q cpuqueue -J "mast_${culture}" -n 15 \
         -R "span[ptile=10]" -R "rusage[mem=20]" -W 24:00 \
         -o "${log_dir}/mast.${culture}.%J.stdout" -eo "${log_dir}/mast.${culture}.%J.stderr" \
         "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/CRC_ZFP36L2.092023/Organoid/scripts/organoid_analysis_pipeline_scripts/10a_run_mast.by_culture.R ${culture}"
done