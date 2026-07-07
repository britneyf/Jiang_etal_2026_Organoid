# By: Britney Forsyth

# This is the project directory
project=CRC_ZFP36L2.092023/Organoid

# Make log directory and output directory
log_dir=/data/chanjlab/$project/logs_noreplicates/gseapy/
out_dir=/data/chanjlab/$project/output_noreplicates/gseapy/
mkdir -p $out_dir $log_dir

# Directory where .csv files are located
data_dir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/mast/mast_nonormalization"

# Loop over each .csv file in the specified directory
for file_path in "${data_dir}"/*.csv; do
    # Extract base filename without extension
    mast_file=$(basename "${file_path}")

    # Submit job using bsub
    bsub -sla llSC4 -R A100 -q cpuqueue -J "gsea.${mast_file}" -n 15 \
         -R "span[ptile=10]" -R "rusage[mem=20]" -W 24:00 \
         -o "${log_dir}/gsea.${mast_file}.%J.stdout" -eo "${log_dir}/gsea.${mast_file}.%J.stderr" \
         "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna_clone; /home/forsythb/anaconda3/envs/scvi/bin/python /data/chanjlab/CRC_ZFP36L2.092023/Organoid/scripts/organoid_analysis_pipeline_scripts/11a_run_gseapy.py ${mast_file}"
done