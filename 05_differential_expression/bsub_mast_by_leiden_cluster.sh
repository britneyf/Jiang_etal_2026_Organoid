# By: Britney Forsyth

project=CRC_ZFP36L2.092023/Organoid

# Make log directory and output directory
log_dir=/data/chanjlab/$project/logs_new/mast/
out_dir=/data/chanjlab/$project/output_new/mast_25/
mkdir -p $out_dir $log_dir

cd $log_dir

# Loop through unique leiden values
unique_leiden_values=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11")

for leiden_value in "${unique_leiden_values[@]}"; do
    bsub -J "mast_leiden_$leiden_value" -n 15 \
         -R "span[ptile=10]" -R "rusage[mem=20]" -W 24:00 \
         -o "$log_dir/mast_leiden_$leiden_value.%J.stdout" -eo "$log_dir/mast_leiden_$leiden_value.%J.stderr" \
         "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/run_mast_by_leiden_cluster.R $leiden_value"
done


# # Loop through unique leiden values
# bsub -J mast_leiden -n 15 -R span[ptile=10] -R rusage[mem=30] -W 10:00 -o mast_leiden.%J.stdout -eo mast_leiden.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/run_mast_by_leiden_cluster.R"