# Set the directories
log_dir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/logs_noreplicates/hotspot
mkdir -p $log_dir

cd $log_dir

bsub -sla llSC4 -R A100 -q cpuqueue -J hotspot_patient_nontumor -n 30 -R rusage[mem=5] -W 10:00 -o hotspot_patient_nontumor.%J.stdout -eo hotspot_patient_nontumor.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; /home/forsythb/anaconda3/envs/scrna/bin/python /data/chanjlab/CRC_ZFP36L2.092023/Organoid/scripts/organoid_analysis_pipeline_scripts/run_hotspot_patient.py"

cd -
