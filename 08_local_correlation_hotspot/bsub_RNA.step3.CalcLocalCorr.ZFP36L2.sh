project=CRC_ZFP36L2.092023/Organoid
data=080524
log_dir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/logs_noreplicates/calclocalcorr/$date"
mkdir -p $log_dir
# Output directory
odir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_noreplicates/calclocalcorr/$date"
mkdir -p $odir
cd $log_dir

cd $log_dir

bsub -q cpuqueue -J LocalCorr -n 30 -R rusage[mem=5] -W 10:00 -o LocalCorr.%J.stdout -eo LocalCorr.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/CalcLocalCorr.ZFP36L2.py"

cd -
