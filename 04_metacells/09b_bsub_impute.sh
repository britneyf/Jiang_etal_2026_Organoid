'''
DESCRIPTION:
Bsub script to impute data.
'''
# Define the project, log, and output directory
project="CRC_ZFP36L2.092023/Organoid"

# Log directory
log_dir="/data/chanjlab/$project/logs_noreplicates/impute_stratify/"
mkdir -p $log_dir

# Output directory
odir="/data/chanjlab/$project/output_noreplicates/impute"
mkdir -p $odir
cd $log_dir

# Input directory
idir="/data/chanjlab/$project/output_noreplicates/stratify"

# Read in adata and set up bsub parameters
for adata in $idir/*h5ad; do 
    patient=`basename $adata | sed 's/.h5ad//g'`
    
    bsub -q cpuqueue -sla llSC4 -R A100 -J impute.$patient -n 20 -R span[ptile=20] -R rusage[mem=10] -W 2:00 -o impute.$patient.%J.stdout -eo impute.$patient.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/09_run_impute.py $adata $patient"
done
cd -