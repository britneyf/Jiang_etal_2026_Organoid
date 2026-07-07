project="CRC_ZFP36L2.092023/Organoid"

mdir=`pwd`

# Log directory
log_dir="/data/chanjlab/$project/logs_new/impute_stratify_for_milopy/"
mkdir -p $log_dir

# Output directory
odir="/data/chanjlab/$project/output_new/milopy"
mkdir -p $odir

cd $log_dir

for adata in $odir/stratified*h5ad; do 
    patient=`basename $adata | sed 's/.h5ad//g' | sed 's/stratified_//g'`
    
   bsub -q cpuqueue -sla llSC4 -R A100 -J impute.$patient -n 20 -R span[ptile=20] -R rusage[mem=10] -W 2:00 -o impute.$patient.%J.stdout -eo impute.$patient.%J.stderr "conda activate multiome_env; /home/chanj3/anaconda3/envs/multiome_env/bin/python $mdir/GetImpute.py $adata"

done

cd -
