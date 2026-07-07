'''
AUTHOR: 
Britney T. Forsyth

DESCRIPTION:
Bsub script to run SCRAN. Performs normalization of adata using counts matrix, gene, and barcodes file from concatenation.
Uses a deconvolution approach to partition cells into pools and performs normalization across cells in each pool. 
'''

# Define project directory, log directory, and output directory
project=CRC_ZFP36L2.092023/Organoid
log_dir=/data/chanjlab/$project/logs_noreplicates/scran/
odir=/data/chanjlab/$project/output_noreplicates/scran/
mkdir -p $odir $log_dir

# Use the concatenated matrix, gene, and barcodes files from the concatenation step as input
mtx_file=/data/chanjlab/$project/output_noreplicates/concatenated/counts.RNA.combined.mtx
g_file=/data/chanjlab/$project/output_noreplicates/concatenated/counts.RNA.combined.genes.csv
bc_file=/data/chanjlab/$project/output_noreplicates/concatenated/counts.RNA.combined.barcodes.csv

# Also define the path to the output file
#ofile=/data/chanjlab/$project/output_noreplicates/scran/scran.matrix.mtx

# Set up the bsub parameters
bsub -J scran -n 10 -R span[ptile=10] -R rusage[mem=20] -W 5:00 -o scran.%J.stdout -eo scran.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate r_4.3.1; ~/anaconda3/envs/r_4.3.1/bin/Rscript --vanilla /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/05_run_scran.R $mtx_file $g_file $bc_file"