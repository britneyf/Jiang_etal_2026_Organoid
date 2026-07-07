# By Britney T Forsyth
# Description: Python script to run Cell Bender.

#!/bin/bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate cellbender

#date=020624
batch_name=$1
path_name=$2
expected_cells=$3

# Read the batch name, files, and expected cells
batch_name=$(echo "$batch_name" | sed 's#//*#/#g')
path_name=$(echo "$path_name" | sed 's#//*#/#g')

# Debugging: Print normalized paths and filenames
echo "Batch Name: $batch_name"
echo "Path Name: $path_name"

# The output directory is /data/chanjlab/CRC_ZFP36L2.092023/Organoid/out_cellbender/$batch_name/cellbender
# $project is CRC_ZFP36L2.092023/Organoid
odir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/cellbender/$batch_name/

# In the output directory, we make a folder for input and output directory
mkdir -p $odir/input $odir/output

# Copy all files from path_name to the input directory
cp "$path_name"/* "$odir/input"

# Output files will be in /data/chanjlab/CRC_ZFP36L2.092023/Organoid/out_cellbender/$batch_name/cellbender/output/cellbender.$batch_name.combined.h5
ofile=$odir/output/cellbender.$batch_name.combined.h5

cd $odir/output
module load gcc/10.2.0 
module load cuda/11.7 
conda activate cellbender

# If the input directory is empty, copy the h5 file from the specified location
#if [ "$(ls -A $odir/input)" ]; then
    :
#else
    #cp $h5_file $odir/input
#fi

/home/forsythb/anaconda3/envs/cellbender/bin/cellbender remove-background \
--input $odir/input \
--output $ofile \
--cuda \
--expected-cells $expected_cells