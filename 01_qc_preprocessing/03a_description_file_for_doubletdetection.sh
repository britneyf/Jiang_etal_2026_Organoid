#!/bin/bash
# By [Your Name]
# Description: Collects paths to sample matrix files and CellBender files for doublet detection

# Set the base directory
base_dir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/raw.multiplex_copy"

# Set the CellBender output directory
cellbender_output_dir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/cellbender"

# Set the output directory
output_dir="/data/chanjlab/forsythb/organoid_analysis_pipeline_scripts"

# Name the output file
output_file="$output_dir/samplematrix_cellbenderpaths_for_doubletdetection.txt"

# Create the header in the output file
echo -e "batch_name\tsample_name\tsample_matrix_path\tbatch_path_to_cellbender" > "$output_file"

# Iterate through each batch folder
for batch_folder in "$base_dir"/*; do
    if [ -d "$batch_folder" ]; then
        # Extract batch name from the folder path
        batch_name=$(basename "$batch_folder")

        # Iterate through each demultiplex folder within the batch folder
        demultiplex_folder="$batch_folder/demultiplex"
        if [ -d "$demultiplex_folder" ]; then
            for sample_folder in "$demultiplex_folder"/*; do
                if [ -d "$sample_folder" ]; then
                    # Extract sample name from the folder path
                    sample_name=$(basename "$sample_folder")

                    # Construct the path to the sample matrix file
                    sample_matrix_path="$sample_folder/count/sample_filtered_feature_bc_matrix/"

                    # Construct the path to the CellBender file
                    cellbender_output_path="$cellbender_output_dir/$batch_name/output/cellbender.$batch_name.combined.h5"

                    # Write the information to the output file
                    echo -e "$batch_name\t$sample_name\t$sample_matrix_path\t$cellbender_output_path" >> "$output_file"
                fi
            done
        fi
    fi
done
