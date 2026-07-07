#!/bin/bash
# By Britney T Forsyth
# Description: Generates a batch description text file with information about sample batches, their paths, and expected cell counts

# Set the base directory
base_dir="/data/chanjlab/CRC_ZFP36L2.092023/Organoid/raw.multiplex_copy"

# Set the output directory
output_dir="/data/chanjlab/forsythb/organoid_analysis_pipeline_scripts"

# Name the output file
output_file="$output_dir/batch_description.txt"

# Create the header in the output file
echo -e "batch_name\tpath_name\texpected_cells" > "$output_file"

# Iterate through each batch folder
for batch_folder in "$base_dir"/*; do
    if [ -d "$batch_folder" ]; then
        # Extract batch name from the folder path
        batch_name=$(basename "$batch_folder")

        # Construct the path name
        path_name="$batch_folder/raw_feature_bc_matrix"

        # Initialize the expected cells count
        expected_cells=0

        # Iterate through each sample folder in the demultiplex folder
        demultiplex_folder="$batch_folder/demultiplex"
        if [ -d "$demultiplex_folder" ]; then
            for sample_folder in "$demultiplex_folder"/*; do
                if [ -d "$sample_folder" ]; then
                    # Construct the path to the barcodes.tsv.gz file
                    barcodes_file="$sample_folder/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"

                    # Count the lines in the gzip file
                    lines=$(zcat "$barcodes_file" | wc -l)

                    # Add to the expected cells count
                    expected_cells=$((expected_cells + lines))
                fi
            done
        fi

        # Write the information to the output file
        echo -e "$batch_name\t$path_name\t$expected_cells" >> "$output_file"
    fi
done