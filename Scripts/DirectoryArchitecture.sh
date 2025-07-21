#!/bin/bash
set -e

#Script will create all directories needed for HiChIP and ChIP analysis and different conda enviorments with the necessary tools installed

# Initialize Conda
eval "$(conda shell.bash hook)"

# Function to check if a directory exists and create if it does not. Grants directories all permissions and ubuntu group ownership
create_directory() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        sudo mkdir -m 777 -p "$1"
        sudo chown ubuntu:ubuntu "$1"
    else
        echo "Directory already exists: $1"
    fi
}

# Prompt the user for the main directory
read -rp "Enter the path to the MAIN_DIR: " MAIN_DIR

# Check if input is empty
if [ -z "$MAIN_DIR" ]; then
    echo "Error: No directory path entered."
    exit 1
fi

# Ensure the main directory itself exists
create_directory "$MAIN_DIR"

# Define subdirectories relative to the main folder
SUBDIRS=("logs" "data" "output/reports" "tmp/cache")


# Loop through and create each subdirectory
for dir in "${SUBDIRS[@]}"; do
    full_path="$MAIN_DIR/$dir"
    create_directory "$full_path"
done

# List of subdirectories relative to MAIN_DIR
SUBDIRS=(
    "0.BlackList"
    "0.GenomeAssembly"
    "1.RawData"
    "1.RawData/HiChIP"
    "1.RawData/ChIP"
    "2.FASTAQC"
    "2.FASTAQC/HiChIP"
    "2.FASTAQC/ChIP"
    "3.TRIM/HiChIP"
    "3.TRIM/ChIP"
    "4.ChIP_Alignment"
    "4.HiChIP_Alignment"
    "4.HiChIP_Alignment/Outputs"
    "5.MACS2"
    "5.MACS2/SORT"
    "5.MACS2/Permissive"
    "5.MACS2/IDR"
    "6.FitHiChIP_Output"
    "7.Deeptool_Matrix"
    "7.Deeptool_Matrix/Coverage"
    "7.Deeptool_Matrix/Outputs"
    "7.Deeptool_Matrix/Graphs"
    "7.Deeptool_Matrix/Sorted_Lists"
    "7.Deeptool_Matrix/Target_Regions"
    "8.JupyterLab"
    "9.HOMER"
    "9.HOMER/Outputs"
    "tmp"
)

# Loop through and create each subdirectory
for sub in "${SUBDIRS[@]}"; do
    create_directory "$MAIN_DIR/$sub"
done

echo "All required directories created under: $MAIN_DIR"

# Function to check if a Conda environment exists and create it if not with respective applications installed
create_conda_env() {
    if ! conda env list | grep -q "$1"; then
        conda create -y -n "$1"
    fi
    conda activate "$1"
    conda install -y -c bioconda "$2"
    conda deactivate
}
    
# Create and set up Conda environments
create_conda_env "DovetailHiChIP" "trim-galore fastqc multiqc"
create_conda_env "MACS2" "macs2 idr homer bedtools"

echo "Setup of conda enviroments complete."

