#!/usr/bin/bash
set -e

#Script will create all directories needed for HiChIP and ChIP analysis.

# Initialize Conda
eval "$(conda shell.bash hook)"

# Function to check if a directory exists and create if it does not. Grants directories all permissions and ubuntu group ownership
create_directory() {
    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        sudo mkdir -m 777 -p "$1"
        sudo chown ubuntu:ubuntu "$1" #set for correct user
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

# List of subdirectories relative to MAIN_DIR
SUBDIRS=(
    "0.BlackList"
    "0.GenomeAssembly"
    "1.RawData"
    "1.RawData/HiChIP"
    "1.RawData/ChIP"
    "2.FASTQC"
    "2.FASTQC/HiChIP"
    "2.FASTQC/ChIP"
    "3.TRIM"
    "3.TRIM/HiChIP"
    "3.TRIM/ChIP"
    "4.ChIP_Alignment"
    "4.HiChIP_Alignment"
    "4.HiChIP_Alignment/Outputs"
    "4.ChIP_Alignment/Outputs"
    "5.MACS2"
    "5.MACS2/SORT"
    "5.MACS2/Permissive"
    "5.MACS2/IDR"
    "5.MACS2/D1"
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
    "10.ROSE"
    "tmp"
)

# Loop through and create each subdirectory
for sub in "${SUBDIRS[@]}"; do
    create_directory "$MAIN_DIR/$sub"
done

echo "All required directories created under: $MAIN_DIR"
