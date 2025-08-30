#!/usr/bin/bash
set -e

#Script will create all directories needed for HiChIP and ChIP analysis.

# Initialize Conda
eval "$(conda shell.bash hook)"

# Function to check if a directory exists and create if it does not. Grants directories all permissions and ubuntu group ownership

create_directory() {

    user="${SUDO_USER:-$USER}"
    group="$(id -gn "$user")"

    if [ ! -d "$1" ]; then
        echo "Creating directory: $1"
        sudo mkdir -m 777 -p "$1"
    else
        echo "Directory already exists: $1"
    fi

    #Enforce ownership and permission
    sudo chown "$user:$group" "$1" 
    sudo mkdir -m 777 -p "$1"    
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
    "4.Alignment"
    "4.Alignment/HiChIP"
    "4.Alignment/ChIP"
    "5.MACS2"
    "5.MACS2/HiChIP"
    "5.MACS2/HiChIP/D1_Peaks"
    "5.MACS2/HiChIP/SORT"
    "5.MACS2/HiChIP/Permissive"
    "5.MACS2/HiChIP/IDR"
    "5.MACS2/ChIP"
    "5.MACS2/ChIP/Peaks"
    "5.MACS2/ChIP/SORT"
    "5.MACS2/ChIP/Permissive"
    "5.MACS2/ChIP/IDR"
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
