#!/bin/bash

#Script will create all directories needed for HiChIP analysis and different conda enviorments with the necessary tools installed

# Initialize Conda
eval "$(conda shell.bash hook)"

# Function to check if a directory exists and create if it does not. Grants directories all permissions and ubuntu group ownership
create_directory() {
    if [ ! -d "$1" ]; then
        sudo mkdir -m 777 "$1"
        sudo chown ubuntu:ubuntu "$1"
    fi
}

# Check and create directories in mounted Volume designated /mnt
create_directory "/mnt/0.BlackList"
create_directory "/mnt/0.GenomeAssembly" 
create_directory "/mnt/1.RawData"
create_directory "/mnt/2.FASTAQC"
create_directory "/mnt/3.TRIM"
create_directory "/mnt/4.HiChIP_Alignment"
create_directory "/mnt/4.HiChIP_Alignment/Outputs"
create_directory "/mnt/5.MACS2"
create_directory "/mnt/5.MACS2/SORT"
create_directory "/mnt/5.MACS2/Permissive"
create_directory "/mnt/5.MACS2/IDR"
create_directory "/mnt/6.FitHiChIP_Output"


# Function to check if a Conda environment exists and create it if not with respective applications installed
create_conda_env() {
    if ! conda env list | grep -q "$1"; then
        conda create -n "$1"
    fi
    conda activate "$1"
    conda install -y -c bioconda "$2"
    conda deactivate
}
    
# Create and set up Conda environments
create_conda_env "DovetailHiChIP"
create_conda_env "FitHiChIP"
create_conda_env "TRIM" "trim-galore"
create_conda_env "FASTAQC" "fastqc multiqc"
create_conda_env "MACS2" "macs2 idr bedtools"
create_conda_env "COOLBOX" "coolbox"
create_conda_env "HOMER" "homer"


echo "Setup complete."
