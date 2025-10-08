#!/usr/bin/env bash
set -e
shopt -s nullglob  # Make globbing return empty array if no match

###########################################################################
## HiChIPTools_install.sh
## Sets up necessary Conda environments and tools for TCF3::HLF HiChIP
###########################################################################

# Initialize Conda reliably
source ~/anaconda3/etc/profile.d/conda.sh

# Install Java (>=17)
echo "Installing Java..."
sudo apt update
sudo apt install -y openjdk-17-jre-headless

# Function to check if a Conda environment exists and create it if not
create_conda_env() {
    if ! conda env list | grep -qE "^\s*$1\s"; then
        echo "Creating Conda environment: $1"
        conda create -y -n "$1"
    else
        echo "Conda environment '$1' already exists."
    fi

    echo "Installing packages in $1: $2"
    conda activate "$1"
    conda install -y -c bioconda $2
    conda deactivate
}

# Create and set up Conda environments
create_conda_env "DovetailHiChIP" "trim-galore fastqc multiqc"
create_conda_env "MACS2" "macs2 idr homer bedtools"
create_conda_env "Picard" "picard"

echo "Setup of Conda environments complete."