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
    conda install -y -c conda-forge -c bioconda $2
    conda deactivate
}

# Create and set up Conda environments
create_conda_env "DovetailHiChIP" "trim-galore=0.6.6 fastqc multiqc deeptools=3.5.1 matplotlib=3.8.0"
create_conda_env "MACS2" "python=3.7.11 macs2=2.2.6 idr=2.0.4.2 homer bedtools=2.30.0"
create_conda_env "Picard" "picard=2.25.7"

echo "Setup of Conda environments complete."