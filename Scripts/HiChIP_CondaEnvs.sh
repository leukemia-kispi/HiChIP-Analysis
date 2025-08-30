#!/usr/bin/bash
set -e

#To create the different conda enviorments used during HiChIP analsyis with some of the necessary tools installed

## Make sure Java >= 11 is installed**
echo "Installing Java."

sudo apt install openjdk-17-jre-headless -y


# Initialize Conda
eval "$(conda shell.bash hook)"

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
create_conda_env "DovetailHiChIP" "trim-galore fastqc multiqc"  #Additional tools will be installed by following the Dovetail HiChIP setup described in HiChIP.md. Including bwa-mem, samtools, pairtools, bedtools, deeptools,
create_conda_env "MACS2" "macs2 idr homer bedtools"
create_conda_env "Picard" "picard"

echo "Setup of conda environments complete."

