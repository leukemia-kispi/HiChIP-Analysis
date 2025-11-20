#!/usr/bin/env bash

################################################################################
## These are the step to install all the dependencies
## used in the examples of this repo
## on Ubuntu 20.04 modefied from Original script to enforce versions
################################################################################
sudo apt-get -y update && sudo -y apt-get upgrade


## Install GCC make and python and pip
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt-get update
sudo apt-get install -y build-essential python3.10 python3.10-distutils python3-pip
# Ensure Python 3.10.8
sudo apt-get install -y python3.10=3.10.8-1+focal1
# Ensure pip 23.0.1
python3.10 -m pip install --upgrade pip==23.0.1


## This will put python3 as the default binairies call by python
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1
sudo update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

################################################################################
## Install bedtools
################################################################################
sudo apt-get -y install bedtools=2.30.0

################################################################################
## Install deeptools and pairtools
################################################################################
## these are required for pyBigWig
sudo apt-get install -y zlib1g-dev libcurl4

pip3 install \
    numpy==1.5.0 \
    pysam==0.21.0 \
    tabulate==0.9.0 \
    scipy==1.10.1 \
    py2bit==0.3.0 \
    matplotlib==3.7.1 \
    pyBigWig==0.3.18 \
    deeptools==3.5.1 \
    pandas==2.0.0


pip3 install pairtools==1.0.2

################################################################################
## Install bwa
################################################################################
git clone https://github.com/lh3/bwa.git
cd bwa; make -j $(nproc)
ls
./bwa
sudo cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin
cd

################################################################################
## Install samtools
################################################################################
# First install the following dependencies:
sudo apt-get install -y libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev

## Then download, compile and assemble samtools
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar xvf samtools-1.11.tar.bz2 

cd samtools-1.11/
./configure 
make -j $(nproc)
sudo make install
cd

###############
##install liblz4
##############
sudo apt-get install liblz4-tool
