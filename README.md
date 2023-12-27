# HiChIP-Analysis
Dovetail, FITHiChIP and MACS2 documentations to analyze HiChIP data generated with Dovetail MNase-HiChIP kit.

This is a setup and command execution guide to create the necessary working environments in conda and install the tools needed for HiChIP analysis of TCF3::HLF HiChIP data generated from the HAL-01 cell line.

The repository will guide through:
- Generation of a Linux operating Virtual Machine with the necessary base to proceed with installation of all needed tools, as well as directory architecture.

- Installation of Dovetail HiChIP pipeline and FitHiChIP with necessary adaptation for TCF3::HLF HiChIP analysis

- Installation of tools for downstream analysis and vizualization of data

## Original Documentations
Dovetail HiChIP
https://hichip.readthedocs.io/en/latest/

FitHiChIP
https://ay-lab.github.io/FitHiChIP/html/index.html

bwa
https://github.com/lh3/bwa

## General Setup of Virtual Machine
Initial setup of a new instance, running on Science Cloud UZH (https://cloud.s3it.uzh.ch/auth/login/?next=/).Example was done in an Ubuntu 18.04 image from scratch. VM with 16 core CPU and 128GB RAM was created for data processing.

Access via the Ubuntu terminal on your working platform. Download Ubuntu application for your laptop to connect to the instance. Instructions for launching an instance and connecting is not covered here. Go to https://docs.s3it.uzh.ch/cloud/training/training_handout/ for details. Proceed with below instructions.

**Virtual machine update**

Ensure the new virtual machine/instance is up to date and upgraded. Confirm with defaults when prompted. This may take few mininutes. 

```
sudo apt update
sudo apt upgrade -y  
```

**Install Anaconda3 to set up Conda environments and access Conda archives**

Install tool for transferring data from or to a server using URLs. Download the bash file for installation of Anaconda3. Verify the checksum of selected install version
```
sudo apt install curl -y 
curl -O https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
sha256sum Anaconda3-2022.10-Linux-x86_64.sh
```
The expected checksum for above example is: 
e7ecbccbc197ebd7e1f211c59df2e37bc6959d081f2235d387e08c9026666acd

Other version available at https://repo.anaconda.com/archive/ .

Execute the the bash file
```
bash Anaconda3-2022.10-Linux-x86_64.sh
```

Follow the license agreement, type "yes" when prompted. Press Enter to confirm the installation location "[/home/ubuntu/anaconda3]". Type 'no' for initialization prompt.

**Make sure the PATH to Anaconda is added**
```
echo $PATH #check current PATH
export PATH=$PATH:/home/ubuntu/anaconda3/bin
echo $PATH #check new PATH
```

**Initialize Anaconda3**
```
conda init
```
Logout "ctr+D" and re-enter the instance. "(base)" should show before the shell name.
Check the installed version of Anaconda and confirm it's working.
```
conda --version
```
Now you can delet the Anaconda3-2022.10-Linux-x86_64.sh file.
```
sudo rm Anaconda3-2022.10-Linux-x86_64.sh
```

**Ensure conda is updated type "y" if promted**
```
conda update -n base -c defaults conda
```

**Ensure all channels are added to conda**
```
conda config --add channels defaults
conda config --add channels conda-forge
```
**Run DirectoryArchitecture.sh**
```
git clone https://github.com/ValdemarP267/HiChIP-Analysis/tree/main
```
Make DirectoryArchitecture.sh script executable and run it:
```
chmod +x ./HiChIP-Analysis/DirectoryArchitecture.sh
./HiChIP-Analysis/DirectoryArchitecture.sh
```
This will create all the conda environments and directory architecture for executing the HiChIP Analsysis.
