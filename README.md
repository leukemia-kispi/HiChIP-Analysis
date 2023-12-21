# HiChIP-Analysis
Dovetail, FITHiChIP and MACS2 documentations to analyze HiChIP data generated with Dovetail MNase-HiChIP kit.

This is a setup and command execution guide to create the necessary working environments in conda to install the tools needed for HiChIP analysis

## Original Documentations
Dovetail HiChIP
https://hichip.readthedocs.io/en/latest/

FitHiChIP
https://ay-lab.github.io/FitHiChIP/html/index.html

## General Setup of Virtual Machine
Initial setup of a new instance, running on Science Cloud UZH (https://cloud.s3it.uzh.ch/auth/login/?next=/).Example was done in an Ubuntu 22.04 image from scratch. Consider amount of cpu and RAM needed to process your data when selecting flavor of the image.

Access via the Ubuntu terminal on your working platform. Download Ubuntu application for your laptop to connect to the instance. Instructions for launching an instance and connecting is not covered here. Go to https://docs.s3it.uzh.ch/cloud/training/training_handout/ for details. Proceed with below instructions.

**Virtual machine update**

Ensure the new virtual machine/instance is up to date and upgraded. Confirm with defaults when prompted. This may take few mininutes. sudo preceeding the commands ensures their executions with root permissions. 

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
**Make sure Python >= 3.6 is installed. Install Java >= 11 and pip if needed**
```
sudo apt install python3 python3-pip -y
sudo apt install openjdk-17-jre-headless -y
```
Confirm any prompts with defaults.

**Ensure conda is updated type "y" if promted**
```
conda update -n base -c defaults conda
```

**Ensure all channels are added to conda**
```
conda config --add channels defaults
conda config --add channels conda-forge
```
VM ready for use