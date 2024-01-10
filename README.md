# HiChIP-Analysis

Documentation to analyze HiChIP data generated with Dovetail MNase-HiChIP kit.

This is a setup and execution guide to create working environments in conda and install the tools needed for analysis of TCF3::HLF HiChIP data generated from the HAL-01 cell line.
May be used for analysis of othere dataset with carefull adaptation of provided scripts.

This guide will take you through:

- Generation of a Linux operating Virtual Machine with the necessary base to proceed with installation of all needed tools, as well as directory architecture.

- Installation of Dovetail HiChIP pipeline and setup of Docker to pull docker image for FitHiChIP.

- Installation of MACS2 and IDR for peak calling and validation.

- Installation of tools for downstream analysis and vizualization of data such as deepTools, HOMER and coolbox.

## Original Documentations

Dovetail HiChIP
https://hichip.readthedocs.io/en/latest/

FitHiChIP
https://ay-lab.github.io/FitHiChIP/html/index.html

bwa
https://github.com/lh3/bwa

HiC-Pro
https://github.com/nservant/HiC-Pro

coolbox
https://gangcaolab.github.io/CoolBox/index.html

deepTools
https://deeptools.readthedocs.io/en/develop/index.html

Homer
http://homer.ucsd.edu/homer/index.html

>[!NOTE]
>If issues during installation occure or during execution of any of the tools refere to these documents for eventual troubleshooting.

## General Setup of Virtual Machine

Original work was done by running Virtual Machines on Science Cloud UZH (https://cloud.s3it.uzh.ch/auth/login/?next=/). Examples were done in an VM of a Ubuntu 20.04 system. The VM was created with 32 core CPU and 128GB RAM for data processing.

Instructions for launching an instance and connecting is not covered here. Go to https://docs.s3it.uzh.ch/cloud/training/training_handout/ for details. For basic setup procedure proceed with below instructions.

**Virtual machine update**

Ensure the new virtual machine/instance is up to date and upgraded. Confirm with defaults when prompted. This may take few mininutes. 

```
sudo apt update
sudo apt upgrade -y  
```

**Install Anaconda3 to set up conda environments and access conda archives**

Esure curl is installed for transferring data from or to a server using URLs. Download the bash file for installation of Anaconda3. Verify the checksum of selected install version

```
sudo apt install curl -y 
curl -O https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
sha256sum Anaconda3-2023.09-0-Linux-x86_64.sh
```

The expected checksum for above example is: 
6c8a4abb36fbb711dc055b7049a23bbfd61d356de9468b41c5140f8a11abd851

Other version available at https://repo.anaconda.com/archive/.

Execute the bash file

```
bash Anaconda3-2023.09-0-Linux-x86_64.sh
```

Follow the license agreement, type "yes" when prompted. Press Enter to confirm the installation location "[/home/ubuntu/anaconda3]". Type 'no' for initialization prompt.

Make sure the PATH to anaconda3 is added

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
Check the installed version of anaconda and confirm it's working.

```
conda --version
```

Now you can delet the Anaconda3-2023.09-0-Linux-x86_64.sh file.

```
sudo rm Anaconda3-2023.09-0-Linux-x86_64.sh
```

**Ensure conda is updated type "y" if promted (Optional)**

```
conda update -n base -c defaults conda
```

**Ensure all needed channels are added to conda**

```
conda config --add channels defaults
conda config --add channels conda-forge
```

## Setup of directory architecture and conda environments

With conda installed we are ready to setup the directory architecture we will use for our analsysis. In addition we create conda environments were our software will be installed. This is needed as some dependecies required may need to be of a certain versions for the tools to work without conflicts.

Clone source code of this repository into selected directory. 

```
cd /home/ubuntu/
git clone https://github.com/ValdemarP267/HiChIP-Analysis.git
```
Make sure whole folder has permission and run script DirectoryArchitecture.sh
The script assumes work is to be saved into a mounted volume called "/mnt".

```
sudo chmod 777 -R ./HiChIP-Analysis
./HiChIP-Analysis/DirectoryArchitecture.sh
```

To avoid memory issues, some of the pipeline steps require writing temporary files into a temp folder. Running the DirectoryArchitecture.sh will create this folder. Temporary files may take up to x3 of the space that the fastq.gz files are taking, make sure the working volume is big enough.

## Install Docker Engine needed for FitHiChIP tool

Docker is a platform designed to help developers build, share, and run container applications.  

**Uninstall old versions, conflicting packages**

The unofficial packages to uninstall are:

- docker.io
- docker-compose
- docker-doc
- podman-docker

Run the following command to uninstall all conflicting packages

```
for pkg in docker.io docker-doc docker-compose podman-docker containerd runc; do sudo apt-get remove $pkg; done
```

**Instalaltion using the apt repository. Set up Docker's apt repository**

Add Docker's official GPG key

```
sudo apt-get update
sudo apt-get install ca-certificates curl gnupg
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg
```

**Add the repository to apt sources**

```
echo \
  "deb [arch="$(dpkg --print-architecture)" signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  "$(. /etc/os-release && echo "$VERSION_CODENAME")" stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
```

To install the latest version, run:

```
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
```

**Verify that the Docker Engine installation is successful by running the hello-world image**

This command downloads a test image and runs it in a container. When the container runs, it prints a confirmation message and exits.

```
sudo docker run hello-world
```

**Linux post-installation steps for Docker Engine enable root user group permissions for docker**

```
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```