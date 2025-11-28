## General Setup

The original execution of the workflow were done on a virtual machine (VM) running the Ubuntu 20.04 operating system. The VM was created with 32-core CPU and 128 GB RAM for data processing.

Follow the instructions below for the basic setup procedure.

**Update**

Ensure the new virtual machine/instance running Linux/Ubuntu is up to date and upgraded. Confirm with defaults when prompted. This may take few minutes. 

```
sudo apt update && sudo apt upgrade -y  
```

**Install Anaconda3 to enable set up of conda environments and access conda archives**

Ensure curl is installed for transferring data using URLs. Download the bash installer for Anaconda3 and verify the checksum of the selected version.

```
sudo apt install curl -y 
curl -O https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
sha256sum Anaconda3-2023.09-0-Linux-x86_64.sh
```

The expected checksum for above example is: 
6c8a4abb36fbb711dc055b7049a23bbfd61d356de9468b41c5140f8a11abd851

Other version available at https://repo.anaconda.com/archive/.

**Execute the Anaconda installer**

```
bash Anaconda3-2023.09-0-Linux-x86_64.sh
```

Follow the license agreement, type "yes" when prompted. Press Enter to confirm the installation location "[/home/$USER/anaconda3]". Type 'no' for initialization prompt.

Make sure the PATH to Anaconda3 is added.

```
echo $PATH #check current PATH
export PATH=$PATH:/home/$USER/anaconda3/bin
echo $PATH #check new PATH
```

**Initialize Anaconda3**

```
conda init
```

Logout "ctr+D" and re-enter the instance. "(base)" should appear before the shell promt.
Check that Anaconda is working.

```
conda --version
```

You can now delete the installer.

```
sudo rm Anaconda3-2023.09-0-Linux-x86_64.sh
```

**Ensure conda is updated type "y" if promted (Optional)**

```
conda update -n base -c defaults conda
```

**Add required conda channels**

```
conda config --add channels defaults
conda config --add channels conda-forge
```
## Setup of Directory Architecture and Conda Environments (DovetailHiChIP, MACS2, PICARD)

With conda installed, we can now set up the directory structure for the analysis and create conda environments for the required software. Different tools may require specific versions of dependencies to avoid conflicts.

Clone the HiChIP-Analysis repository into a selected directory:

```
cd /home/<$USER>/
git clone https://github.com/ValdemarP267/HiChIP-Analysis.git
```
Ensure the whole folder has proper permissions. Run the DirectoryArchitecture.sh script and specify the main working directory when prompted.

```
sudo chmod 777 -R ./HiChIP-Analysis
bash ./HiChIP-Analysis/Scripts/DirectoryArchitecture.sh
```

To avoid memory issues, some pipeline steps write temporary files into a temp directory. The script DirectoryArchitecture.sh creates this folder along with all other required directories. Temporary files may occupy up to 3× the size of the fastq.gz files, so ensure enough storage space is available.

Next execute the HiChIPTools_install.sh

```
bash ./HiChIP-Analysis/HiChIPTools_install.sh
```

Confirm any prompts with defaults.

>[!NOTE]
>Running the script HiChIPTools_install.sh should create three conda environments:
> - DovetailHiChIP with trim-galore, fastqc and multiqc installed. Additional tools are installed following [HiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/HiChIP.md)
> - MACS2 with macs2, idr installed
> - Picard with picard installed

## Install Docker Engine (required for FitHiChIP)

Docker is a platform designed to help developers build, share, and run container applications.  Follow instructions below or visit official site for [Docker](https://docs.docker.com/engine/install/ubuntu/). 

**Uninstall old versions, conflicting packages**

Remove the following unofficial packages if installed:

- docker.io
- docker-compose
- docker-doc
- podman-docker

Run:

```
for pkg in docker.io docker-doc docker-compose podman-docker containerd runc; do sudo apt-get remove $pkg; done
```

**Instalaltion using the apt repository. Set up Docker's apt repository**

Install dependencies and add Docker’s GPG key:

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

This command downloads a test image and runs it in a container. It will prints a confirmation message and exits.

```
sudo docker run hello-world
```

**Linux post-installation steps for Docker Engine enable root user group permissions for Docker**

```
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

With above steps done you are ready to proceed with:

- [HiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/HiChIP.md) - Setup instructions for the Dovetail pipeline for alignments and pre-processing of HiChIP output.
- [FitHiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/FitHiChIP.md) - Setup instructions for FitHiChIP for loop calling.
- [Coolbox.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/Coolbox.md) - Setup of Coolbox visualization toolkit for genomic data imported to Jupyter Notebook for creating visuals and browse. 