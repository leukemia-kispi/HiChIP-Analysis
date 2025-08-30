## General Setup

Original analysis were done on a virtual machine running Ubuntu 20.04 operating system. The VM was created with 32 core CPU and 128GB RAM for data processing.

For basic setup procedure proceed with below instructions.

**Update**

Ensure the new virtual machine/instance is up to date and upgraded. Confirm with defaults when prompted. This may take few mininutes. 

```
sudo apt update
sudo apt upgrade -y  
```

**Install Anaconda3 to set up conda environments and access conda archives**

Ensure curl is installed for transferring data from or to a server using URLs. Download the bash file for installation of Anaconda3. Verify the checksum of selected install version

```
sudo apt install curl -y 
curl -O https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
sha256sum Anaconda3-2023.09-0-Linux-x86_64.sh
```

The expected checksum for above example is: 
6c8a4abb36fbb711dc055b7049a23bbfd61d356de9468b41c5140f8a11abd851

Other version available at https://repo.anaconda.com/archive/.

**Execute the bash file**

```
bash Anaconda3-2023.09-0-Linux-x86_64.sh
```

Follow the license agreement, type "yes" when prompted. Press Enter to confirm the installation location "[/home/<USERNAME>/anaconda3]". Type 'no' for initialization prompt.

Make sure the PATH to anaconda3 is added

```
echo $PATH #check current PATH
export PATH=$PATH:/home/<USERNAME>/anaconda3/bin
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
## Setup of directory architecture and DovetailHiChIP, MACS2, PICARD conda environments

With Conda installed, we are now ready to set up the directory structure for our analysis. Additionally, we will create Conda environments where the required software will be installed. Different tools may require specific versions of their dependencies in order to function properly without conflicts.

Clone the source code of the HiChIP-Analysis repository into a selected directory. 

```
cd /home/<USERNAME>/
git clone https://github.com/ValdemarP267/HiChIP-Analysis.git
```
Make sure whole folder has permission. Run script DirectoryArchitecture.sh. Specefiy the main working directory when promted.

```
sudo chmod 777 -R ./HiChIP-Analysis
bash ./HiChIP-Analysis/DirectoryArchitecture.sh
```

To avoid memory issues, some of the pipeline steps require writing temporary files into a temp folder. Running the DirectoryArchitecture.sh will create this folder in addition to other directories that will be used during the workflow. Temporary files may take up to x3 of the space that fastq.gz files do, make sure the working volume is sufficient.

Next execute the CondaEnv.sh

```
bash ./HiChIP-Analysis/CondaEnv.sh
```

Confirm any prompts with defaults.

>[!NOTE]
>Running the script CondaEnv.sh should create three conda environments:
> - DovetailHiChIP with trim-galore, fastqc and multiqc installed
> - MACS2 with macs2, idr installed
> - Picard with picard isntalled

## Install Docker Engine needed for FitHiChIP tool

Docker is a platform designed to help developers build, share, and run container applications.  Follow instructions below or visit offical site for [Docker](https://docs.docker.com/engine/install/ubuntu/).

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

This command downloads a test image and runs it in a container. It will prints a confirmation message and exits.

```
sudo docker run hello-world
```

**Linux post-installation steps for Docker Engine enable root user group permissions for docker**

```
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

With above steps done you are ready to proceed with:

>[HiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/HiChIP.md) - Setup instructions for the Dovetail pipleline for aligments and pre-processing of HiChIP output.
>[FitHiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/FitHiChIP.md) - Setup instructions for FitHiChIP for loop calling.
>[Coolbox.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/Coolbox.md) - Setup of Coolbox visualization toolkit for genomic data imported to Jupyter Notebook for creating visuals and browse. 