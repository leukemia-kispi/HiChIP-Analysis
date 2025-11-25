# FitHiChIP Loop calling

## FitHiChIP setup

FitHiChIP will be run from Docker image. The instructions in the [GeneralSetup.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/GeneralSetup.md) will guide you through docker installation. This simplefies FitHiChIP installation as dependecies dont need to be installed from scratch.

Select the directory to contain the FitHiChIP source code, and clone it.

```
cd /home/<$USER>
git clone https://github.com/ay-lab/FitHiChIP.git
sudo chmod 777 -R FitHiChIP
```

Data inputs needed to run FitHiChIP:

- Pairs files generated during HiChIP aligment as described in [HiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/HiChIP.md). These have to be  filtered and converted to HiC-Pro format (done by the HiChIP_mergedRep.sh or HiChIP_singleRep.sh)
- MACS2 called peaks from relevant ChIP-seq data or D1 MACS2 called peaks from primary algimnents file generated from HiChIP aligment data. Described in [HiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/HiChIP.md)
- Config file(example provided in this repository) specefying file locations and parameters

## HiC-Pro Installation

HiC-Pro needs to be installed on the system before the FitHiChIP docker can be run.

Clone the HiC-Pro source code and give it permissions

```
cd /home/$USER
git clone https://github.com/nservant/HiC-Pro.git
sudo chmod 777 -R HiC-Pro
```

Use the environment.yml file to create a conda environment with HiC-Pro and its dependecies installed. Activate the environment.
Set the environment on the Path to be able to run it when using docker for FitHiChIP

```
conda env create -f HiC-Pro/environment.yml -p HiCPro_ENV
conda activate /home/$USER/HiCPro_ENV/
```

Make sure unzip is installed in /home/$USER/HiCPro_ENV/.

```
sudo apt install unzip
```

Complete installation by entering the HiC-Pro directory and add the "[/home/$USER/HiCPro_ENV/]" path to PREFIX in config-install.txt file using the nano editor tool. This will be the installation location for HiC-Pro and some of the dependecies. Make sure $USER is the username on the system.

```
cd /home/$USER/HiC-Pro
nano config-install.txt
```

Complete the installation and add HiC-Pro to the $PATH.

```
make configure
make install
export PATH="/home/$USER/HiCPro_ENV/HiC-Pro/bin/":$PATH
```

## Running FitHiChIP Loop calling

**Filter pairs**

The "paritools select" command in the DovtailHiChIP conda environment will filter and include unique and rescue mapped read pairs. Done by the HiChIP_mergedRep_Alignment.sh or HiChIP_singleRep_Alignment.sh scripts.

**HiCPro Valid Pairs Files**

The paired files generated during aligment needs to be converted into valid HiC-Pro files to proceed with FitHiChIP. Done by the HiChIP_mergedRep_Alignment.sh or HiChIP_singleRep_Alignment.sh scripts.

**MACS2 D1 Peak calling**

The published ChIP-seq data for TCF3::HLF in HAL-01 did not provide clear and abundant peaks to let us define primary peaks within our HiChIP data using FitHiChIP. We instead used the TCF3::HLF HiChIP data itself to generate a file with primary aligments.

Running the MACS2_D1_HiChIP.sh script will generate several macs2.narrowpeak files for both merged replicate samples and single replicate samples. These can also be used for validation of high confidence peaks using the IDR tools that should also be installed.

Lastly the important output used for FitHiChIP loop calling are the HiChIP_<CellLine>_<conditions>_merged_nodd_hicpro_mapped.filtered.pairs.gz files.

**FitHiChIP Loop Calling**
Run FitHiChIP via bash script using docker image.
Make sure the configfile is properly setup and that all input files are on the specefied paths. Read the documentation from [Dovetail HiChIP](https://hichip.readthedocs.io/en/latest/) and [FitHiChIP](https://ay-lab.github.io/FitHiChIP/html/index.html) for details.

```
conda activate /home/$USER/HiCPro_ENV/
cd /home/$USER/FitHiChIP/
bash ./FitHiChIP_Docker.sh -C /home/$USER/HiChIP-Analysis/<$CONFIGFILE>
