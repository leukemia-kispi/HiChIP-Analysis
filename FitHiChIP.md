# FitHiChIP Loop calling

## FitHiChIP setup

FitHiChIP will be run from Docker image. The instructions in the [README.md](https://github.com/ValdemarP267/HiChIP-Analysis) will guide you through docker installation. This simplefies FitHiChIP installation as dependecies dont need to be installed from scratch.

Select the directory to contain the FitHiChIP source code, and clone it.

```
cd /home/ubuntu
git clone https://github.com/ay-lab/FitHiChIP.git
sudo chmod 777 -R FithHiChIP
```

Data inputs needed to run FitHiChIP:

- Pairs files generated during HiChIP aligment [HiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/HiChIP.md) converted to HiC-Pro format (See below)
- MACS2 called peaks from relevant ChIP-seq data or D1 MACS2 called peaks from primary algimnents created from HiChIP aligment data by following instructions in [HiChIP.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/HiChIP.md)
- Config file(example provided in this repository) specefying file locations and parameters

## HiC-Pro Installation

HiC-Pro needs to be installed on the system before the FitHiChIP docker can be run.

Clone the HiC-Pro source code and give it permissions

```
cd /home/ubuntu
git clone https://github.com/nservant/HiC-Pro.git
sudo chmod 777 -R HiC-Pro
```

Use the environment.yml file to create a conda environment with HiC-Pro and its dependecies installed. Activate the environment.
Set the environment on the Path to be able to run it when using docker for FitHiChIP

```
conda env create -f HiC-Pro/environment.yml -p HiCPro_ENV
conda activate /home/ubuntu/HiCPro_ENV/
```

Make sure unzip is installed in /home/ubuntu/HiCPro_ENV/.

```
sudo apt install unzip
```

Complete installation by entering the HiC-Pro directory and add the "[/home/ubuntu/HiCPro_ENV/]" path to PREFIX in config-install.txt file using the nano editor tool. This will be the installation location for HiC-Pro and some of the dependecies.

```
cd /home/ubuntu/HiC-Pro
nano config-install.txt
```

Complete the installation and add HiC-Pro to the $PATH.

```
make configure
make install
export PATH="/home/ubuntu/HiCPro_ENV/HiC-Pro/bin/":$PATH
```

## Running FitHiChIP Loop calling

**Filter pairs**

The "paritools select" command in the DovtailHiChIP conda environment will filter and include unique and rescue mapped read pairs.

```
conda activate DovetailHiChIP
cd /mnt/4.HiChIP_Alignment
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs -o JoinedRep_TCF3_HLF_hg38_nodd_mapped.filtered.pairs
```

Same can be done for single replicate generated files.

```
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' Rep1_TCF3_HLF_hg38_nodd_mapped.pairs -o Rep1_TCF3_HLF_hg38_nodd_mapped.filtered.pairs
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' Rep2_TCF3_HLF_hg38_nodd_mapped.pairs -o Rep2_TCF3_HLF_hg38_nodd_mapped.filtered.pairs
```

**HiCPro Valid Pairs Files**

The paired files generated during aligment needs to be converted into valid HiC-Pro files to proceed with FitHiChIP.

```
conda activate DovetailHiChIP
cd /mnt/4.HiChIP_Alignment
grep -v '#' JoinedRep_TCF3_HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > JoinedRep_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz
```

Same for single replicate generated files.

```
grep -v '#' Rep1_TCF3_HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > Rep1_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz
grep -v '#' Rep1_TCF3_HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > Rep2_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz
```

**MACS2 D1 Peak calling**

The published ChIPseq data for TCF3::HLF in HAL-01 did not provide clear and abundant peaks to let us define primary peaks within our HiChIP data using FitHiChIP. We instead used the TCF3::HLF HiChIP data itself to generate primary aligment.

Assuming both the TCF3_HiChIP_fusedrep.sh and TCF3_HiChIP_singleRep.sh scripts are run, it is possible to adjust the MACS2 settings to improve MACS2 peak calling and for better confidence in the peaks.

Running the MACS2.sh script in the MACS2 conda environment will generate several macs2.narrowpeak files for both fused replicate samples and single replicate samples. These ared used for validation of high confidence peaks using the IDR tools that should also be installed.

Lastly the important output used for FitHiChIP loop calling are the Oracle annotated files.

**FitHiChIP Loop Calling**
Run FitHiChIP via bash script using docker image.
Make sure the configfile is properly setup and that all input files are on the specefied paths.

```
conda activate /home/ubuntu/HiCPro_ENV/
cd /home/ubuntu/FitHiChIP/
bash ./FitHiChIP_Docker.sh -C /home/ubuntu/HiChIP-Analysis/configfile_CB_TCF3_JoinedRep_5kb_50kb_3M
```
