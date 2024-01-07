# FitHiChIP Loop calling

FitHiChIP will be run from Docker image. If instructions in the [README.md](https://github.com/ValdemarP267/HiChIP-Analysis) will guide you through docker installation. This ensure to run FitHiChIP without having to install all dependecies from scratch.

Select the directory to contain the FitHiChIP source code, and clone it

```
conda activate FitHiChIP
cd /home/ubuntu
git clone https://github.com/ay-lab/FitHiChIP.git
sudo chmod 777 -R FithHiChIP
```

With the Outputs from above you will need:

- The Pairs files converted to HiC-Pro format.
- MACS2 called peaks from relevant ChIP-seq data or do MACS2 call peaks from primary algimnents in HiChIP data
- Config file(example provided in this repository) specefying file locations and parameters

**Filter pairs**

Using partools select in the DovtailHiChIP conda environment to filter and include unique and rescue mapped read pairs.

```
conda activate DovetailHiChIP
cd /mnt/4.HiChIP_Alignment
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs -o JoinedRep_TCF3_HLF_hg38_nodd_mapped.filtered.pairs
```

Same for single replicate run

```
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' Rep1_TCF3_HLF_hg38_nodd_mapped.pairs -o Rep1_TCF3_HLF_hg38_nodd_mapped.filtered.pairs
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' Rep2_TCF3_HLF_hg38_nodd_mapped.pairs -o Rep2_TCF3_HLF_hg38_nodd_mapped.filtered.pairs
```

**HiCPro Valid Pairs Files**

The paired files generated during Aligment needs to be converted into valid HiC-Pro files to proceed  with FitHiChIP this can be done with the following commands.

```
conda activate DovetailHiChIP
cd /mnt/4.HiChIP_Alignment
grep -v '#' JoinedRep_TCF3_HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > JoinedRep_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz
```

Same for single replicate run

```
grep -v '#' Rep1_TCF3_HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > Rep1_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz
grep -v '#' Rep1_TCF3_HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > Rep2_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz
```

**MACS2 D1 Peak calling**

As we do not have good ChIPseq data for TCF3::HLF in HAL-01 that lets us define primary peaks within our HiChIP data using FitHiChIP. We instead have to you the HiChIP data itself to generate primary aligment.

Assuming you run both the TCF3_HiChIP_fusedrep.sh and TCF3_HiChIP_singleRep.sh scripts. We can take it a step further and adjust the MACS2 settings to improve MACS2 peak calling and improve our confidence in the peaks.

Running the MACS2.sh script in the MACS2 conda environment will generate several macs2.narrowpeak files for both fused replicate samples and single replicate samples. These ared used for validation of high confidence peaks by doing the idr.

Lastly the important output is the Oracle file that is to be used for FitHiChIP calling. When testing the file showed better representation for the placement of relevant peaks as observed when overlayed with bigwig files.

**FitHiChIP Loop Calling**
Run FitHiChIP via bash script using docker image.
Make sure the configfile is properly setup and that all input files are on the specefied paths.

```
cd /home/ubuntu/FitHiChIP/
bash ./FitHiChIP_Docker.sh -C /home/ubuntu/HiChIP-Analysis/configfile_CB_TCF3_JoinedRep_5kb_Merge_50kb_3M
```
