# FitHiChIP Loop calling

FitHiChIP will be run from Docker image. If instructions in the [README.md](https://github.com/ValdemarP267/HiChIP-Analysis) will guide you through docker installation. This ensure to run FitHiChIP without having to install all dependecies from scratch.

Select the directory to contain the FitHiChIP source code, and clone it

```
conda activate FitHiChIP
cd /home/ubutu
git clone https://github.com/ay-lab/FitHiChIP.git
sudo chmod 777 -R FithHiChIP
```

With the Outputs from above you will need:

- The Pairs files converted to HiC-Pro format
- MACS2 called peaks from relevant ChIP-seq data or do MACS2 call peaks from primary algimnents in HiChIP data.
- Config file(example provided in this repository) specefying file locations and parameters

**Filter pairs**

```
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' JoinedRep_TCF3_HLF_hg38_nodd_hicpro_mapped.pairs -o JoinedRep_TCF3-HLF_hg38_nodd_mapped.filtered.pairs
```

**HiCPro Valid Pairs Files**

```
grep -v '#' JoinedRep_TCF3-HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > JoinedRep_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz
```

**FitHiChIP Loop Calling**
Run FitHiChIP via bash script using docker image.

```
cd /home/ubuntu/FitHiChIP/
bash ./FitHiChIP_Docker.sh -C /home/ubuntu/HiChIP-Analysis/configfile_CB_TCF3_JoinedRep_5kb_Merge_50kb_3M
```
