#!/bin/bash

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="/mnt/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="/mnt/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="/mnt/0.BlackList/hg38-blacklist.v2.bed "
# Set output directories
OUTPUT_HICHIP_ALIGN="/mnt/4.HiChIP_Alignment"
OUTPUT_HICHIP_SUB="/mnt/4.HiChIP_Alignment/Outputs"
OUTPUT_MACS2="/mnt/5.MACS2"
OUTPUT_MACS2_SORT="/mnt/5.MACS2/SORT"
OUTPUT_MACS2_Permissive="/mnt/5.MACS2/Permissive"
OUTPUT_MACS2_IDR="/mnt/5.MACS2/IDR"
#Inputs
MAPPED_BAM="JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam"
MAPPED_BAM_Rep1="Rep1_TCF3_HLF_hg38_nodd_mapped.PT.bam"
MAPPED_BAM_Rep2="Rep2_TCF3_HLF_hg38_nodd_mapped.PT.bam"
MAPPED_BLF_BAM="BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam"
MAPPED_BLF_BAM_Rep1="BLF_Rep1_TCF3_HLF_hg38_nodd_mapped.PT.bam"
MAPPED_BLF_BAM_Rep2="BLF_Rep2_TCF3_HLF_hg38_nodd_mapped.PT.bam"
#Outputs
PRIMARY_ALN="BLF_JoinedRep_TCF3_HLF_hg38_nodd_primary.aln.bed"
PRIMARY_ALN_Rep1="BLF_Rep1_TCF3_HLF_hg38_nodd_primary.aln.bed"
PRIMARY_ALN_Rep2="BLF_Rep2_TCF3_HLF_hg38_nodd_primary.aln.bed"
MACS2_JoinedRep="BLF_JoinedRep_TCF3_HLF_gs_hg38_nodd_p9.macs2"
MACS2_JoinedRep_Oracle="BLF_JoinedRep_TCF3_HLF_Oracle.macs2"
MACS2_JoinedRep_SORT="BLF_JoinedRep_TCF3_HLF_OracleSort.macs2"
MACS2_Rep1="BLF_Rep1_TCF3-HLF_gs_hg38_p9.macs2"
MACS2_Rep2="BLF_Rep2_TCF3-HLF_gs_hg38_p9.macs2"
MACS2_Rep1_Permissive="BLF_Rep1_TCF3_HLF_Permissive.macs2"
MACS2_Rep2_Permissive="BLF_Rep2_TCF3_HLF_Permissive.macs2"
MACS2_Rep1_SORT="BLF_Rep1_TCF3_HLF_PermissiveSort.macs2"
MACS2_Rep2_SORT="BLF_Rep2_TCF3_HLF_PermissiveSort.macs2"
MACS2_Rep1_OracleSet="BLF_Rep1_TCF3_HLF_OracleSet.macs2"
MACS2_Rep2_OracleSet="BLF_Rep2_TCF3_HLF_OracleSet.macs2"

# Initialize Conda
eval "$(conda shell.bash hook)"

# Activate Conda Environment named MACS2
CONDA_ENV="MACS2"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

#Remove black list.
#bedtools intersect -v -abam $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -b $BLACKLIST > $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM
#bedtools intersect -v -abam $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM_Rep1 -b $BLACKLIST > $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM_Rep1
#bedtools intersect -v -abam $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM_Rep2 -b $BLACKLIST > $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM_Rep2

echo "Black list regions removed"

#MACS2 Peak calling on BlackList filtered BAM.
samtools view -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -h -F 0x900 | bedtools bamtobed -i stdin > $OUTPUT_HICHIP_SUB/$PRIMARY_ALN
samtools view -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM_Rep1 -h -F 0x900 | bedtools bamtobed -i stdin > $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep1
samtools view -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM_Rep2 -h -F 0x900 | bedtools bamtobed -i stdin > $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep2

echo "Primary alignment .bed file generated"

macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN -p 0.000000001 -g 2913022398 -n $OUTPUT_MACS2/$MACS2_JoinedRep
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep1 -p 0.000000001 -g 2913022398 -n $OUTPUT_MACS2/$MACS2_Rep1
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep2 -p 0.000000001 -g 2913022398 -n $OUTPUT_MACS2/$MACS2_Rep2

echo "Standard MACS2 peak calling done"

#Oracle File for IDR 
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN --keep-dup 10 --min-length 300 -p 0.000000001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2/$MACS2_JoinedRep_Oracle
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep1 --keep-dup 10 --min-length 300 -p 0.000000001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2/$MACS2_Rep1_OracleSet
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep2 --keep-dup 10 --min-length 300 -p 0.000000001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2/$MACS2_Rep2_OracleSet
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep1 --keep-dup 10 --min-length 300 -p 0.00001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2_Permissive/$MACS2_Rep1_Permissive
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep2 --keep-dup 10 --min-length 300 -p 0.00001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2_Permissive/$MACS2_Rep2_Permissive
sudo chmod 777 -R $OUTPUT_MACS2
sort -k8,8nr $OUTPUT_MACS2/$MACS2_JoinedRep_Oracle > $OUTPUT_MACS2_SORT/$MACS2_JoinedRep_SORT
sort -k8,8nr $OUTPUT_MACS2_Permissive/$MACS2_Rep1_Permissive > $OUTPUT_MACS2_SORT/$MACS2_Rep1_SORT
sort -k8,8nr $OUTPUT_MACS2_Permissive/$MACS2_Rep2_Permissive > $OUTPUT_MACS2_SORT/$MACS2_Rep2_SORT

idr --samples $OUTPUT_MACS2_SORT/$MACS2_Rep1_SORT $OUTPUT_MACS2_SORT/$MACS2_Rep2_SORT --peak-list $OUTPUT_MACS2_SORT/$MACS2_JoinedRep_SORT --input-file-type narrowPeak --rank p.value --output-file $OUTPUT_MACS2_IDR/OraclePeaks_HAL01_TCF3_cle-idr --plot --log-output-file $OUTPUT_MACS2_IDR/OraclePeak_HAL01_TCF3_cle.idr.log

echo "All MACS2 Oracle runs Complete"