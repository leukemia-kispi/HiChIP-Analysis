#!usr/bin/bash

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="/mnt/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="/mnt/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="/mnt/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for read before and after trimming, *fg.gz files
FASTQ_DIR="/mnt/1.RawData"
# Set output directories
OUTPUT_DIR_TRIM="/mnt/3.TRIM"
OUTPUT_CHIP_ALIGN="/mnt/4.ChIP_Alignment"
OUTPUT_CHIP_SUB="/mnt/4.ChIP_Alignment/Outputs"
BIGWIG_Coverage="/mnt/7.Deeptool_Matrix/Coverage/"


# Initialize Conda
eval "$(conda shell.bash hook)"

# Activate Conda Environment named DovetailHiChIP
CONDA_ENV="DovetailHiChIP"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

# Array containing replicate numbers found in filenames generated with TrimGalore.
NUMBERS=("1" "2") # Replace with your actual replicate numbers

# Flag to check if trimming needs to be performed initially set to false
perform_trimming=false

# Loop through each pair of FASTQ files if working with paired-end read files
for num in "${NUMBERS[@]}"; do
    # Check if trimmed files already exist for the current replicate
    if [ ! -f "$OUTPUT_DIR_TRIM/*rep${num}_R1_val_1.fq.gz" ] || [ ! -f "$OUTPUT_DIR_TRIM/*rep${num}_R2_val_2.fq.gz" ]; then
        perform_trimming=true
        break  # No need to check other replicates once one is found missing
    fi
done

# Perform trimming only if the flag is set to true
if [ "$perform_trimming" = true ]; then
    for num in "${NUMBERS[@]}"; do
        # Set path to input FASTQ files using wildcard pattern
        READ1="$FASTQ_DIR/*_rep${num}_R1.fastq.gz"
        READ2="$FASTQ_DIR/*_rep${num}_R2.fastq.gz"

        # Trim samples and generate new FastQC files for all replicates
        trim_galore --fastqc --phred33 --length 30 --output_dir $OUTPUT_DIR_TRIM -j 4 --paired $READ1 $READ2

        echo "Trimming for sample $num completed."
    done
else
    echo "Trimming not needed as output files already exist."
fi

# Set path to input FASTQ files using wildcard pattern
# BLF referse to black list filtered file
CHIP_R1="$OUTPUT_DIR_TRIM/*_R1_val_1.fq.gz"
CHIP_R2="$OUTPUT_DIR_TRIM/*_R2_val_2.fq.gz"
MAPPED_BAM="*_rep${num}_cle_sort.bam"
MAPPED_BLF_BAM="BLF_*_rep${num}_cle_sort.bam"
MAPPED_BLF_BAM_DUPFLAG="BLF_*_rep${num}_cle_sort_dupsflag.bam"
MAPPED_BLF_TXT_DUP="BLF_*_rep${num}_cle_sort_dups.txt"
MAPPED_BLF_BAM_DD="BLF_*_rep${num}_cle_sort_dd.bam"
MAPPED_MERGED_BLF_BAM_DD="BLF_*_merged_cle_sort_dd.bam"
BIGWIG_BLF_DD="BLF_*_rep${num}_cle_sort_dd.bw"
$BIGWIG_MERGED_BLF_DD="BLF_*_merged_cle_sort_dd.bw"

#Create sorted BAM files with grep to remove alignments to alternative contigs, unlocalized sequence, or unplaced sequence.#####################
bwa mem -5 -T25 -t32 $REF_FASTA $CHIP_R1 $CHIP_2 | samtools view -hS | grep -v chrUn | grep -v random | grep -v _alt | samtools view -bS -@32 | samtools sort -@32 -o $OUTPUT_CHIP_ALIGN/$MAPPED_BAM 


#Remove BlackList region & index
bedtools intersect -v -abam $OUTPUT_CHIP_ALIGN/$MAPPED_BAM -b $BLACKLIST > $OUTPUT_CHIP_ALIGN/$MAPPED_BLF_BAM
samtools index $OUTPUT_CHIP_ALIGN/$MAPPED_BLF_BAM


# Activate Conda Environment named Picard
CONDA_ENV="Picard"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

#Picard MarkDuplicates and Remove
java -jar /home/ubuntu/miniconda3/envs/Picard/share/picard-2.25.7-0/picard.jar MarkDuplicates -I $OUTPUT_CHIP_ALIGN/$MAPPED_BLF_BAM -O $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DUPFLAG -M $OUTPUT_CHIP_SUB/$MAPPED_BLF_TXT_DUP --REMOVE_DUPLICATES false



# Activate Conda Environment named DovetailHiChIP
CONDA_ENV="ChIP"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

#Generate final deduplicated BAM file and indexing
samtools view -b -F 1024 $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DUPFLAG > $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD
samtools index $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD

###MERGE Replicate Files################
samtools merge -f $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD_REP1 $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD_REP2   ####have to update file annotations 
samtools index $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD

# Activate Conda Environment named deeptool
CONDA_ENV="deeptool"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

#Enrichment for IGV
bamCoverage -b $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD -o $BIGWIG_Coverage/$BIGWIG_BLF_DD --effectiveGenomeSize 2913022398 -bl $BLACKLIST--normalizeUsing RPKM -p max -bs 10 \
--extendReads --ignoreForNormalization M

bamCoverage -b  $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD -o $BIGWIG_Coverage/$BIGWIG_MERGED_BLF_DD --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 \
--extendReads --ignoreForNormalization M

# Activate Conda Environment named macs2
CONDA_ENV="macs2"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

#Call Peaks with permissive settings.
macs2 callpeak -t BLF_ChIP_HAL01_H3K27ac_Rep1_cle_sort_dd.bam -c BLF_ChIP_HAL01_control_Input_Rep1_cle_sort_dd.bam -f BAMPE -g hs -p 0.05 -B --outdir /mnt/RawChIP_Histon/MACS2 -n BLF_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2


macs2 callpeak -t BLF_ChIP_HAL01_H3K27ac_Rep2_cle_sort_dd.bam -c BLF_ChIP_HAL01_control_Input_Rep2_cle_sort_dd.bam -f BAMPE -g hs -p 0.05 -B --outdir /mnt/RawChIP_Histon/MACS2 -n BLF_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05macs2


#Enter Merged Folder
macs2 callpeak -t  $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD -c  $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD_INPUT -f BAMPE -g hs -p 0.000000001 -B --outdir /mnt/RawChIP_Histon/MACS2 -n HAL01_H3K27ac_merged_cle_sort_dd_p9macs2

macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN --keep-dup 10 --min-length 300 -p 0.000000001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2/$MACS2_JoinedRep_Oracle


#Sort Peak files

sort -k8,8nr BLF_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak > /mnt/RawChIP_Histon/IDR/Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak


#IDR without Oracle Peak file
idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05macs2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file HAL-01_H3K27ac_cle-idr --plot --log-output-file HAL-01_H3K27ac_cle.idr.log

idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05bmacs2_peaks.broadPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05bmacs2_peaks.broadPeak  --input-file-type broadPeak --rank p.value --output-file HAL-01_H3K27ac_cle_broad-idr --plot --log-output-file HAL-01_H3K27ac_cle_broad.idr.log
	                                                                     

#IDR with Oracle Peak file
idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05macs2_peaks.narrowPeak --peak-list Sort_HAL01_H3K27ac_merged_cle_sort_dd_p9macs2_peaks.narrowPeak  --input-file-type narrowPeak --rank p.value --output-file Oracle_HAL-01_H3K27ac_cle-idr --plot --log-output-file Oracle_HAL-01_H3K27ac_cle.idr.log

idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05bmacs2_peaks.broadPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05bmacs2_peaks.broadPeak --peak-list Sort_HAL01_H3K27ac_merged_cle_sort_dd_p9bmacs2_peaks.broadPeak  --input-file-type broadPeak --rank p.value --output-file Oracle_HAL-01_H3K27ac_cle_broad-idr --plot --log-output-file Oracle_HAL-01_H3K27ac_cle_broad.idr.log