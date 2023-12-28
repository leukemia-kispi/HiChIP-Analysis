#!/bin/bash

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="/mnt/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="/mnt/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="/mnt/0.BlackList/hg38-blacklist.v2.bed "
# Set Path for read *fg.gz files
FASTQ_DIR="/mnt/1.RawData"
# Set output directories
OUTPUT_DIR_TRIM="/mnt/3.TRIM"
OUTPUT_HICHIP_ALIGN="/mnt/4.HiChIP_Alignment"
OUTPUT_HICHIP_SUB="/mnt/4.HiChIP_Alignment/"
# Thread usage
cores= 32
#Thread usage for pairtools dedup and split processes
cores2= 16
# Set Path to temporary directory
TEMP="/mnt/tmp"

create_directory "/mnt/5.MACS2"
create_directory "/mnt/6.FitHiChIP_Output"

# Initialize Conda
eval "$(conda shell.bash hook)"

#TrimGalor.
conda activate TRIM

# Array containing your replicates found in file names generated with TrimGalore. Expected filename structure *_rep${num}_1.fastq.gz.
NUMBERS=("1" "2") #replace  with your replicate numbers

#Loop through each pair of FASTQ files if working with paire-end read files
for num in "${NUMBERS[@]}"; do
	# Replace NUM_SAMPLES i<= with actual number of samples
    # Set path to input FASTQ files using wildcard pattern
	READ1="$FASTQ_DIR"/"*_rep${num}_R1.fastq.gz"
	READ2="$FASTQ_DIR"/"*_rep${num}_R2.fastq.gz"

# Trim samples and generate new fastqc files for all replicates
trim_galore --fastqc --phred33 --length 50 --output_dir $OUTPUT_DIR_TRIM -j 4 --paired $READ1 $READ2

 echo "Trimming for sample $num completed."
done

echo "All Trimming completed."

cd $OUTPUT_DIR_TRIM
multiqc .

#DovetailHiChIP
conda activate DovetailHiChIP

#Fuse Fasta files
# Concatenate R1 fastq files
cat "$OUTPUT_DIR_TRIM"/*_R1_val_1.fq.gz > JoinedFastq_R1.fq.gz

# Concatenate R2 fastq files
cat "$OUTPUT_DIR_TRIM"/*_R2_val_2.fq.gz > JoinedFastq_R2.fq.gz

# Alignment
cd $OUTPUT_HICHIP_ALIGN
HiChiP_R1="/mnt/3.TRIM/JoinedFastq_R1.fq.gz"
HiChiP_R2="/mnt/3.TRIM/JoinedFastq_R1.fq.gz"
MAPPED_PAIRS="JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs"
MAPPED_BAM="JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam"

bwa mem -5SP -T0 -t$cores $REF_FASTA $HiChiP_R1 $HiChiP_R2 | \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores2 --nproc-out $cores2 --chroms-path $REF_GENOME | \
pairtools sort --tmpdir=$TEMP--nproc $cores|\
#pairtools dedup --nproc-in $cores2 --nproc-out $cores2 --mark-dups --output-stats <stats.txt>|\
pairtools split --nproc-in $cores2 --nproc-out $cores2 --output-pairs $MAPPED_PAIRS --output-sam -|\
samtools view -bS -@$cores | \
samtools sort -@$cores -o $MAPPED_BAM;samtools index $MAPPED_BAM

echo "HiCHIP Aligmnent Complete"

#QC compare ChIP-seq TCF3-HLF_FLAG
bash enrichment_stats.sh -g /mnt/GRCh38_assembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -b /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -p /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_cle-idr.bed -t 16 -x HiChIP-TCF3-HLF
bash enrichment_stats.sh -g /mnt/GRCh38_assembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -b /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -p /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -t 16 -x HiChIP-TCF3-HLF_bw175

#QC Plot ChIP-seq TCF3-HLF_FLAG
python3 plot_chip_enrichment_bed.py -bam /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -peaks /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_cle-idr.bed -output HiChIP_TCF3_HLF_ChIP_FLAG_enrichment.png
python3 plot_chip_enrichment_bed.py -bam /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -peaks /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -output HiChIP_TCF3_HLF_ChIP_FLAG_bw175_enrichment.png

#Enrichment for IGV
conda activate deeptool
bamCoverage -b JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -o BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.bw --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

#ContacMaps
conda activate Juicebox
java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/My_HiC-Pro/juicer_tools_1.22.01.jar pre --threads 32 JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs JoinedRep_TCF3-HLF_HAL01_hg38_nodd_contact_map.hic /mnt/GRCh38_assembly/GRCh38_no_alt_ref.genome

#call 1D peaks with MACS2
conda activate MACS2

#Remove black list
bedtools intersect -v -abam JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -b /mnt/BlastListPeaks/hg38-blacklist.v2.bed > BLF_JoinedRep_TCF3-HLF_HAL01_hg38_nodd_mapped.PT.bam

#MACS2 Peak calling on BlackList filtered BAM.
samtools view -b BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -h -F 0x900 | bedtools bamtobed -i stdin > BLF_JoinedRep_TCF3_HLF_hg38_nodd_primary.aln.bed

macs2 callpeak -t BLF_JoinedRep_TCF3_HLF_hg38_nodd_primary.aln.bed -p 0.000000001 -n BLF_JoinedRep_TCF3_HLF_hg38_nodd_p9.macs2
macs2 callpeak -t BLF_JoinedRep_TCF3-HLF_HAL01_hg38_nodd_primary.aln.bed -p 0.000000001 -g 2913022398 -n BLF_JoinedRep_TCF3_HLF_gs_hg38_nodd_p9.macs2

#Oracle File for IDR and IDR
macs2 callpeak -t BLF_JoinedRep_TCF3-HLF_HAL01_hg38_nodd_primary.aln.bed --keep-dup 10 --min-length 300 -p 0.000000001 --bw 300 --mfold 5 50 -n BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2
sort -k8,8nr  BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2_peaks.narrowPeak > /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Oracle_MACS2/Sort_BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2_peaks.narrowPeak 
sort -k8,8nr  BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak > /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak
sort -k8,8nr  BLF_TCF3-HLF_HAL01_rep2_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak > /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep2_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak

idr --samples /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep2_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak --peak-list /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Oracle_MACS2/Sort_BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file Oracle_HAL-01_TCF3_JoinedRep_cle-idr --plot --log-output-file Oracle_HAL-01_TCF3_JoinedRep_cle.idr.log

#FitHiChIP Loop calling
conda activate /home/ubuntu/My_HiCPro_ENV/

export PATH="/home/ubuntu/My_HiC-Pro/HiC-Pro_3.0.0/bin/":$PATH

BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak
#FitHiChIP Loop Calling
#Filter pairs
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' JoinedRep_TCF3_HLF_hg38_nodd_hicpro_mapped.pairs -o JoinedRep_TCF3-HLF_hg38_nodd_mapped.filtered.pairs

#HiCPro Valid Pairs Files
grep -v '#' JoinedRep_TCF3-HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > JoinedRep_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz

cd /home/ubuntu/My_HiC-Pro/FitHiChIP/
bash ./FitHiChIP_Docker.sh -C /home/ubuntu/My_HiC-Pro/FitHiChIP/configfile_CB_TCF3_JoinedRep_5kb_Merge_50kb_3M