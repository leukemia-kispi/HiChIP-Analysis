#!/bin/bash

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="/mnt/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="/mnt/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="/mnt/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for read before and after trimming, *fg.gz files
FASTQ_DIR="/mnt/1.RawData/H3K27ac"
# Set output directories
OUTPUT_DIR_TRIM="/mnt/3.TRIM/H3K27ac"
OUTPUT_HICHIP_ALIGN="/mnt/4.HiChIP_Alignment/H3K27ac"
OUTPUT_HICHIP_SUB="/mnt/4.HiChIP_Alignment/H3K27ac/Outputs"
# Thread usage
cores=32
#Thread usage for pairtools dedup and split processes
cores2=16
# Set Path to temporary directory
TEMP="/mnt/tmp"

# Initialize Conda
eval "$(conda shell.bash hook)"

# Activate Conda Environment named DovetailHiChIP
CONDA_ENV="DovetailHiChIP"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

# Set path to input FASTQ files using wildcard pattern
#READ1="$FASTQ_DIR/*_R1.fastq.gz"
#READ2="$FASTQ_DIR/*_R2.fastq.gz"

# Trim samples and generate new FastQC files for all replicates
#trim_galore --fastqc --phred33 --length 50 --output_dir $OUTPUT_DIR_TRIM -j 4 --paired $READ1 $READ2

#echo "Trimming for sample H3K27ac completed."

# Alignment Output directory
cd $OUTPUT_HICHIP_ALIGN

# Set path to input FASTQ files using wildcard pattern
HICHIP_R1="$OUTPUT_DIR_TRIM/*_R1_val_1.fq.gz"
HICHIP_R2="$OUTPUT_DIR_TRIM/*_R2_val_2.fq.gz"
MAPPED_PAIRS="_H3K27ac_hg38_dd_mapped.pairs"
MAPPED_BAM="_H3K27ac_hg38_dd_mapped.PT.bam"
MAPPED_BLF_BAM="BLF_H3K27ac_hg38_dd_mapped.PT.bam"

bwa mem -5SP -T0 -t$cores $REF_FASTA $HICHIP_R1 $HICHIP_R2 | \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores2 --nproc-out $cores2 --chroms-path $REF_GENOME | \
pairtools sort --tmpdir=$TEMP --nproc $cores | \
pairtools dedup --nproc-in $cores2 --nproc-out $cores2 --mark-dups --output-stats H3K27ac_stats.txt | \
pairtools split --nproc-in $cores2 --nproc-out $cores2 --output-pairs $MAPPED_PAIRS --output-sam -|\
samtools view -bS -@$cores | \
samtools sort -@$cores -o $MAPPED_BAM;samtools index $MAPPED_BAM

echo "HiCHIP Aligmnent Complete for H3K27ac"

cd /home/ubuntu

#Remove black list
bedtools intersect -v -abam $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -b $BLACKLIST > $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM

#QC compare ChIP-seq TCF3-HLF_FLAG
bash /home/ubuntu/HiChiP/enrichment_stats.sh -g $REF_FASTA -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -p /home/ubuntu/HiChIP_Analysis/ChIP-Seq/ChIP_Seq_HAL01_H3K27ac_merged_cle_dd_q0.01macs2_peaks.narrowPeak.bed -t $cores2 -x $OUTPUT_HICHIP_SUB/HiChIP_H3K27ac

#QC Plot ChIP-seq TCF3-HLF_FLAG
python3 /home/ubuntu/HiChiP/plot_chip_enrichment_bed.py -bam $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -peaks /home/ubuntu/HiChIP_Analysis/ChIP-Seq/ChIP_Seq_HAL01_H3K27ac_merged_cle_dd_q0.01macs2_peaks.narrowPeak.bed -output $OUTPUT_HICHIP_SUB/HiChIP_H3K27ac_enrichment.png

echo "HiCHIP Aligmnent QC Complete for H3K27ac"

#Enrichment for IGV
bamCoverage -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -o $OUTPUT_HICHIP_SUB/BLF_H3K27ac_hg38_dd_mapped.bw --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

echo "Generated Bigwig file Complete for H3K27ac"

#ContacMaps
java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/HiChiP/juicer_tools_1.22.01.jar pre --threads $cores $OUTPUT_HICHIP_ALIGN/$MAPPED_PAIRS $OUTPUT_HICHIP_SUB/H3K27ac_HAL01_hg38_dd_contact_map.hic $REF_GENOME

echo "Generated .hic file Complete H3K27ac"
