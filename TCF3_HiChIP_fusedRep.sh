#!/bin/bash

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="/mnt/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="/mnt/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="/mnt/0.BlackList/hg38-blacklist.v2.bed "
# Set Path for read before and after trimming, *fg.gz files
FASTQ_DIR="/mnt/1.RawData"
HiChIP_R1="/mnt/3.TRIM/JoinedFastq_R1.fq.gz"
HiChIP_R2="/mnt/3.TRIM/JoinedFastq_R1.fq.gz"
# Set output directories
OUTPUT_DIR_TRIM="/mnt/3.TRIM"
OUTPUT_HICHIP_ALIGN="/mnt/4.HiChIP_Alignment"
OUTPUT_HICHIP_SUB="/mnt/4.HiChIP_Alignment/Outputs"
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
        trim_galore --fastqc --phred33 --length 50 --output_dir $OUTPUT_DIR_TRIM -j 4 --paired $READ1 $READ2

        echo "Trimming for sample $num completed."
    done
else
    echo "Trimming not needed as output files already exist."
fi

#Fuse Fasta files
# Concatenate R1 fastq files
cat $OUTPUT_DIR_TRIM/*_R1_val_1.fq.gz > JoinedFastq_R1.fq.gz

# Concatenate R2 fastq files
cat $OUTPUT_DIR_TRIM/*_R2_val_2.fq.gz > JoinedFastq_R2.fq.gz

# Alignment Output directory
cd $OUTPUT_HICHIP_ALIGN
MAPPED_PAIRS="JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs"
MAPPED_BAM="JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam"

# Alignment, dedup skipped
bwa mem -5SP -T0 -t$cores $REF_FASTA $HICHIP_R1 $HICHIP_R2 | \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores2 --nproc-out $cores2 --chroms-path $REF_GENOME | \
pairtools sort --tmpdir=$TEMP --nproc $cores | \
#pairtools dedup --nproc-in $cores2 --nproc-out $cores2 --mark-dups --dry-run --output-stats JoinedRep_stats.txt | \
pairtools split --nproc-in $cores2 --nproc-out $cores2 --output-pairs $MAPPED_PAIRS --output-sam -|\
samtools view -bS -@$cores | \
samtools sort -@$cores -o $MAPPED_BAM;samtools index $MAPPED_BAM

echo "HiCHIP Aligmnent Complete"

#QC compare ChIP-seq TCF3-HLF_FLAG
bash ./HiChiP/enrichment_stats.sh -g $REF_FASTA -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM  -p ./HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -t $cores2 -x $OUTPUT_HICHIP_SUB/HiChIPJoinedFastq-TCF3-HLF_bw175

#QC Plot ChIP-seq TCF3-HLF_FLAG
python3 ./HiChiP/plot_chip_enrichment_bed.py -bam $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -peaks ./HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -output $OUTPUT_HICHIP_SUB/HiChIPJoinedFastq_TCF3_HLF_ChIP_FLAG_bw175_enrichment.png

echo "HiCHIP Aligmnent QC Complete"

#Enrichment for IGV
bamCoverage -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -o $OUTPUT_HICHIP_SUB/BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.bw --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

echo "Generated Bigwig file Complete"

#ContacMaps
java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/HiChiP/juicer_tools_1.22.01.jar pre --threads $cores $OUTPUT_HICHIP_ALIGN/JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs $OUTPUT_HICHIP_SUB/JoinedRep_TCF3-HLF_HAL01_hg38_nodd_contact_map.hic $REF_GENOME

echo "Generated .hic file Complete"
