#!/usr/bin/env bash

#############################
### MAIN FILE PATHS #########
#############################

# Prompt the user for the main directory
read -rp "Enter the path to the MAIN_DIR: " MAIN_DIR

# Check if input is empty
if [ -z "$MAIN_DIR" ]; then
    echo "Error: No directory path entered."
    exit 1
fi

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="/mnt/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="/mnt/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="/mnt/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for reads *fg.gz files
FASTQ_DIR="/mnt/1.RawData"
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



################################################
### ADJUST THESE TO MATCH SAMPLE NAMES #########
################################################

# Array containing Cell lines, replicate numbers and conditions found in filenames and defining samples.
# Expected Sample nomenclature follows this pattern HiChIP_<CellLine>_<conditions>_Rep<NUMBERS>_suffix.
CellLine=("HAL01") #Replace with your actual Cell Line 
conditions=("TCF3") #Replace with your actual conditions
NUMBERS=("1" "2") # Replace with your actual replicate numbers


##########################################################
### INITIALIZE DOVETAILHICHIP CONDA ENVIRONMENT #########
##########################################################

#Environment should have following installed: bwa-mem, samtools, pairtools, fastqc, bedtools, multiqc, trim-galore, deeptools

# Initialize Conda
eval "$(conda shell.bash hook)"

# Echo current Conda environment
echo "Current Conda environment: $CONDA_DEFAULT_ENV"

# Only activate if not already in the target environment
CONDA_ENV="DovetailHiChIP"
if [[ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV" ]]; then
    echo "Activating Conda environment: $CONDA_ENV"
    conda activate "$CONDA_ENV"
else
    echo "Already in the target Conda environment."
fi




######################
### TRIMMING #########
######################

# Flag to check if trimming needs to be performed initially set to false
perform_trimming=false

# Loop through each pair of FASTQ files if working with paired-end read files
for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do
        for num in "${NUMBERS[@]}"; do
            # Check if trimmed files already exist for the current replicate
            if [ ! -f "$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep${num}_R1_val_1.fq.gz" ] || [ ! -f "$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep${num}_R2_val_2.fq.gz" ]; then
                perform_trimming=true
                break  # No need to check other replicates once one is found missing
            fi
        done
    done
done

# Perform trimming only if the flag is set to true
if [ "$perform_trimming" = true ]; then
    for cell in "${CellLine[@]}"; do 
        for cond in "${conditions[@]}"; do
            for num in "${NUMBERS[@]}"; do
                # Set path to input FASTQ files using wildcard pattern
                READ1="$FASTQ_DIR/HiChIP_${cell}_${cond}_Rep${num}_R1.fastq.gz"
                READ2="$FASTQ_DIR/HiChIP_${cell}_${cond}_Rep${num}_R2.fastq.gz"
        
                echo "Trimming pair $READ1 and $READ2"
            
                # Trim samples and generate new FastQC files for all replicates
                trim_galore --fastqc --phred33 --length 30 --output_dir $OUTPUT_DIR_TRIM -j 4 --paired $READ1 $READ2

                echo "Trimming for sample $num completed."
            done
        done
    done
else
    echo "Trimming not needed as output files already exist."
fi

#NOTE FOR VALIDIP Consider adding MULTIQc step after Trimming with promt to possibly skip this step#######

#Fuse Fasta files
# Concatenate R1 fastq files
cat $OUTPUT_DIR_TRIM/*rep1_R1_val_1.fq.gz $OUTPUT_DIR_TRIM/*rep2_R1_val_1.fq.gz > JoinedFastq_R1.fq.gz

# Concatenate R2 fastq files
cat $OUTPUT_DIR_TRIM/*rep1_R2_val_2.fq.gz $OUTPUT_DIR_TRIM/*rep2_R2_val_2.fq.gz > JoinedFastq_R2.fq.gz

#Give permission to new files
sudo chmod 777 $OUTPUT_DIR_TRIM/JoinedFastq_R1.fq.gz
sudo chmod 777 $OUTPUT_DIR_TRIM/JoinedFastq_R2.fq.gz

echo "Fusion of FASTA replicates done"

# Alignment Output directory
cd $OUTPUT_HICHIP_ALIGN

# Input/Output files for alignment
# BLF referse to black list filtered file
HICHIP_R1="$OUTPUT_DIR_TRIM/JoinedFastq_R1.fq.gz"
HICHIP_R2="$OUTPUT_DIR_TRIM/JoinedFastq_R2.fq.gz"
MAPPED_PAIRS="JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs"
MAPPED_BAM="JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam"
MAPPED_BLF_BAM="BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam"

# Alignment, dedup skipped
bwa mem -5SP -T0 -t$cores $REF_FASTA $HICHIP_R1 $HICHIP_R2 | \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores2 --nproc-out $cores2 --chroms-path $REF_GENOME | \
pairtools sort --tmpdir=$TEMP --nproc $cores | \
#pairtools dedup --nproc-in $cores2 --nproc-out $cores2 --mark-dups --dry-run --output-stats JoinedRep_stats.txt | \
pairtools split --nproc-in $cores2 --nproc-out $cores2 --output-pairs $MAPPED_PAIRS --output-sam -|\
samtools view -bS -@$cores | \
samtools sort -@$cores -o $MAPPED_BAM;samtools index $MAPPED_BAM

echo "HiCHIP Aligmnent Complete"

cd /home/ubuntu

#Remove black list
bedtools intersect -v -abam $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -b $BLACKLIST > $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM
samtools index $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM

#QC compare ChIP-seq TCF3-HLF_FLAG
bash /home/ubuntu/HiChiP/enrichment_stats.sh -g $REF_FASTA -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -p /home/ubuntu/HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3HLF_FLAG_bw175_cle-idr.bed -t $cores2 -x $OUTPUT_HICHIP_SUB/HiChIPJoinedFastq-TCF3HLF_bw175

#QC Plot ChIP-seq TCF3-HLF_FLAG
python3 /home/ubuntu/HiChiP/plot_chip_enrichment_bed.py -bam $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -peaks /home/ubuntu/HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3HLF_FLAG_bw175_cle-idr.bed -output $OUTPUT_HICHIP_SUB/HiChIPJoinedFastq_TCF3HLF_ChIP_FLAG_bw175_enrichment.png

echo "HiCHIP Aligmnent QC Complete"

#Enrichment for IGV
bamCoverage -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -o $OUTPUT_HICHIP_SUB/BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.bw --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

echo "Generated Bigwig file Complete"

#ContacMaps
java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/HiChiP/juicer_tools_1.22.01.jar pre --threads $cores $OUTPUT_HICHIP_ALIGN/$MAPPED_PAIRS $OUTPUT_HICHIP_SUB/JoinedRep_TCF3HLF_HAL01_hg38_nodd_contact_map.hic $REF_GENOME

echo "Generated .hic file Complete"
