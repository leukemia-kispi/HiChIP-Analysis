#!/usr/bin/env bash
set -e
shopt -s nullglob # make globbing return empty array if no match

#Using this script assumes the script DirectoryArchitecture&CondaEnv.sh was run beforhand to create the DirectoryArchitecture and generate conda environments with needed tools.

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
REF_FASTA="$MAIN_DIR/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="$MAIN_DIR/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="$MAIN_DIR/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for reads *fg.gz files
FASTQ_DIR="$MAIN_DIR/1.RawData/HiChIP"
FASTQC="$MAIN_DIR/2.FASTQC/HiChIP"
# Set output directories
OUTPUT_DIR_TRIM="$MAIN_DIR/3.TRIM/HiChIP"
OUTPUT_HICHIP_ALIGN="$MAIN_DIR/4.Alignment/HiChIP"
OUTPUT_HICHIP_SUB="$MAIN_DIR/4.Alignment/HiChIP/Outputs"
# Thread usage
cores=32
#Thread usage for pairtools dedup and split processes
cores2=16
# Set Path to temporary directory
TEMP="$MAIN_DIR/tmp"

################################################
### ADJUST THESE TO MATCH SAMPLE NAMES #########
################################################

# Array containing Cell lines, replicate numbers and conditions found in filenames and defining samples.
# Expected Sample nomenclature follows this pattern HiChIP_<CellLine>_<conditions>_Rep<NUMBERS>_suffix.
CellLine=("HAL01") #Replace with your actual Cell Line 
conditions=("TCF3HLF") #Replace with your actual conditions
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

#############################
### MULITQC #################
#############################
echo
echo "Start MULTIQC summary of FASTQC output for trimmed fastq files"
# Promt to procceed or skip IDR
read -rp "Do you want to proceed MULTIQC Analysis (y/n): " confirm
if [[ "$confirm" == "y" ]]; then  
    multiqc $OUTPUT_DIR_TRIM
    echo "MULTIQC done"
else    
    echo " MULTIQC skipped"
fi

####################################
### MERGING OF FASTQ FILES #########
####################################

#Fastq files merged before alignment. Improve depth and downstream analysis.

# Concatenate R1 fastq files
for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do
        cat $OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep1_R1_val_1.fq.gz $OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep2_R1_val_1.fq.gz > $OUTPUT_DIR_TRIM/MergeFastq_${cell}_${cond}_R1.fq.gz
    done
done

# Concatenate R2 fastq files
for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do
        cat $OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep1_R2_val_2.fq.gz $OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep2_R2_val_2.fq.gz > $OUTPUT_DIR_TRIM/MergeFastq_${cell}_${cond}_R2.fq.gz
    done        
done

#Give permission to new files
for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do
        sudo chmod 777 $OUTPUT_DIR_TRIM/MergeFastq_${cell}_${cond}_R1.fq.gz
        sudo chmod 777 $OUTPUT_DIR_TRIM/MergeFastq_${cell}_${cond}_R2.fq.gz
    done    
done

echo "Merging of fastq replicate files done"

####################################
### DOVETAIL ALIGNMENT #############
####################################

for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do
        # Input files for alignment
        HICHIP_R1="$OUTPUT_DIR_TRIM/MergeFastq_${cell}_${cond}_R1.fq.gz"
        HICHIP_R2="$OUTPUT_DIR_TRIM/MergeFastq_${cell}_${cond}_R2.fq.gz"

        # Output files
        # BLF referse to black list filtered files
        MAPPED_PAIRS="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_mapped.pairs"
        MAPPED_BAM="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_mapped.PT.bam"
        MAPPED_BLF_BAM="$OUTPUT_HICHIP_ALIGN/BLF_Merge_${cell}_${cond}_hg38_nodd_mapped.PT.bam"

        # Alignment, dedup skipped. Can be included by removing comment mark.
        bwa mem -5SP -T0 -t$cores $REF_FASTA $HICHIP_R1 $HICHIP_R2 | \
        pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores2 --nproc-out $cores2 --chroms-path $REF_GENOME | \
        pairtools sort --tmpdir=$TEMP --nproc $cores | \
        #pairtools dedup --nproc-in $cores2 --nproc-out $cores2 --mark-dups --dry-run --output-stats Merge_${cell}_${cond}_stats.txt | \
        pairtools split --nproc-in $cores2 --nproc-out $cores2 --output-pairs $MAPPED_PAIRS --output-sam -|\
        samtools view -bS -@$cores | \
        samtools sort -@$cores -o $MAPPED_BAM;samtools index $MAPPED_BAM
        
        #Remove black listed regions
        bedtools intersect -v -abam $MAPPED_BAM -b $BLACKLIST > $MAPPED_BLF_BAM
        samtools index $MAPPED_BLF_BAM
    done  
done

echo "HiChIP Alignment Complete"

####################
### COVERAGE #######
####################

for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do
        # Input files for generating coverage files
        MAPPED_BAM="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_mapped.PT.bam"

        # Output file
        BIGWIG_OUT="$OUTPUT_HICHIP_SUB/BLF_Merge_${cell}_${cond}_hg38_nodd_mapped.bw"

        # Enrichment for IGV
        bamCoverage -b $MAPPED_BAM -o $BIGWIG_OUT --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

        echo "Generated Bigwig file $MAPPED_BAM "
    done
done

#########################
### CONTACT FILES #######
#########################
for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do

        # Input files
        MAPPED_PAIRS="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_mapped.pairs"
        # Output files
        CONTACT_MAP="$OUTPUT_HICHIP_SUB/Merge_${cell}_${cond}_hg38_nodd_contact_map.hic"

        # ContacMaps
        java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/HiChiP/juicer_tools_1.22.01.jar pre --threads "$cores" "$MAPPED_PAIRS" "$CONTACT_MAP" "$REF_GENOME"

        echo "Generated $CONTACT_MAP"
    done
done