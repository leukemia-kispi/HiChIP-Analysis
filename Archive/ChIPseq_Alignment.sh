#!/usr/bin/env bash
set -e #Exit immediately if any command exits with a non-zero status
set -u #Treat unset variables as an error when substituting.
shopt -s nullglob # make globbing return empty array if no match

# Using this script assumes DirectoryArchitecture.sh was execute beforhand

# Prompt the user for the main directory
read -rp "Enter the path to the MAIN_DIR: " MAIN_DIR
read -rp "How many cores do you assign: " CORES

# Check if input is empty
if [ -z "$MAIN_DIR" ]; then
    echo "Error: No directory path entered."
    exit 1
fi

#Check if input for CPU cores is empty
if [ -z "$CORES" ]; then
    echo "Error: Decide thread assigment."
    exit 1
fi

####################################
### SETUP LOGGING ##################
####################################
LOG_DIR="$MAIN_DIR/logs"
mkdir -p "$LOG_DIR"

# Capture all stdout and stderr into run.log while still showing it on screen
exec > >(tee -a "$LOG_DIR/run.log") 2>&1

echo "====================================="
echo "  ChIP-seq Pipeline Log Started"
echo "  Date: $(date)"
echo "  Main Directory: $MAIN_DIR"
echo "  Using $CORES cores"
echo "  Log file: $LOG_DIR/run.log"
echo "====================================="

#############################
### MAIN FILE PATHS #########
#############################

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="$MAIN_DIR/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
BLACKLIST="$MAIN_DIR/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for read files, *fastq.gz and ID mapping file
FASTQ_DIR="$MAIN_DIR/1.RawData/ChIP" #Enusure read files are uploaded to this directory
MAPPING_FILE="$MAIN_DIR/1.RawData/ChIP/ChIP_IDs.txt" #Provide your own list if running different samples
# Set output directories
OUTPUT_DIR_TRIM="$MAIN_DIR/3.TRIM/ChIP"
OUTPUT_CHIP_ALIGN="$MAIN_DIR/4.Alignment/ChIP"
OUTPUT_CHIP_SUB="$MAIN_DIR/4.Alignment/ChIP/Outputs"
BIGWIG_COVERAGE="$MAIN_DIR/7.Deeptool_Matrix/Coverage"

#Ensure there are no hidden spaces in the ID.txt file
sed -i 's/\r$//' "$MAPPING_FILE"

#############################
### MAPPING NEW IDs #########
#############################

#Assuming the fastqc files were for HAL-01 ChIP-seq for histone marks, downloaded from the European Nucleotide Archive deposited under accession number ERP109232. Use the Histone_ChIP_IDs.txt file to update the annotations.
# Accession IDs = acc
# Paired_Read = read_num
# True sample name = newname
# MAPPING_FILES is a .txt file with mapping_file format (acc newname):
#     ERR2618839  ChIP_HAL01_H3K27ac_Rep1
#     ERR2618840  ChIP_HAL01_H3K27ac_Rep2

echo
echo "=== DRY RUN: Planned renames ==="
# Promt to procceed or skip sample ID mapping step.
read -rp "Do you need to map new IDs using new_ID .txt file (y/n): " confirm
if [[ "$confirm" == "y" ]]; then
    while IFS=$'\t' read -r acc newname; do
        # Skip empty lines or lines starting with #
        [[ -z "$acc" || "$acc" == \#* ]] && continue

        for read_num in 1 2; do
            #Variable for old and new file names
            old_file="$FASTQ_DIR/${acc}_${read_num}.fastq.gz"
            new_file="$FASTQ_DIR/${newname}_R${read_num}.fastq.gz"
            # Check if old file exists and if renamed file exists does not overwrite.
            if [[ -f "$old_file" ]]; then
                if [[ -f "$new_file" ]]; then
                    echo "SKIP: $(basename "$new_file") already exists — will not overwrite"
                else
                    echo "$(basename "$old_file") → $(basename "$new_file")"
                fi
            else
                echo "WARNING: $(basename "$old_file") not found"
            fi
        done
    done < "$MAPPING_FILE"

    #Ask user if they wish to proceed with file name changes.
    echo
    read -rp "Do you want to apply these changes? (y/n): " confirm
    if [[ "$confirm" == "y" ]]; then
        while IFS=$'\t' read -r acc newname; do
            #Ignore blank lines and comment lines in the mapping file.
            [[ -z "$acc" || "$acc" == \#* ]] && continue

            for read_num in 1 2; do
                old_file="$FASTQ_DIR/${acc}_${read_num}.fastq.gz"
                new_file="$FASTQ_DIR/${newname}_R${read_num}.fastq.gz"

                if [[ -f "$old_file" && ! -f "$new_file" ]]; then
                    mv "$old_file" "$new_file"
                    echo "Renamed: $(basename "$old_file") → $(basename "$new_file")"
                elif [[ -f "$new_file" ]]; then
                    echo "SKIP: $(basename "$new_file") already exists — not overwritten"
                fi
            done
        done < "$MAPPING_FILE"
    else
        echo "No changes applied."
    fi
else
 echo "Skip mapping new IDs"
fi

################################################
### EXTRACT ANNOTATIONS FROM SAMPLE NAME #######
################################################

# Initialize arrays
CellLine=()
conditions=()
NUMBERS=()

# Read mapping file and extract parts of newname
while IFS=$'\t' read -r acc newname; do
    [[ -z "$acc" || "$acc" == \#* ]] && continue

    base=${newname#ChIP_}      # Remove "ChIP_" prefix
    cell=$(echo "$base" | cut -d'_' -f1)
    cond=$(echo "$base" | cut -d'_' -f2)
    rep=$(echo "$base" | grep -oP 'Rep\K[0-9]+')
    
    # Append unique values into arrays
    [[ " ${CellLine[*]} " != *" $cell "* ]] && CellLine+=("$cell")
    [[ " ${conditions[*]} " != *" $cond "* ]] && conditions+=("$cond")
    [[ " ${NUMBERS[*]} " != *" $rep "* ]] && NUMBERS+=("$rep")
done < "$MAPPING_FILE"

echo "CellLine: ${CellLine[@]}"
echo "conditions: ${conditions[@]}"
echo "NUMBERS: ${NUMBERS[@]}"

##########################################################
### INITIALIZE DOVETAILHICHIP CONDA ENVIRONMENT #########
##########################################################

#Environment should have following installed: bwa-mem, samtools, pairtools, fastqc, bedtools, multiqc, trim-galore, deeptools

# Initialize Conda reliably
source ~/anaconda3/etc/profile.d/conda.sh

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
            if [ ! -f "$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R1_val_1.fq.gz" ] || [ ! -f "$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R2_val_2.fq.gz" ]; then
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
                READ1="$FASTQ_DIR/ChIP_${cell}_${cond}_Rep${num}_R1.fastq.gz"
                READ2="$FASTQ_DIR/ChIP_${cell}_${cond}_Rep${num}_R2.fastq.gz"
        
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

#######################
### ALIGNMENT #########
#######################

#Loop through each pair of trimmed FASTQ files to complete alignment
for cell in "${CellLine[@]}"; do     
    for cond in "${conditions[@]}"; do
        for num in "${NUMBERS[@]}"; do
            # Set path to trimmed FASTQ files usign wildcard pattern
            TRIMMED_R1_FILES="$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R1_val_1.fq.gz"
            TRIMMED_R2_FILES="$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R2_val_2.fq.gz"
            # Set output file name for BAM files.
            MAPPED_BAM="$OUTPUT_CHIP_ALIGN/ChIP_${cell}_${cond}_Rep${num}_cle_sort.bam"
    
            #Create sorted BAM files with grep to remove alignments to alternative contigs, unlocalized sequence, or unplaced sequence.#####################
            bwa mem -5 -T25 -t"$CORES" "$REF_FASTA" "$TRIMMED_R1_FILES" "$TRIMMED_R2_FILES" | samtools view -hS | grep -v chrUn | grep -v random | grep -v _alt | samtools view -bS -@"$CORES" | samtools sort -@"$CORES" -o "$MAPPED_BAM" 

            # Set output file name for filtered BAM files. BLF referse to black list filtered file
            MAPPED_BLF_BAM="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort.bam"

            #Remove BlackList region & index
            bedtools intersect -v -abam "$MAPPED_BAM" -b "$BLACKLIST" > "$MAPPED_BLF_BAM"
            samtools index "$MAPPED_BLF_BAM"
        
            echo "Alignment, blacklist removal and indexing done for sample $(basename "$MAPPED_BLF_BAM") "
        done
    done
done

##########################################################
### INITIALIZE PICARD CONDA ENVIRONMENT #################
##########################################################

# Echo current Conda environment
echo "Current Conda environment: $CONDA_DEFAULT_ENV"

# Only activate if not already in the target environment
CONDA_ENV="Picard"
if [[ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV" ]]; then
    echo "Activating Conda environment: $CONDA_ENV"
    conda activate "$CONDA_ENV"
else
    echo "Already in the target Conda environment."
fi

######################################
### DUPLICATION FLAG: PICARD #########
######################################
for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do
        for num in "${NUMBERS[@]}"; do
            # Set input file names
            MAPPED_BLF_BAM="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort.bam"

            # Set output file name for dedup files.
            MAPPED_BLF_BAM_DUPFLAG="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dupsflag.bam"
            MAPPED_BLF_TXT_DUP="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dups.txt"
    
            #Picard MarkDuplicates with universal pathing for the picard.jar
            java -jar "$(conda run -n Picard bash -c 'echo $CONDA_PREFIX')/share/$(ls "$(conda run -n Picard bash -c 'echo $CONDA_PREFIX')/share" | grep picard | sort -V | tail -n 1)/picard.jar" MarkDuplicates -I "$MAPPED_BLF_BAM" -O "$MAPPED_BLF_BAM_DUPFLAG" -M "$MAPPED_BLF_TXT_DUP" --REMOVE_DUPLICATES false
        
            echo "Flagging duplicates in $(basename "$MAPPED_BLF_BAM_DUPFLAG")  done"
        done
    done
done

#################################################################
### INITIALIZE DOVETAILHICHIP CONDA ENVIRONMENT ################
#################################################################

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

##############################################
### DEDUPLICATION: CHIP SAMPLES ##############
##############################################
for cell in "${CellLine[@]}"; do      
    for cond in "${conditions[@]}"; do
        for num in "${NUMBERS[@]}"; do
            # Set input file names
            MAPPED_BLF_BAM_DUPFLAG="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dupsflag.bam"
            # Set output file name for dedup files.
            MAPPED_BLF_BAM_DD="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dd.bam"

            # Generate final deduplicated BAM files (remove flagged duplicates) and indexing
            samtools view -b -F 1024 "$MAPPED_BLF_BAM_DUPFLAG" > "$MAPPED_BLF_BAM_DD"
            samtools index "$MAPPED_BLF_BAM_DD"
            echo "Duplicates removed in $(basename "$MAPPED_BLF_BAM_DD") done"
        done    
    done
done

# Merge replicates only if both exist
# Set output file name for dedup merged files.
for cell in "${CellLine[@]}"; do   
    for cond in "${conditions[@]}"; do
        # Set file name for black list filterd, cleaned, sorted and deduplicated replicate BAM files
        MAPPED_BLF_BAM_DD_REP1="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep1_cle_sort_dd.bam" 
        MAPPED_BLF_BAM_DD_REP2="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep2_cle_sort_dd.bam" 
        # Set output file name for merged  BAM files
        MAPPED_MERGED_BLF_BAM_DD="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd.bam"

        if [[ -f "$MAPPED_BLF_BAM_DD_REP1" && -f "$MAPPED_BLF_BAM_DD_REP2" ]]; then
            echo "Merging replicates for ${cell}_${cond}..."    
            # MERGE Replicate File and index
            samtools merge -f "$MAPPED_MERGED_BLF_BAM_DD" "$MAPPED_BLF_BAM_DD_REP1" "$MAPPED_BLF_BAM_DD_REP2"
            samtools index "$OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD"
            echo "Merged file created: $(basename "$MAPPED_MERGED_BLF_BAM_DD")"
        else
            echo "Skipping merge for ${cell}_${cond}: missing replicate(s)"
            [[ ! -f "$MAPPED_BLF_BAM_DD_REP1" ]] && echo "  - Missing: $(basename "$MAPPED_BLF_BAM_DD_REP1")"
            [[ ! -f "$MAPPED_BLF_BAM_DD_REP2" ]] && echo "  - Missing: $(basename "$MAPPED_BLF_BAM_DD_REP2")"
        fi
    done
done

##############################
### GENOMIC COVERAGE #########
##############################
for cell in "${CellLine[@]}"; do   
    for cond in "${conditions[@]}"; do
        for num in "${NUMBERS[@]}"; do
            # Set input file name
            MAPPED_BLF_BAM_DD="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dd.bam"
            # Set output file name for bigwig file
            BIGWIG_BLF_DD="$BIGWIG_COVERAGE/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dd.bw"

            # Generate Enrichment coverage files for IGV visualization for each replicate
            bamCoverage -b "$MAPPED_BLF_BAM_DD" -o "$BIGWIG_BLF_DD" --effectiveGenomeSize 2913022398 -bl "$BLACKLIST" --normalizeUsing RPKM -p max -bs 10 \
            --extendReads --ignoreForNormalization M
        done    
    done
done

# Generate Enrichment coverage files for IGV visualization for merged files if they exist.
for cell in "${CellLine[@]}"; do   
    for cond in "${conditions[@]}"; do
        # Set input file name
        MAPPED_MERGED_BLF_BAM_DD="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd.bam"
        #Set output file name for bigwig file
        BIGWIG_MERGED_BLF_DD="$BIGWIG_COVERAGE/BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd.bw"

        if [[ -f "$MAPPED_MERGED_BLF_BAM_DD" ]]; then
            echo "Generating coverage for merged sample ${cell}_${cond}..."
            bamCoverage -b  "$MAPPED_MERGED_BLF_BAM_DD" -o "$BIGWIG_MERGED_BLF_DD" --effectiveGenomeSize 2913022398 -bl "$BLACKLIST" --normalizeUsing RPKM -p max -bs 10 \
            --extendReads --ignoreForNormalization M
            echo "Coverage file created: $(basename "$BIGWIG_MERGED_BLF_DD")"
        else
            echo "Skipping merged coverage: merged BAM not found for ${cell}_${cond}"
        fi
    done
done

echo "ChIP-seq Aligment, Clean-up and Coverage file generation completed"