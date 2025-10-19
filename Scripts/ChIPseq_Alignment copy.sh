#!/usr/bin/env bash
set +e                 # Don’t exit on error — handle manually
set -o pipefail        # But fail properly on pipe errors
shopt -s nullglob # make globbing return empty array if no match

# Using this script assumes DirectoryArchitecture.sh was execute beforhand.
# Expected to provide pair-end read ChIP-seq fastq files.

# Check if GNU Parallel is installed otherwise exis
if ! command -v parallel &>/dev/null; then
    echo "Error: GNU parallel not found. Please install with 'sudo apt install parallel'"
    exit 1
fi

# Prompt the user for the main directory
read -rp "Enter the path to the MAIN_DIR: " MAIN_DIR
read -rp "How many cores does the machine have: " CORES
read -rp "How many threads do you assign per job (recommend 4): " TOOL_THREADS


# Check if input is empty
[[ -z "$MAIN_DIR" ]] && { echo "No directory entered"; exit 1; }
[[ -z "$CORES" ]] && { echo "Thread count missing"; exit 1; }
[[ -z "$TOOL_THREADS" ]] && { echo "Thread assigment per job missing"; exit 1; }


# Compute number of parallel jobs (at least 1)
JOBS=$(( CORES / TOOL_THREADS ))
if [[ $JOBS -lt 1 ]]; then JOBS=1; fi

####################################
### SETUP LOGGING ##################
####################################
LOG_DIR="$MAIN_DIR/logs"
mkdir -p "$LOG_DIR"

# Capture all stdout and stderr into run.log while still showing it on screen
exec > >(tee -a "$LOG_DIR/run.log") 2>&1

echo "==== Starting ChIP-seq Pipeline ====="
echo "====================================="
echo "  ChIP-seq Pipeline Log Started"
echo "  Date: $(date)"
echo "  Main Directory: $MAIN_DIR"
echo "  Using $CORES cores"
echo "  Threads per job: $TOOL_THREADS"
echo "  Parallel jobs (JOBS): $JOBS"
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

# Ensure there are no hidden spaces in the ID.txt file
if [[ -f "$MAPPING_FILE" ]]; then
    sed -i 's/\r$//' "$MAPPING_FILE"
else
    echo "Warning: mapping file $MAPPING_FILE not found. Continuing but arrays will be empty."
fi

#############################
### MAPPING NEW IDs #########
#############################

# Paired_Read = read_num
# True sample name = newname
# MAPPING_FILES is a .txt file with mapping_file format (acc newname):
# Example from European Nucleotide Archive deposited under accession number ERP109232. Use the Histone_ChIP_IDs.txt file to update the annotations.
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
if [[ -f "$MAPPING_FILE" ]]; then
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
fi

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

# Trimming funciton
trim_sample() {
    local cell=$1
    local cond=$2
    local num=$3
     # Set path to input FASTQ files
    R1="$FASTQ_DIR/ChIP_${cell}_${cond}_Rep${num}_R1.fastq.gz"
    R2="$FASTQ_DIR/ChIP_${cell}_${cond}_Rep${num}_R2.fastq.gz"
     # Set path to output 
    OUT_R1="$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R1_val_1.fq.gz"
    OUT_R2="$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R2_val_2.fq.gz"

    # Perform trimming only if the expected output files are missing
    if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
        echo "Trimming already done for ${cell}_${cond}_Rep${num}, skipping."
        return
    fi

    echo "Trimming ${cell}_${cond}_Rep${num}..."
    # Trim samples and generate new FastQC files for all replicates
    trim_galore --fastqc --phred33 --length 30 --output_dir "$OUTPUT_DIR_TRIM" -j "$TOOL_THREADS" --paired "$R1" "$R2"
   
}

export -f trim_sample
export FASTQ_DIR OUTPUT_DIR_TRIM TOOL_THREADS

# Run in parallel
parallel -j "$JOBS" trim_sample ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

#######################
### ALIGNMENT #########
#######################

# Alignment funciton
align_sample() {
    local cell=$1
    local cond=$2
    local num=$3
    # Set path to trimmed input FASTQ files
    TRIM_R1="$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R1_val_1.fq.gz"
    TRIM_R2="$OUTPUT_DIR_TRIM/ChIP_${cell}_${cond}_Rep${num}_R2_val_2.fq.gz"
    # Set path to aligned output files
    BAM="$OUTPUT_CHIP_ALIGN/ChIP_${cell}_${cond}_Rep${num}_cle_sort.bam"
    # Set output file name for filtered BAM files. BLF referse to black list filtered file
    BLF_BAM="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort.bam"
    
    #Perform alignment and skip if files already exist
    if [[ -f "$BLF_BAM" ]]; then
        echo "Alignment already done for ${cell}_${cond}_Rep${num}, skipping."
        return
    fi

    echo "Aligning ${cell}_${cond}_Rep${num}..."
    #Create sorted BAM files with grep to remove alignments to alternative contigs, unlocalized sequence, or unplaced sequence.#####################
    bwa mem -5 -T25 -t"$TOOL_THREADS" "$REF_FASTA" "$TRIM_R1" "$TRIM_R2" \
        | samtools view -hS \
        | grep -v chrUn | grep -v random | grep -v _alt \
        | samtools view -bS -@"$TOOL_THREADS"  \
        | samtools sort -@"$TOOL_THREADS"  -o "$BAM"

    # Blacklist filtering
    bedtools intersect -v -abam "$BAM" -b "$BLACKLIST" > "$BLF_BAM"
    samtools index "$BLF_BAM"
}

export -f align_sample
export OUTPUT_CHIP_ALIGN REF_FASTA BLACKLIST TOOL_THREADS OUTPUT_DIR_TRIM

parallel -j "$JOBS" align_sample ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

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

# Picard funciton
mark_duplicates() {
    local cell=$1
    local cond=$2
    local num=$3
    # Set path to aligned input files
    BLF_BAM="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort.bam"
    # Set path to flagged output files
    BLF_DUPFLAG_BAM="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dupsflag.bam"
    BLF_DUPFLAG_TXT="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dups.txt"
    
    # Flag duplicates and skip if files exist
    if [[ -f "$BLF_DUPFLAG_BAM" ]]; then
        echo "Duplicate marking already done for ${cell}_${cond}_Rep${num}, skipping."
        return
    fi

    java -jar "$(conda run -n Picard bash -c 'echo $CONDA_PREFIX')/share/$(ls "$(conda run -n Picard bash -c 'echo $CONDA_PREFIX')/share" | grep picard | sort -V | tail -n 1)/picard.jar" \
    MarkDuplicates -I "$BLF_BAM" -O "$BLF_DUPFLAG_BAM" -M "$BLF_DUPFLAG_TXT" --REMOVE_DUPLICATES false
}

export -f mark_duplicates
export OUTPUT_CHIP_ALIGN OUTPUT_CHIP_SUB

parallel -j "$JOBS" mark_duplicates ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

#################################################################
### RE-ACTIVATE DOVETAILHICHIP CONDA ENVIRONMENT ################
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

dedup_sample() {
    local cell=$1
    local cond=$2
    local num=$3

    # Input: BAM file with duplicates flagged by Picard
    BLF_DUPFLAG_BAM="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dupsflag.bam"
    # Output: deduplicated BAM
    BLF_DD_BAM="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dd.bam"

    # Skip if deduplication already done
    if [[ -f "$BLF_DD_BAM" ]]; then
        echo "Deduplication already done for ${cell}_${cond}_Rep${num}, skipping."
        return
    fi

    if [[ ! -f "$BLF_DUPFLAG_BAM" ]]; then
        echo "Input flagged BAM not found for ${cell}_${cond}_Rep${num}, skipping dedup."
        return
    fi

    echo "Deduplicating ${cell}_${cond}_Rep${num}..."
    # Remove flagged duplicates (F 1024) and index
    samtools view -b -F 1024 "$BLF_DUPFLAG_BAM" > "$BLF_DD_BAM"
    samtools index "$BLF_DD_BAM"
}

export -f dedup_sample
export OUTPUT_CHIP_SUB

# Run in parallel
parallel -j "$JOBS" dedup_sample ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

############################################
### MERGE REPLICATES #######################
############################################

merge_replicates() {
    local cell=$1
    local cond=$2
    # Set input replicate files for merging
    REP1="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep1_cle_sort_dd.bam"
    REP2="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep2_cle_sort_dd.bam"
    # Set output merged file
    MERGED="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd.bam"

    # Skip if merged files already exist
    if [[ -f "$MERGED" ]]; then
        echo "Merged BAM already exists for ${cell}_${cond}, skipping."
        return
    fi

    # Merge replicate files if both exist
    if [[ -f "$REP1" && -f "$REP2" ]]; then
        echo "Merging replicates for ${cell}_${cond}..."
        samtools merge -f "$MERGED" "$REP1" "$REP2"
        samtools index "$MERGED"
    else
        echo "Cannot merge ${cell}_${cond}: missing replicate(s)"
        [[ ! -f "$REP1" ]] && echo "  - Missing: $(basename "$REP1")"
        [[ ! -f "$REP2" ]] && echo "  - Missing: $(basename "$REP2")"
    fi
}

export -f merge_replicates
parallel -j "$JOBS" merge_replicates ::: "${CellLine[@]}" ::: "${conditions[@]}"

##############################
### GENOMIC COVERAGE #########
##############################

coverage_sample() {
    local cell=$1
    local cond=$2
    local num=$3

    # Input deduplicated BAM
    BLF_DD_BAM="$OUTPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dd.bam"
    # Output BigWig coverage
    BLF_BIGWIG="$BIGWIG_COVERAGE/BLF_ChIP_${cell}_${cond}_Rep${num}_cle_sort_dd.bw"

    # Skip if coverage file already exists
    if [[ -f "$BLF_BIGWIG" ]]; then
        echo "Coverage already generated for ${cell}_${cond}_Rep${num}, skipping."
        return
    fi

    if [[ ! -f "$BLF_DD_BAM" ]]; then
        echo "Dedup BAM not found for ${cell}_${cond}_Rep${num}, skipping coverage."
        return
    fi

    echo "Generating coverage for ${cell}_${cond}_Rep${num}..."
    bamCoverage -b "$BLF_DD_BAM" -o "$BLF_BIGWIG" \
        --effectiveGenomeSize 2913022398 -bl "$BLACKLIST" \
        --normalizeUsing RPKM -p "$TOOL_THREADS" -bs 10 \
        --extendReads --ignoreForNormalization M
}

coverage_merged() {
    local cell=$1
    local cond=$2

    BLF_MERGED_BAM="$OUTPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd.bam"
    BLF_BIGWIG_MERGED="$BIGWIG_COVERAGE/BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd.bw"

    # Skip if merged BAM or output BigWig does not exist
    if [[ ! -f "$BLF_MERGED_BAM" ]]; then
        echo "Merged BAM not found for ${cell}_${cond}, skipping coverage."
        return
    fi
    if [[ -f "$BLF_BIGWIG_MERGED" ]]; then
        echo "Coverage already generated for merged ${cell}_${cond}, skipping."
        return
    fi

    echo "Generating coverage for merged ${cell}_${cond}..."
    bamCoverage -b "$BLF_MERGED_BAM" -o "$BLF_BIGWIG_MERGED" \
        --effectiveGenomeSize 2913022398 -bl "$BLACKLIST" \
        --normalizeUsing RPKM -p "$TOOL_THREADS" -bs 10 \
        --extendReads --ignoreForNormalization M
}

export -f coverage_sample coverage_merged
export OUTPUT_CHIP_SUB OUTPUT_CHIP_ALIGN BIGWIG_COVERAGE BLACKLIST TOOL_THREADS

# Replicate coverage
parallel -j "$JOBS" coverage_sample ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

# Merged coverage
parallel -j "$JOBS" coverage_merged ::: "${CellLine[@]}" ::: "${conditions[@]}"

echo "ChIP-seq Aligment, Clean-up and Coverage file generation completed"