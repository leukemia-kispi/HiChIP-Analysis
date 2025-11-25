#!/usr/bin/env bash
set -euo pipefail # But fail properly on pipe errors
shopt -s nullglob # make globbing return empty array if no match

#Using this script assumes the script DirectoryArchitecture&CondaEnv.sh was run beforhand to create the DirectoryArchitecture and generate conda environments with needed tools.

# Check if GNU Parallel is installed otherwise exis
if ! command -v parallel &>/dev/null; then
    echo "Error: GNU parallel not found. Please install with 'sudo apt install parallel'"
    exit 1
fi

# Prompt the user for the main directory
read -rp "Enter the path to the MAIN_DIR: " MAIN_DIR
read -rp "How many cores does the machine have: " CORES
read -rp "How many threads do you assign per job (recommend 8,16,32): " TOOL_THREADS


# Check if input is empty
[[ -z "$MAIN_DIR" ]] && { echo "No directory entered"; exit 1; }
[[ -z "$CORES" ]] && { echo "Thread count missing"; exit 1; }
[[ -z "$TOOL_THREADS" ]] && { echo "Thread assigment per job missing"; exit 1; }

# Ensure TOOL_THREADS is integer and at least 1
TOOL_THREADS=$(( TOOL_THREADS ))
if [[ $TOOL_THREADS -lt 1 ]]; then
    TOOL_THREADS=1
fi

# Helper: safe half-threads (>=1)
half_threads() {
    local t=$(( TOOL_THREADS / 2 ))
    if [[ $t -lt 1 ]]; then t=1; fi
    echo "$t"
}

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

echo "==== Starting HiChIP-seq Pipeline ====="
echo "====================================="
echo "  HiChIP-seq Pipeline Log Started"
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
REF_GENOME="$MAIN_DIR/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="$MAIN_DIR/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for read before and after trimming, *fg.gz files
FASTQ_DIR="$MAIN_DIR/1.RawData/HiChIP"
# Provide your own list if running different samples
MAPPING_FILE="$MAIN_DIR/1.RawData/HiChIP/HiChIP_IDs.txt" 
# Set output directories
OUTPUT_DIR_TRIM="$MAIN_DIR/3.TRIM/HiChIP"
OUTPUT_HICHIP_ALIGN="$MAIN_DIR/4.Alignment/HiChIP"
OUTPUT_HICHIP_SUB="$MAIN_DIR/4.Alignment/HiChIP/Outputs"
BIGWIG_COVERAGE="$MAIN_DIR/7.Deeptool_Matrix/Coverage"
# Set Path to temporary directory
TEMP="$MAIN_DIR/tmp"

# Ensure output directories exist
mkdir -p "$OUTPUT_DIR_TRIM" "$OUTPUT_HICHIP_ALIGN" "$OUTPUT_HICHIP_SUB" "$BIGWIG_COVERAGE" "$TEMP"

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
        base=${newname#HiChIP_}      # Remove "HiChIP_" prefix
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
    local cell="$1"
    local cond="$2"
    local rep="$3"

    # Set path to input FASTQ files
    local R1="$FASTQ_DIR/HiChIP_${cell}_${cond}_Rep${rep}_R1.fastq.gz"
    local R2="$FASTQ_DIR/HiChIP_${cell}_${cond}_Rep${rep}_R2.fastq.gz"
    # Set path to output 
    local OUT_R1="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep${rep}_R1_val_1.fq.gz"
    local OUT_R2="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep${rep}_R2_val_2.fq.gz"

    # Perform trimming only if the expected output files are missing
    if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
        echo "Trimming already done for ${cell}_${cond}_Rep${rep}, skipping."
        return 0
    fi

    echo "Trimming ${cell}_${cond}_Rep${rep}..."
    # Trim samples and generate new FastQC files for all replicates
    trim_galore --fastqc --phred33 --length 50 --output_dir "$OUTPUT_DIR_TRIM" -j "$TOOL_THREADS" --paired "$R1" "$R2"
   
}

export -f trim_sample
export FASTQ_DIR OUTPUT_DIR_TRIM TOOL_THREADS

# Run in parallel
parallel -j "$JOBS" trim_sample ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

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
#Assuming only 2 replicates that are merged.

concatenate() {
    local cell="$1"
    local cond="$2"

# Input files
local READ1_REP1="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep1_R1_val_1.fq.gz"
local READ1_REP2="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep2_R1_val_1.fq.gz"
local READ2_REP1="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep1_R2_val_2.fq.gz"
local READ2_REP2="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_Rep2_R2_val_2.fq.gz"
# Output files 
local READ1_MERGED="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_merged_R1_val_1.fq.gz"
local READ2_MERGED="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_merged_R2_val_2.fq.gz"

# Concatenate R1 fastq files
cat "$READ1_REP1" "$READ1_REP2" > "$READ1_MERGED"

# Concatenate R2 fastq files
cat "$READ2_REP1" "$READ2_REP2" > "$READ2_MERGED"

#Give permission to new files
sudo chmod 777 "$READ1_MERGED"
sudo chmod 777 "$READ2_MERGED"

echo "Merging of fastq replicate files done ${cell}_${cond}_Rep${rep}..."

}

export -f concatenate
export  OUTPUT_DIR_TRIM

parallel -j "$JOBS" concatenate ::: "${CellLine[@]}" ::: "${conditions[@]}"

#######################
### ALIGNMENT #########
#######################

# Alignment funciton
align_sample() {
    local cell="$1"
    local cond="$2"

    # Set path to trimmed input FASTQ files
    local TRIM_R1="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_merged_R1_val_1.fq.gz"
    local TRIM_R2="$OUTPUT_DIR_TRIM/HiChIP_${cell}_${cond}_merged_R2_val_2.fq.gz"
    
    # Set path to aligned output files
    local MAPPED_PAIRS="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_mapped.pairs"
    local MAPPED_BAM="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_mapped.PT.bam"
    # Set output file name for filtered BAM files. BLF referse to black list filtered file
    local BLF_BAM="$OUTPUT_HICHIP_ALIGN/BLF_HiChIP_${cell}_${cond}_merged_nodd_mapped.PT.bam"
    
    #Perform alignment and skip if files already exist
    if [[ -f "$BLF_BAM" ]]; then
        echo "Alignment already done for ${cell}_${cond}, skipping."
        return 0
    fi
    
    echo "Aligning ${cell}_${cond}_..."

    # compute safe nproc values (>=1)
    local nproc_in
    nproc_in="$(half_threads)"

    # Alignment, dedup skipped. Can be included by removing comment mark.
    bwa mem -5SP -T0 -t"$TOOL_THREADS" "$REF_FASTA" "$TRIM_R1" "$TRIM_R2" | \
    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in "$nproc_in" --nproc-out "$nproc_in" --chroms-path "$REF_GENOME" | \
    pairtools sort --tmpdir="$TEMP" --nproc "$TOOL_THREADS" | \
    pairtools split --nproc-in "$nproc_in" --nproc-out "$nproc_in" --output-pairs $MAPPED_PAIRS --output-sam -|\
    samtools view -bS -@"$TOOL_THREADS" | \
    samtools sort -@"$TOOL_THREADS" -o "$MAPPED_BAM";samtools index "$MAPPED_BAM"

    # Blacklist filtering
    bedtools intersect -v -abam "$MAPPED_BAM" -b "$BLACKLIST" > "$BLF_BAM"
    samtools index "$BLF_BAM"
}

export -f align_sample
export OUTPUT_HICHIP_ALIGN REF_FASTA TEMP REF_GENOME BLACKLIST TOOL_THREADS OUTPUT_DIR_TRIM

parallel -j "$JOBS" align_sample ::: "${CellLine[@]}" ::: "${conditions[@]}"

####################
### FILTER PAIRS ###
####################

# The "paritools select" command in the DovtailHiChIP conda environment will filter
# and include unique and rescue mapped read pairs.

filter_pairs() {
    local cell="$1"
    local cond="$2"

    # Input files
    local MAPPED_PAIRS="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_mapped.pairs"

    # Output files
    local MAPPED_PAIRS_FILTERED="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_mapped.filtered.pairs"

    # Perform filtering if  output files are missing
    if [[ -f "$MAPPED_PAIRS_FILTERED" ]]; then
        echo "Filter pairs file already generated for ${cell}_${cond}, skipping."
        return 0
    fi

    pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' "$MAPPED_PAIRS" -o "$MAPPED_PAIRS_FILTERED"

    echo "Filterd paired files for ${cell}_${cond}."

}

export -f filter_pairs
export OUTPUT_HICHIP_ALIGN

parallel -j "$JOBS" filter_pairs ::: "${CellLine[@]}" ::: "${conditions[@]}"


################################
### HICPRO VALID PAIRS FILES ### 
################################

# The paired files generated during aligment needs to be converted
# into valid HiC-Pro files to proceed with FitHiChIP.

convert_pairs() {
    local cell="$1"
    local cond="$2"
   

    # Input files
    local MAPPED_PAIRS_FILTERED="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_mapped.filtered.pairs"

    # Output files
    local MAPPED_PAIRS_HICPRO="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_hicpro_mapped.filtered.pairs.gz"

    # Perform conversion if output files are missing
    if [[ -f "$MAPPED_PAIRS_HICPRO" ]]; then
        echo "Filter pairs file already generated for ${cell}_${cond}, skipping."
        return 0
    fi

    grep -v '#' "$MAPPED_PAIRS_FILTERED"| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > "$MAPPED_PAIRS_HICPRO"
            
    echo "Converted paired files for ${cell}_${cond} into HiC-Pro format"
}

export -f convert_pairs
export OUTPUT_HICHIP_ALIGN

parallel -j "$JOBS" convert_pairs ::: "${CellLine[@]}" ::: "${conditions[@]}"

####################
### COVERAGE #######
####################

coverage() {
    local cell="$1"
    local cond="$2"
    local rep="$3"

    # Input files for generating coverage files
    local MAPPED_BAM="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_mapped.PT.bam"

    # Output file
    local BIGWIG_OUT="$OUTPUT_HICHIP_SUB/HiChIP_${cell}_${cond}_merged_nodd_mapped_RPKM.bw"
    
    # Perform coverage analysis if output files are missing
    if [[ -f "$BIGWIG_OUT" ]]; then
        echo "Coverage already generated for ${cell}_${cond}, skipping."
        return 0
    fi

    # Skip if missing input files
    if [[ ! -f "$MAPPED_BAM" ]]; then
        echo "Input bam file not found — skipping coverage for ${cell}_${cond}."
        return 0
    fi

    # Enrichment for IGV
    bamCoverage -b "$MAPPED_BAM" -o "$BIGWIG_OUT" --effectiveGenomeSize 2913022398 -bl "$BLACKLIST" --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

    echo "Generated BigWig file for ${cell}_${cond}"

}

export -f coverage
export OUTPUT_HICHIP_ALIGN OUTPUT_HICHIP_SUB BLACKLIST

parallel -j "$JOBS" coverage ::: "${CellLine[@]}" ::: "${conditions[@]}"

#########################
### CONTACT FILES #######
#########################

contact() {
    local cell="$1"
    local cond="$2"


    # Input files
    local MAPPED_PAIRS="$OUTPUT_HICHIP_ALIGN/HiChIP_${cell}_${cond}_merged_nodd_mapped.pairs"
    # Output files
    local CONTACT_MAP="$OUTPUT_HICHIP_SUB/HiChIP_${cell}_${cond}_merged_nodd_contact_map.hic"
    
    # Perform generation of contac maps if output files are missing
    if [[ -f "$CONTACT_MAPT" ]]; then
        echo "ontact map already exists for ${cell}_${cond}, skipping."
        return 0
    fi
    # Skip if missing input file
    if [[ ! -f "$MAPPED_PAIRS" ]]; then
        echo "Mapped pairs input missing: ${cell}_${cond} — skipping contact map."
        return 0
    fi

    echo "Generating contact map for ${cell}_${cond}_Rep${rep}..."
    # ContacMaps
    java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/HiChiP/juicer_tools_1.22.01.jar pre --threads "$TOOL_THREADS" "$MAPPED_PAIRS" "$CONTACT_MAP" "$REF_GENOME"

echo "Generated contact maps for ${cell}_${cond}."

}

export -f contact
export OUTPUT_HICHIP_ALIGN OUTPUT_HICHIP_SUB TOOLS_THREADS REF_GENOME

parallel -j "$JOBS" contact ::: "${CellLine[@]}" ::: "${conditions[@]}"

echo "HiChIP aligmnet for replicates, coverage file and contact map file generation complete"
