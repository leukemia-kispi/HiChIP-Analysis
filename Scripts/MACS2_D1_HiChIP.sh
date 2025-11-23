#!/usr/bin/env bash
set -euo pipefail # But fail properly on pipe errors
shopt -s nullglob # make globbing return empty array if no match

# Using this script assumes DirectoryArchitecture.sh was execute beforhand.
# This MACS2 pipeline expects inputs originating from pair-end read files.

# Check if GNU Parallel is installed otherwise exits
if ! command -v parallel &>/dev/null; then
    echo "Error: GNU parallel not found. Please install with 'sudo apt install parallel'"
    exit 1
fi

# Prompt the user for the main directory
read -rp "Enter the path to the MAIN_DIR: " MAIN_DIR
read -rp "How many cores does the machine have: " CORES
read -rp "How many threads do you assign per job (recommend 1): " TOOL_THREADS


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

echo "==== Starting MACS2 Pipeline ====="
echo "====================================="
echo "  MACS2 D1 Peak calling Log Started"
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

# Set input directories
INPUT_HICHIP_ALIGN="$MAIN_DIR/4.Alignment/HiChIP"
INPUT_HICHIP_SUB="$MAIN_DIR/4.Alignment/HiChIP/Outputs"
# Provide your own list if running different samples
MAPPING_FILE="$MAIN_DIR/1.RawData/HiChIP/HiChIP_IDs.txt" 
# Set output directories
OUTPUT_MACS2="$MAIN_DIR/5.MACS2/HiChIP"
OUTPUT_MACS2_SORT="$MAIN_DIR/5.MACS2/HiChIP/SORT"
OUTPUT_MACS2_PERMISSIVE="$MAIN_DIR/5.MACS2/HiChIP/Permissive"
OUTPUT_MACS2_IDR="$MAIN_DIR/5.MACS2/HiChIP/IDR"

# Ensure output directories exist
mkdir -p "$INPUT_HICHIP_ALIGN" "$INPUT_HICHIP_SUB" "$OUTPUT_MACS2" "$OUTPUT_MACS2_SORT" "$OUTPUT_MACS2_PERMISSIVE" "$OUTPUT_MACS2_IDR"

# Ensure there are no hidden spaces in the ID.txt file
if [[ -f "$MAPPING_FILE" ]]; then
    sed -i 's/\r$//' "$MAPPING_FILE"
else
    echo "Warning: mapping file $MAPPING_FILE not found. Continuing but arrays will be empty."
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

#################################################################
### INITIALIZE MACS2 CONDA ENVIRONMENT #########################
#################################################################

# Initialize Conda
eval "$(conda shell.bash hook)"

# Echo current Conda environment
echo "Current Conda environment: $CONDA_DEFAULT_ENV"

# Only activate if not already in the target environment
CONDA_ENV="MACS2"
if [[ "$CONDA_DEFAULT_ENV" != "$CONDA_ENV" ]]; then
    echo "Activating Conda environment: $CONDA_ENV"
    conda activate "$CONDA_ENV"
else
    echo "Already in the target Conda environment."
fi

################################
### MACS2 PEAK CALLING #########
################################

# Parameters can be adapted for personl needs
EFFECTIVE_GENOME=2913022398
PERMISSIVE_PVAL=0.00001
STRICT_PVAL=0.000000001

########################################################
### GENERATE PRIMARY ALIGNMENT FILES AND MACS2 #########
########################################################

echo
echo "Initiating conversion of HiChIP alignment files into primary alignments followed by MACS2 peak calling"

primary_aln_macs2(){
    local cell="$1"
    local cond="$2"

    #Set input files
    local MAPPED_MERGED="$INPUT_HICHIP_ALIGN/BLF_HiChIP_${cell}_${cond}_merged_nodd_mapped.PT.bam"
    #Set output files
    local PRIMARY_ALN_MERGED="$INPUT_HICHIP_ALIGN/BLF_HiChIP_${cell}_${cond}_merged_nodd_primary.aln.bed"
    local MACS2_MERGED="$OUTPUT_MACS2/BLF_HiChIP_${cell}_${cond}_merged_nodd_p9macs2"

    samtools view -b "$MAPPED_MERGED" -h -F 0x900 | bedtools bamtobed -i stdin > "$PRIMARY_ALN_MERGED"

    echo "Generated $PRIMARY_ALN_MERGED proceed with macs2 peak calling"

    macs2 callpeak -t "$PRIMARY_ALN_MERGED" --keep-dup 10 --min-length 300 -p "$STRICT_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 -n "$MACS2_MERGED"
   
    echo "Peak calling finished: $MACS2_MERGED"

}

export -f primary_aln_macs2
export INPUT_HICHIP_ALIGN OUTPUT_MACS2 EFFECTIVE_GENOME STRICT_PVAL

# Launch permissive runs in parallel (per replicate)
parallel --ungroup --line-buffer -j "$JOBS" primary_aln_macs2 ::: "${CellLine[@]}" ::: "${conditions[@]}"

echo
read -rp "Do you want to proceed with primary alignnments for individual replicates files (y/n): " confirm
if [[ "$confirm" == "y" ]]; then
    primary_aln_macs2_rep(){
        local cell="$1"
        local cond="$2"
        local rep="$3"
    
        #Set input files
        local MAPPED_REP="$INPUT_HICHIP_SUB/BLF_HiChIP_${cell}_${cond}_Rep${rep}_nodd_mapped.PT.bam"
        #Set output files
        local PRIMARY_ALN_REP="$INPUT_HICHIP_SUB/BLF_HiChIP_${cell}_${cond}_Rep${rep}_nodd_primary.aln.bed"
        local MACS2_REP="$OUTPUT_MACS2/BLF_HiChIP_${cell}_${cond}_Rep${rep}_nodd_p9macs2"

        samtools view -b "$MAPPED_REP" -h -F 0x900 | bedtools bamtobed -i stdin > "$PRIMARY_ALN_REP"

        echo "Generated $PRIMARY_ALN_REP proceed with macs2 peak calling"
    
        macs2 callpeak -t "$PRIMARY_ALN_REP" --keep-dup 10 --min-length 300 -p "$STRICT_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 -n "$MACS2_REP"

        echo "Peak calling finished: $MACS2_REP"
    }

    export -f primary_aln_macs2_rep
    export INPUT_HICHIP_SUB OUTPUT_MACS2 EFFECTIVE_GENOME STRICT_PVAL

    # Launch permissive runs in parallel (per replicate)
    parallel --ungroup --line-buffer -j "$JOBS" primary_aln_macs2_rep ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

else
    echo "Skip primary alignment for replicates"
fi

################################################
### OPTIONAL: PERMISSIVE MACS2 and IDR #########
################################################
echo
echo "IDR Analysis"
# Promt to procceed or skip IDR
read -rp "Do you want to proceed IDR Analysis (y/n): " confirm
if [[ "$confirm" == "y" ]]; then

    primary_aln_IDR(){
        local cell="$1"
        local cond="$2"


        # Set input files
        local PRIMARY_ALN_MERGED="$INPUT_HICHIP_ALIGN/BLF_HiChIP_${cell}_${cond}_merged_nodd_primary.aln.bed"
        local PREFIX_MERGED="BLF_HiChIP_${cell}_${cond}_merged_nodd_p9macs2"
        local PRIMARY_ALN_REP1="$INPUT_HICHIP_SUB/BLF_HiChIP_${cell}_${cond}_Rep1_nodd_primary.aln.bed"
        local PRIMARY_ALN_REP2="$INPUT_HICHIP_SUB/BLF_HiChIP_${cell}_${cond}_Rep2_nodd_primary.aln.bed"

        # Set output files
        local MACS2_PERMISSIVE_REP1="$OUTPUT_MACS2_PERMISSIVE/BLF_HiChIP_${cell}_${cond}_Rep1_nodd_p5macs2"
        local MACS2_PERMISSIVE_REP2="$OUTPUT_MACS2_PERMISSIVE/BLF_HiChIP_${cell}_${cond}_Rep2_nodd_p5macs2"
        local IDR_OUT="$OUTPUT_MACS2_IDR/${cell}_${cond}_cle_narrowPeak.idr"
        local IDR_LOG="$OUTPUT_MACS2_IDR/${cell}_${cond}_cle_narrowPeak.idr.log"

        # Perform alignment and skip if files already exist
        if [[ -f "$IDR_OUT" ]]; then
            echo "IDR already done for ${cell}_${cond}, skipping."
            return
        fi

        # PERMISSIVE MACS2 on replicate files
        local MACS2_PERMISSIVE_REP1="$OUTPUT_MACS2_PERMISSIVE/BLF_HiChIP_${cell}_${cond}_Rep1_nodd_p5macs2"
        local MACS2_PERMISSIVE_REP2="$OUTPUT_MACS2_PERMISSIVE/BLF_HiChIP_${cell}_${cond}_Rep2_nodd_p5macs2"
        macs2 callpeak -t "$PRIMARY_ALN_REP1" --keep-dup 10 --min-length 300 -p "$PERMISSIV_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 -n "$MACS2_PERMISSIVE_REP1"
        macs2 callpeak -t "$PRIMARY_ALN_REP2" --keep-dup 10 --min-length 300 -p "$PERMISSIV_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 -n "$MACS2_PERMISSIVE_REP2"

        # Sort peak files
        sort -k8,8nr "${MACS2_MERGED}_peaks.narrowPeak" > "$OUTPUT_MACS2_SORT/Sort_${PREFIX_MERGED}_peaks.narrowPeak"
        sort -k8,8nr "${MACS2_PERMISSIVE_REP1}_peaks.narrowPeak"  > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PERMISSIVE_REP1}_peaks.narrowPeak"
        sort -k8,8nr "${MACS2_PERMISSIVE_REP1}_peaks.narrowPeak"  > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PERMISSIVE_REP2}_peaks.narrowPeak"

        idr --samples $OUTPUT_MACS2_SORT/$MACS2_Rep1_SORT $OUTPUT_MACS2_SORT/$MACS2_Rep2_SORT --peak-list $OUTPUT_MACS2_SORT/$MACS2_JoinedRep_SORT --input-file-type narrowPeak --rank p.value --output-file "$IDR_OUT" --plot --log-output-file "$IDR_LOG"

        echo "IDR for HiChIP primary alignment HiChIP_${cell}_${cond} completed"
    }

    export -f primary_aln_IDR
    export INPUT_HICHIP_ALIGN INPUT_HICHIP_SUB OUTPUT_MACS2 OUTPUT_MACS2_PERMISSIVE EFFECTIVE_GENOME PERMISSIV_PVAL

    # Launch permissive runs in parallel (per replicate)
    parallel --ungroup --line-buffer -j "$JOBS" primary_aln_IDR ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

else
    echo "Skipping IDR"
fi

echo "Pipeline complete."