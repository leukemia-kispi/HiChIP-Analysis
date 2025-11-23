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
PERMISSIVE_PVAL=1e-5
STRICT_PVAL=1e-9

########################################################
### GENERATE PRIMARY ALIGNMENT FILES AND MACS2 #########
########################################################
# ---------- merged ----------
echo
echo "Initiating conversion of HiChIP alignment files into primary alignments followed by MACS2 peak calling"

primary_aln_macs2(){
    local cell="$1"
    local cond="$2"

    #Set input files
    local MAPPED_MERGED="$INPUT_HICHIP_ALIGN/BLF_HiChIP_${cell}_${cond}_merged_nodd_mapped.PT.bam"
    #Set output files
    local PRIMARY_OUT="$INPUT_HICHIP_ALIGN/BLF_HiChIP_${cell}_${cond}_merged_nodd_primary.aln.bed"
    local PREFIX="BLF_HiChIP_${cell}_${cond}_merged_nodd_p${STRICT_PVAL}macs2"

    if [[ ! -f "$MAPPED_MERGED" ]]; then
        echo "Missing merged BAM: $MAPPED_MERGED — skipping ${cell}_${cond}"
        return 0
    fi

    echo "Generating primary aln: $PRIMARY_OUT"
    samtools view -b "$MAPPED_MERGED" -h -F 0x900 | bedtools bamtobed -i stdin > "$PRIMARY_OUT"

    if [[ ! -s "$PRIMARY_OUT" ]]; then
        echo "Primary file empty: $PRIMARY_OUT — skipping MACS2"
        return 0
    fi

    echo "Running strict MACS2 on $PRIMARY_OUT"
    macs2 callpeak -t "$PRIMARY_OUT" --keep-dup 10 --min-length 300 -p "$STRICT_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 --outdir "$OUTPUT_MACS2" -n "$PREFIX"
   
    echo "Peak calling finished: $PREFIX"

}

export -f primary_aln_macs2
export INPUT_HICHIP_ALIGN OUTPUT_MACS2 EFFECTIVE_GENOME STRICT_PVAL

if [[ ${#CellLine[@]} -ne 0 && ${#conditions[@]} -ne 0 ]]; then
    # Launch permissive runs in parallel
    parallel --ungroup --line-buffer -j "$JOBS" primary_aln_macs2 ::: "${CellLine[@]}" ::: "${conditions[@]}"
else
    echo "No merged samples found to process."
fi

# ---------- per-replicate ----------
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
        local PRIMARY_OUT="$INPUT_HICHIP_SUB/BLF_HiChIP_${cell}_${cond}_Rep${rep}_nodd_primary.aln.bed"
        local PREFIX="BLF_HiChIP_${cell}_${cond}_Rep${rep}_nodd_p${STRICT_PVAL}macs2"

        if [[ ! -f "$MAPPED_REP" ]]; then
            echo "Missing replicate BAM: $MAPPED_REP — skipping ${cell}_${cond}_Rep${rep}"
            return 0
        fi

        samtools view -b "$MAPPED_REP" -h -F 0x900 | bedtools bamtobed -i stdin > "$PRIMARY_OUT"
    
        if [[ ! -s "$PRIMARY_OUT" ]]; then
            echo "Primary file empty: $PRIMARY_OUT — skipping MACS2"
            return 0
        fi
        echo "Running strict MACS2 on $PRIMARY_OUT"
        macs2 callpeak -t "$PRIMARY_OUT" --keep-dup 10 --min-length 300 -p "$STRICT_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 --outdir "$OUTPUT_MACS2" -n "$PREFIX"

        echo "Peak calling finished: $PREFIX"
    }

    export -f primary_aln_macs2_rep
    export INPUT_HICHIP_SUB OUTPUT_MACS2 EFFECTIVE_GENOME STRICT_PVAL

    if [[ ${#CellLine[@]} -ne 0 && ${#conditions[@]} -ne 0 && ${#NUMBERS[@]} -ne 0 ]]; then
        # Launch permissive runs in parallel (per replicate)
        parallel --ungroup --line-buffer -j "$JOBS" primary_aln_macs2_rep ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"
    else
        echo "No replicate samlpes found."
    fi
else
    echo "Skip primary alignment for replicates"
fi

################################################
### OPTIONAL: PERMISSIVE MACS2 and IDR #########
################################################
echo
echo "IDR Analysis"
# Promt to procceed or skip IDR
read -rp "Do you want to proceed permissive MACS2 on replicates and perform IDR (y/n): " confirm_idr
if [[ "$confirm_idr" == "y" ]]; then

    primary_aln_IDR(){
        local cell="$1"
        local cond="$2"

        # Set input files
        local PRIMARY_OUT_MERGED="$INPUT_HICHIP_ALIGN/BLF_HiChIP_${cell}_${cond}_merged_nodd_primary.aln.bed"
        local MACS2_MERGED="$OUTPUT_MACS2/BLF_HiChIP_${cell}_${cond}_merged_nodd_p${STRICT_PVAL}macs2"
        local PRIMARY_OUT_REP1="$INPUT_HICHIP_SUB/BLF_HiChIP_${cell}_${cond}_Rep1_nodd_primary.aln.bed"
        local PRIMARY_OUT_REP2="$INPUT_HICHIP_SUB/BLF_HiChIP_${cell}_${cond}_Rep2_nodd_primary.aln.bed"

        # Set output files
        local PREFIX_PERMISSIVE_REP1="$OUTPUT_MACS2_PERMISSIVE/BLF_HiChIP_${cell}_${cond}_Rep1_nodd_p${PERMISSIVE_PVAL}macs2"
        local PREFIX_PERMISSIVE_REP2="$OUTPUT_MACS2_PERMISSIVE/BLF_HiChIP_${cell}_${cond}_Rep2_nodd_p${PERMISSIVE_PVAL}macs2"
        local MACS2_SORT_MERGED="$OUTPUT_MACS2_SORT/Sort_BLF_HiChIP_${cell}_${cond}_merged_nodd_p${STRICT_PVAL}macs2_peaks.narrowPeak"
        local MACS2_SORT_REP1="$OUTPUT_MACS2_SORT/Sort_BLF_HiChIP_${cell}_${cond}_Rep1_nodd_p${PERMISSIVE_PVAL}macs2_peaks.narrowPeak"
        local MACS2_SORT_REP2="$OUTPUT_MACS2_SORT/Sort_BLF_HiChIP_${cell}_${cond}_Rep2_nodd_p${PERMISSIVE_PVAL}macs2_peaks.narrowPeak"
        local IDR_OUT="$OUTPUT_MACS2_IDR/HiChIP_${cell}_${cond}_narrowPeak.idr"
        local IDR_LOG="$OUTPUT_MACS2_IDR/HiChIP_${cell}_${cond}_narrowPeak.idr.log"

        if [[ ! -f "$PRIMARY_OUT_REP1" || ! -f "$PRIMARY_OUT_REP2" ]]; then
            echo "Missing primary replicate bed(s) for ${cell}_${cond}, skipping IDR."
            return 0
        fi
       
        # run permissive MACS2 if outputs missing (macs2 writes into OUTPUT_MACS2_PERMISSIVE with basename prefix)
        if [[ ! -f "$OUTPUT_MACS2_PERMISSIVE/$(basename "$PREFIX_PERMISSIVE_REP1")_peaks.narrowPeak" ]]; then
            echo "Permissive MACS2 on rep1 for ${cell}_${cond}"
            macs2 callpeak -t "$PRIMARY_OUT_REP1" --keep-dup 10 --min-length 300 -p "$PERMISSIVE_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$(basename "$PREFIX_PERMISSIVE_REP1")"
        fi
        if [[ ! -f "$OUTPUT_MACS2_PERMISSIVE/$(basename "$PREFIX_PERMISSIVE_REP2")_peaks.narrowPeak" ]]; then
            echo "Permissive MACS2 on rep2 for ${cell}_${cond}"
            macs2 callpeak -t "$PRIMARY_OUT_REP2" --keep-dup 10 --min-length 300 -p "$PERMISSIVE_PVAL" -g "$EFFECTIVE_GENOME" --bw 300 --mfold 5 50 --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$(basename "$PREFIX_PERMISSIVE_REP2")"
        fi

        # Sort peak files
        sort -k8,8nr "${MACS2_MERGED}_peaks.narrowPeak" > "$MACS2_SORT_MERGED"
        sort -k8,8nr "${PREFIX_PERMISSIVE_REP1}_peaks.narrowPeak"  > "$MACS2_SORT_REP1"
        sort -k8,8nr "${PREFIX_PERMISSIVE_REP2}_peaks.narrowPeak"  > "$MACS2_SORT_REP2"
        
        if [[ ! -f "$MACS2_SORT_REP1" || ! -f "$MACS2_SORT_REP2" ]]; then
            echo "Missing sorted permissive peaks for IDR for ${cell}_${cond}; skipping IDR."
            return 0
        fi

        echo "Running IDR (rep1 vs rep2) for ${cell}_${cond}"
        idr --samples "$MACS2_SORT_REP1" "$MACS2_SORT_REP2" --input-file-type narrowPeak --rank p.value --output-file "$IDR_OUT" --plot --log-output-file "$IDR_LOG"
        
        if [[ -f "$MACS2_SORT_MERGED" ]]; then
        idr --samples "$MACS2_SORT_REP1" "$MACS2_SORT_REP2" --peak-list "$MACS2_SORT_MERGED" --input-file-type narrowPeak --rank p.value --output-file "${IDR_OUT%.idr}__withOracle.idr" --plot --log-output-file "${IDR_LOG%.idr.log}_withOracle.idr.log"

        echo "Completed IDR for ${cell}_${cond}"
    }

    export -f primary_aln_IDR
    export INPUT_HICHIP_ALIGN INPUT_HICHIP_SUB OUTPUT_MACS2 OUTPUT_MACS2_PERMISSIVE OUTPUT_MACS2_SORT OUTPUT_MACS2_IDR EFFECTIVE_GENOME PERMISSIVE_PVAL STRICT_PVAL
    
    if [[ ${#CellLine[@]} -ne 0 && ${#conditions[@]} -ne 0 ]]; then
        # Launch permissive runs in parallel (per replicate)
        parallel --ungroup --line-buffer -j "$JOBS" primary_aln_IDR ::: "${CellLine[@]}" ::: "${conditions[@]}"
    else
        echo "No samples found to run IDR."
    fi    
else
    echo "Skipping IDR"
fi

echo "Pipeline complete."