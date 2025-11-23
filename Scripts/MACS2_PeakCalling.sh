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
echo "  MACS2 Peak calling Log Started"
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
INPUT_CHIP_ALIGN="$MAIN_DIR/4.Alignment/ChIP"
INPUT_CHIP_SUB="$MAIN_DIR/4.Alignment/ChIP/Outputs"
# Provide your own list if running different samples
MAPPING_FILE="$MAIN_DIR/1.RawData/ChIP/ChIP_IDs.txt" 
# Set output directories
OUTPUT_MACS2="$MAIN_DIR/5.MACS2/ChIP"
OUTPUT_MACS2_SORT="$MAIN_DIR/5.MACS2/ChIP/SORT"
OUTPUT_MACS2_PERMISSIVE="$MAIN_DIR/5.MACS2/ChIP/Permissive"
OUTPUT_MACS2_IDR="$MAIN_DIR/5.MACS2/ChIP/IDR"

# Ensure output directories exist
mkdir -p "$INPUT_CHIP_ALIGN" "$INPUT_CHIP_SUB" "$OUTPUT_MACS2" "$OUTPUT_MACS2_SORT" "$OUTPUT_MACS2_PERMISSIVE" "$OUTPUT_MACS2_IDR"

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

#################################################################
### INITIALIZE MACS2 CONDA ENVIRONMENT #########################
#################################################################

# Initialize Conda reliably
source ~/anaconda3/etc/profile.d/conda.sh

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
PERMISSIVE_PVAL=0.05
STRICT_PVAL=1e-9

###################################
# MACS2 permissive per replicate  #
###################################

echo
echo "Permissive MACS2 for later IDR Analysis"
read -rp "Do you want to proceed with Permissive MACS2 for later IDR Analysis (y/n): " confirm
if [[ "$confirm" == "y" ]]; then

    macs2_permissive_rep() {
        local cell=$1
        local cond=$2
        local rep=$3

        # Set inpute files
        local IN_BAM="$INPUT_CHIP_SUB/BLF_ChIP_${cell}_${cond}_Rep${rep}_cle_sort_dd.bam"
        local CTR_BAM="$INPUT_CHIP_SUB/BLF_ChIP_${cell}_controlInput_Rep${rep}_cle_sort_dd.bam"
        local PREFIX="BLF_ChIP_${cell}_${cond}_Rep${rep}_cle_sort_dd_p${PERMISSIVE_PVAL}macs2"

        # Check inputs
        if [[ ! -f "$IN_BAM" ]]; then
            echo "Missing input BAM: $IN_BAM — skipping ${cell}_${cond}_Rep${rep}"
            return 0
        fi

        # Check control
        if [[ ! -f "$CTR_BAM" ]]; then
            echo "Missing input BAM: $CTR_BAM — skipping ${cell}_${cond}_Rep${rep}"
            return 0
        fi

        # Narrow peaks
        local NAR_OUT="$OUTPUT_MACS2_PERMISSIVE/${PREFIX}_peaks.narrowPeak"
        if [[ -f "$NAR_OUT" ]]; then
            echo "Permissive narrow peaks exist for ${cell}_${cond}_Rep${rep}, skipping."
        else
            echo "Running MACS2 permissive narrow for ${cell}_${cond}_Rep${rep}..."
            macs2 callpeak -t "$IN_BAM" -c "$CTR_BAM" -f BAMPE -g "$EFFECTIVE_GENOME" -p "$PERMISSIVE_PVAL" -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "${PREFIX}" || { echo "MACS2 narrow failed for ${PREFIX}"; return 1; }
        fi

        # Broad peaks
        local BRO_OUT="$OUTPUT_MACS2_PERMISSIVE/${PREFIX}_peaks.broadPeak"
        if [[ -f "$BRO_OUT" ]]; then
            echo "Permissive broad peaks exist for ${cell}_${cond}_Rep${rep}, skipping."
        else
            echo "Running MACS2 permissive broad for ${cell}_${cond}_Rep${rep}..."
            macs2 callpeak -t "$IN_BAM" -c "$CTR_BAM" -f BAMPE -g "$EFFECTIVE_GENOME" -p "$PERMISSIVE_PVAL" --broad -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "${PREFIX}" || { echo "MACS2 broad failed for ${PREFIX}"; return 1; }
        fi

        # Sort peak files if present
        if [[ -f "$NAR_OUT" ]]; then
            sort -k8,8nr "$NAR_OUT" > "$OUTPUT_MACS2_SORT/Sort_${PREFIX}_peaks.narrowPeak"
        fi
        if [[ -f "$BRO_OUT" ]]; then
            sort -k8,8nr "$BRO_OUT" > "$OUTPUT_MACS2_SORT/Sort_${PREFIX}_peaks.broadPeak"
        fi

        echo "Completed permissive MACS2 for ${cell}_${cond}_Rep${rep}"
    }
    export -f macs2_permissive_rep
    export INPUT_CHIP_SUB OUTPUT_MACS2_PERMISSIVE OUTPUT_MACS2_SORT EFFECTIVE_GENOME PERMISSIVE_PVAL

    # Launch permissive runs in parallel
    parallel --ungroup --line-buffer -j "$JOBS" macs2_permissive_rep ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

else
    echo "Skip Permissive MACS2 run"
fi

###################################
# MACS2 merged (stringent) calls  #
###################################

echo 
echo "Running MACS2 peak calling"
macs2_merged() {
    local cell=$1
    local cond=$2

    # Set inpute files
    local IN_BAM_MERGED="$INPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd.bam"
    local CTR_BAM_MERGED="$INPUT_CHIP_ALIGN/BLF_ChIP_${cell}_controlInput_merged_cle_sort_dd.bam"
    local PREFIX="BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd_p${STRICT_PVAL}macs2"
    
    # Check inputs
    if [[ ! -f "$IN_BAM_MERGED" ]]; then
        echo "Missing input BAM: $IN_BAM_MERGED — skipping ${cell}_${cond}_merged"
        return 0
    fi

    # Check control
    if [[ ! -f "$CTR_BAM_MERGED" ]]; then
        echo "Missing input BAM: $CTR_BAM_MERGED — skipping ${cell}_${cond}_merged"
        return 0
    fi

    # Narrow peaks
    local NAR_OUT="$OUTPUT_MACS2/${PREFIX}_peaks.narrowPeak"
    if [[ -f "$NAR_OUT" ]]; then
        echo "Narrow peaks exist for ${cell}_${cond}_merged, skipping."
    else
        echo "Running MACS2 narrow for ${cell}_${cond}_merged..."
        macs2 callpeak -t "$IN_BAM_MERGED" -c "$CTR_BAM_MERGED" -f BAMPE -g "$EFFECTIVE_GENOME" -p "$STRICT_PVAL" -B --outdir "$OUTPUT_MACS2" -n "${PREFIX}" || { echo "MACS2 narrow failed for ${PREFIX}"; return 1; }
    fi

    # Broad peaks
    local BRO_OUT="$OUTPUT_MACS2/${PREFIX}_peaks.broadPeak"
    if [[ -f "$BRO_OUT" ]]; then
        echo "Broad peaks exist for ${cell}_${cond}_merged, skipping."
    else
        echo "Running MACS2 broad for ${cell}_${cond}_merged..."
        macs2 callpeak -t "$IN_BAM_MERGED" -c "$CTR_BAM_MERGED" -f BAMPE -g "$EFFECTIVE_GENOME" -p "$STRICT_PVAL" --broad -B --outdir "$OUTPUT_MACS2" -n "${PREFIX}" || { echo "MACS2 broad failed for ${PREFIX}"; return 1; }
    fi

    # Sort peak files if present
    if [[ -f "$NAR_OUT" ]]; then
        sort -k8,8nr "$NAR_OUT" > "$OUTPUT_MACS2_SORT/Sort_${PREFIX}_peaks.narrowPeak"
    fi
    if [[ -f "$BRO_OUT" ]]; then
        sort -k8,8nr "$BRO_OUT" > "$OUTPUT_MACS2_SORT/Sort_${PREFIX}_peaks.broadPeak"
    fi

    echo "Completed MACS2 peak calling for ${cell}_${cond}_merged"
}

export -f macs2_merged
export INPUT_CHIP_ALIGN OUTPUT_MACS2 OUTPUT_MACS2_SORT EFFECTIVE_GENOME STRICT_PVAL

# Launch permissive runs in parallel (per replicate)
parallel --ungroup --line-buffer -j "$JOBS" macs2_merged ::: "${CellLine[@]}" ::: "${conditions[@]}"


####################################
### IDR ANALYSIS: OPTIONAL #########
####################################
echo
echo "IDR Analysis"
# Promt to procceed or skip IDR
read -rp "Do you want to proceed IDR Analysis (y/n): " confirm
if [[ "$confirm" == "y" ]]; then

    macs2_IDR(){
        local cell=$1
        local cond=$2

        # Loop over peak types
        for peaktype in narrowPeak broadPeak; do
            # Set inpute files    
            local MACS2_R1="$OUTPUT_MACS2_SORT/Sort_BLF_ChIP_${cell}_${cond}_Rep1_cle_sort_dd_p${PERMISSIVE_PVAL}macs2_peaks.${peaktype}"
            local MACS2_R2="$OUTPUT_MACS2_SORT/Sort_BLF_ChIP_${cell}_${cond}_Rep2_cle_sort_dd_p${PERMISSIVE_PVAL}macs2_peaks.${peaktype}"
            local MACS2_ORACLE="$OUTPUT_MACS2_SORT/Sort_BLF_ChIP_${cell}_${cond}_merged_cle_sort_dd_p${STRICT_PVAL}macs2_peaks.${peaktype}"
           
            # Set output file name for IDR files
            local IDR_OUT="$OUTPUT_MACS2_IDR/${cell}_${cond}_cle_${peaktype}.idr"
            local IDR_LOG="$OUTPUT_MACS2_IDR/${cell}_${cond}_cle_${peaktype}.idr.log"
            
            # ---- File checks ----
            local missing=0
            for f in "$MACS2_R1" "$MACS2_R2" "$MACS2_ORACLE"; do
                if [[ ! -f "$f" ]]; then
                    echo "Missing required file: $f"
                    missing=1
                fi
            done

            if [[ $missing -eq 1 ]]; then
                echo "Skipping IDR for $cell $cond ($peaktype) due to missing files."
                continue  # Will continue with pipeline if you decide to delete files of one peaktype
            fi
            # ----------------------

            # IDR without Oracle Peak file
            echo "Running IDR (rep1 vs rep2) for $cell $cond ($peaktype)..."
            idr --samples "$MACS2_R1" "$MACS2_R2" --input-file-type "$peaktype" --rank p.value --output-file "$IDR_OUT" --plot --log-output-file "$IDR_LOG"
                                                                                
            # IDR with Oracle Peak file
            echo "Running IDR with oracle for $cell $cond ($peaktype)..."
            idr --samples "$MACS2_R1" "$MACS2_R2" --peak-list "$MACS2_ORACLE" --input-file-type "$peaktype" --rank p.value --output-file "${IDR_OUT%.idr}_withOracle.idr" --plot --log-output-file "${IDR_LOG%.idr.log}_withOracle.idr.log"
        
        done
    }    

    export -f macs2_IDR
    export OUTPUT_MACS2_SORT OUTPUT_MACS2_IDR PERMISSIVE_PVAL STRICT_PVAL

    # Launch IDR per condition
    parallel --ungroup --line-buffer -j "$JOBS" macs2_IDR ::: "${CellLine[@]}" ::: "${conditions[@]}"

else
    echo "Skipping IDR"
fi

echo "MACS2 + IDR pipeline completed"