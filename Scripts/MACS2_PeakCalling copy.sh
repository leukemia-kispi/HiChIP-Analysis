#!/usr/bin/env bash
set +e                 # Don’t exit on error — handle manually
set -o pipefail        # But fail properly on pipe errors
shopt -s nullglob # make globbing return empty array if no match

# Using this script assumes DirectoryArchitecture.sh was execute beforhand.
# This MACS2 pipeline expects inputs originating from pair-end read files.

# Check if GNU Parallel is installed otherwise exis
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

# Parameters
EFFECTIVE_GENOME=2913022398
PERMISSIVE_PVAL=0.05
MERGED_PVAL=0.000000001

echo
echo "Permissive MACS2 for later IDR Analysis"
read -rp "Do you want to proceed with Permissive MACS2 for later IDR Analysis (y/n): " confirm
if [[ "$confirm" == "y" ]]; then

    macs2_permissive_rep() {
        local cell=$1
        local cond=$2
        local rep=$3

        # Set inpute files
        IN_BAM_REP1="$INPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_Rep1_cle_sort_dd.bam"
        IN_BAM_REP2="$INPUT_CHIP_ALIGN/BLF_ChIP_${cell}_${cond}_Rep2_cle_sort_dd.bam"
        CTR_BAM_REP1="$INPUT_CHIP_ALIGN/BLF_ChIP_${cell}_controlInput_Rep1_cle_sort_dd.bam"
        CTR_BAM_REP2="$INPUT_CHIP_ALIGN/BLF_ChIP_${cell}_controlInput_Rep1_cle_sort_dd.bam"

        # Check inputs
        if [[ ! -f "$IN_BAM" ]]; then
            echo "Missing input BAM: $IN_BAM — skipping ${cell}_${cond}_Rep1"
            return
        fi

        # Narrow peaks
        local NAR_OUT="$OUTDIR/${PREFIX}_peaks.narrowPeak"
        if [[ -f "$NAR_OUT" ]]; then
            echo "Permissive narrow peaks exist for ${cell}_${cond}_Rep${rep}, skipping."
        else
            echo "Running MACS2 permissive narrow for ${cell}_${cond}_Rep${rep}..."
            macs2 callpeak -t "$IN_BAM" -c "$CTRL_BAM" -f BAMPE -g "$EFFECTIVE_GENOME" -p "$PERMISSIVE_PVAL" -B --outdir "$OUTDIR" -n "${PREFIX}" || { echo "MACS2 narrow failed for ${PREFIX}"; return 1; }
        fi

        # Broad peaks
        local BRO_OUT="$OUTDIR/${PREFIX}_peaks.broadPeak"
        if [[ -f "$BRO_OUT" ]]; then
            echo "Permissive broad peaks exist for ${cell}_${cond}_Rep${rep}, skipping."
        else
            echo "Running MACS2 permissive broad for ${cell}_${cond}_Rep${rep}..."
            macs2 callpeak -t "$IN_BAM" -c "$CTRL_BAM" -f BAMPE -g "$EFFECTIVE_GENOME" -p "$PERMISSIVE_PVAL" --broad -B --outdir "$OUTDIR" -n "${PREFIX}" || { echo "MACS2 broad failed for ${PREFIX}"; return 1; }
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

    # Launch permissive runs in parallel (per replicate)
    parallel --ungroup --line-buffer -j "$JOBS" macs2_permissive_rep ::: "${CellLine[@]}" ::: "${conditions[@]}" ::: "${NUMBERS[@]}"

else
    echo "Skip Permissive MACS2 run"
fi

IN_BAM_MERGED="$"
CTR_BAM_MERGED="BLF_ChIP_${cell}_control_Input_merge_cle_sort_dd.bam"

echo
echo "Permissive MACS2 for later IDR Analysis"
# Promt to procceed or skip Permissive MACS2
read -rp "Do you want to proceed with Permissive MACS2 for later IDR Analysis (y/n): " confirm
if [[ "$confirm" == "y" ]]; then   
    for cell in "${CellLine[@]}"; do 
        for cond in "${conditions[@]}"; do
            # Set inpute files
            MACS2_INPUT_R1="BLF_ChIP_${cell}_${cond}_Rep1_cle_sort_dd.bam"
            MACS2_INPUT_R2="BLF_ChIP_${cell}_${cond}_Rep2_cle_sort_dd.bam"
            CONTROL_R1="BLF_ChIP_${cell}_control_Input_Rep1_cle_sort_dd.bam"
            CONTROL_R2="BLF_ChIP_${cell}_control_Input_Rep2_cle_sort_dd.bam"
            CONTROL_MERGED="BLF_ChIP_${cell}_control_Input_merge_cle_sort_dd.bam"

            # Set output file name for peak files. MACS2 adds its own suffix   
            MACS2_PEAK_R1="${cell}_${cond}_Rep1_cle_sort_dd_p0.05macs2"
            MACS2_PEAK_R2="${cell}_${cond}_Rep2_cle_sort_dd_p0.05macs2"
                
            #Call Peaks with permissive settings <p 0.05> to be used for IDR.
            macs2 callpeak -t "$OUTPUT_CHIP_SUB/$MACS2_INPUT_R1" -c "$OUTPUT_CHIP_SUB/$CONTROL_R1" -f BAMPE  -g 2913022398 -p 0.05 -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$MACS2_PEAK_R1"
            macs2 callpeak -t "$OUTPUT_CHIP_SUB/$MACS2_INPUT_R1" -c "$OUTPUT_CHIP_SUB/$CONTROL_R1" -f BAMPE  -g 2913022398 -p 0.05 --broad -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$MACS2_PEAK_R1"
            macs2 callpeak -t "$OUTPUT_CHIP_SUB/$MACS2_INPUT_R2" -c "$OUTPUT_CHIP_SUB/$CONTROL_R2" -f BAMPE  -g 2913022398 -p 0.05 -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$MACS2_PEAK_R2"
            macs2 callpeak -t "$OUTPUT_CHIP_SUB/$MACS2_INPUT_R2" -c "$OUTPUT_CHIP_SUB/$CONTROL_R2" -f BAMPE  -g 2913022398 -p 0.05 --broad -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$MACS2_PEAK_R2"

            #Sort Peak files
            sort -k8,8nr "$OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R1}_peaks.narrowPeak" > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.narrowPeak"
            sort -k8,8nr "$OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R1}_peaks.broadPeak" > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.broadPeak"
            sort -k8,8nr "$OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R2}_peaks.narrowPeak" > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.narrowPeak"
            sort -k8,8nr "$OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R2}_peaks.broadPeak" > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.broadPeak"

            echo "Completed Permissive MACS2 for ${cell} ${cond} "
        done
    done
else
    echo "Skip Permissive MACS2 run"
fi

for cell in "${CellLine[@]}"; do 
    for cond in "${conditions[@]}"; do
        # Set inpute files
        MACS2_INPUT="BLF_ChIP_${cell}_${cond}_merge_cle_sort_dd.bam"
        CONTROL_MERGED="BLF_ChIP_${cell}_control_Input_merged_cle_sort_dd.bam"

        # Set output file name for peak files. MACS2 adds its own suffix   
        MACS2_PEAK="${cell}_${cond}_merged_cle_sort_dd_p0.05macs2"

        #Enter Merged Folder
        macs2 callpeak -t  "$OUTPUT_CHIP_ALIGN/$MACS2_INPUT" -c  "$OUTPUT_CHIP_ALIGN/$CONTROL_MERGED" -f BAMPE -g 2913022398  -p 0.000000001 -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$MACS2_PEAK"
        macs2 callpeak -t  "$OUTPUT_CHIP_ALIGN/$MACS2_INPUT" -c  "$OUTPUT_CHIP_ALIGN/$CONTROL_MERGED" -f BAMPE -g 2913022398  -p 0.000000001 --broad -B --outdir "$OUTPUT_MACS2_PERMISSIVE" -n "$MACS2_PEAK"

        #Sort Peak files
        sort -k8,8nr "$OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK}_peaks.narrowPeak" > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK}_peaks.narrowPeak"
        sort -k8,8nr "$OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK}_peaks.broadPeak" > "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK}_peaks.broadPeak"

    done
done

####################################
### IDR ANALYSIS: OPTIONAL #########
####################################
echo
echo "IDR Analysis"
# Promt to procceed or skip IDR
read -rp "Do you want to proceed IDR Analysis (y/n): " confirm
if [[ "$confirm" == "y" ]]; then  
    for cell in "${CellLine[@]}"; do 
        for cond in "${conditions[@]}"; do
            # Set input file names
            MACS2_PEAK_R1="${cell}_${cond}_Rep1_cle_sort_dd_p0.05macs2"
            MACS2_PEAK_R2="${cell}_${cond}_Rep2_cle_sort_dd_p0.05macs2"
            MACS2_PEAK_ORACLE="${cell}_${cond}_merged_cle_sort_dd_p0.05macs2"

            # Set output file name for IDR files
            IDR_OUTPUT="${cell}_${cond}"

            for peaktype in "narrowPeak" "broadPeak"; do
                # IDR without Oracle Peak file
                idr --samples "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.${peaktype}" "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.${peaktype}" --input-file-type "$peaktype" --rank p.value --output-file "$OUTPUT_MACS2_IDR/${IDR_OUTPUT}_cle_${peaktype}-idr" --plot --log-output-file "$OUTPUT_MACS2_IDR/${IDR_OUTPUT}_cle_${peaktype}.idr.log"
                                                                                
                # IDR with Oracle Peak file
                idr --samples "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.${peaktype}" "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.${peaktype}" --peak-list "$OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_ORACLE}_peaks.${peaktype}" --input-file-type "$peaktype" --rank p.value --output-file "$OUTPUT_MACS2_IDR/Oracle_${IDR_OUTPUT}_cle_${peaktype}-idr" --plot --log-output-file "$OUTPUT_MACS2_IDR/Oracle_${IDR_OUTPUT}_cle_${peaktype}.idr.log"
            done
        done
    done
else
    echo "Skip IDR"
fi