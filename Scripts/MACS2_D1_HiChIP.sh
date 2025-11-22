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
INPUT_HICHIP_ALIGN="$MAIN_DIR/4.Alignment/HiChIP"
$INPUT_HICHIP_SUB="$MAIN_DIR/4.Alignment/HiChIP/Outputs"
# Provide your own list if running different samples
BLACKLIST="$MAIN_DIR/0.BlackList/hg38-blacklist.v2.bed"
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

########################################################
### GENERATE PRIMARY ALIGNMENT FILES AND MACS2 #########
########################################################


for cell in "${CellLine[@]}"; do 
    for cond in "${conditions[@]}"; do

        #Inputs/Output files
        MAPPED_BLF_BAM_MERGE="BLF_HiChIP_${cell}_${cond}_Merge_hg38_nodd_mapped.PT.bam"
        PRIMARY_ALN_MERGE="BLF_HiChIP_${cell}_${cond}_Merge_hg38_nodd_primary.aln.bed"
        MACS2_Merge="BLF_HiChIP_${cell}_${cond}_Merge_gs_hg38_nodd_p9.macs2"

        samtools view -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -h -F 0x900 | bedtools bamtobed -i stdin > $OUTPUT_HICHIP_SUB/$PRIMARY_ALN

        echo "$PRIMARY_ALN_MERGE generated"

        macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN --keep-dup 10 --min-length 300 -p 0.000000001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2/$MACS2_Merge
        #macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN -p 0.000000001 -g 2913022398 -n $OUTPUT_MACS2/$MACS2_Merge

        echo "$MACS2_Merge peak calling done for "

    done
done

for cell in "${CellLine[@]}"; do 
    for cond in "${conditions[@]}"; do
        for num in "${NUMBERS[@]}"; do

            #Inpute/Output Files
            MAPPED_BLF_BAM_Rep="BLF_HiChIP_${cell}_${cond}_Rep{$num}_hg38_nodd_mapped.PT.bam"
            PRIMARY_ALN_Rep="BLF_HiChIP_${cell}_${cond}_Rep{$num}_hg38_nodd_primary.aln.bed"
            MACS2_REP="BLF_HiChIP_${cell}_${cond}_Rep{$num}_gs_hg38_p9.macs2"

            samtools view -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM_Rep -h -F 0x900 | bedtools bamtobed -i stdin > $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep

            echo "$MAPPED_BLF_BAM_Rep generated"

            PRIMARY_ALN_Rep="BLF_HiChIP_${cell}_${cond}_Rep{$num}_hg38_nodd_primary.aln.bed"
            
            macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep --keep-dup 10 --min-length 300 -p 0.000000001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2/$MACS2_Rep
            #macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep -p 0.000000001 -g 2913022398 -n $OUTPUT_MACS2/$MACS2_Rep

            echo "$MACS2_REP peak calling done"

        done
    done
done

################################################
### OPTIONAL: PERMISSIVE MACS2 and IDR #########
################################################

PRIMARY_ALN_Rep="BLF_HiChIP_${cell}_${cond}_Rep{$num}_hg38_nodd_primary.aln.bed"
MACS2_Rep_Permissive=""

#Oracle File for IDR 
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep --keep-dup 10 --min-length 300 -p 0.00001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2_Permissive/$MACS2_Rep_Permissive
macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN_Rep2 --keep-dup 10 --min-length 300 -p 0.00001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2_Permissive/$MACS2_Rep2_Permissive
sudo chmod 777 -R $OUTPUT_MACS2


sort -k8,8nr $OUTPUT_MACS2/$MACS2_JoinedRep_Oracle"_peaks.narrowPeak" > $OUTPUT_MACS2_SORT/$MACS2_JoinedRep_SORT
sort -k8,8nr $OUTPUT_MACS2_Permissive/$MACS2_Rep1_Permissive"_peaks.narrowPeak"  > $OUTPUT_MACS2_SORT/$MACS2_Rep1_SORT
sort -k8,8nr $OUTPUT_MACS2_Permissive/$MACS2_Rep2_Permissive"_peaks.narrowPeak"  > $OUTPUT_MACS2_SORT/$MACS2_Rep2_SORT

idr --samples $OUTPUT_MACS2_SORT/$MACS2_Rep1_SORT $OUTPUT_MACS2_SORT/$MACS2_Rep2_SORT --peak-list $OUTPUT_MACS2_SORT/$MACS2_JoinedRep_SORT --input-file-type narrowPeak --rank p.value --output-file $OUTPUT_MACS2_IDR/OraclePeaks_HAL01_TCF3_cle-idr --plot --log-output-file $OUTPUT_MACS2_IDR/OraclePeak_HAL01_TCF3_cle.idr.log

echo "All MACS2 Oracle runs Complete"