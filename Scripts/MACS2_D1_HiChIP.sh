#!usr/bin/bash
set -e
shopt -s nullglob # make globbing return empty array if no match

# Using this script assumes DirectoryArchitecture.sh was execute beforhand

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

#Set path to blacklist. 
BLACKLIST="$MAIN_DIR/0.BlackList/hg38-blacklist.v2.bed"
# Set output directories
OUTPUT_HICHIP_ALIGN="$MAIN_DIR/4.Alignment/HiChIP"
OUTPUT_MACS2="$MAIN_DIR/5.MACS2/HiChIP"
OUTPUT_MACS2_SORT="$MAIN_DIR/5.MACS2/HiChIP/SORT"
OUTPUT_MACS2_PERMISSIVE="$MAIN_DIR/5.MACS2/HiChIP/Permissive"
OUTPUT_MACS2_IDR="/mnt/5.MACS2/HiChIP/IDR"

################################################
### ADJUST THESE TO MATCH SAMPLE NAMES #########
################################################

# Array containing Cell lines, replicate numbers and conditions found in filenames and defining samples.
# Expected Sample nomenclature follows this pattern HiChIP_<CellLine>_<conditions>_Rep<NUMBERS>_suffix.
CellLine=("HAL01") #Replace with your actual Cell Line 
conditions=("TCF3HLF") #Replace with your actual conditions
NUMBERS=("1" "2") # Replace with your actual replicate numbers

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