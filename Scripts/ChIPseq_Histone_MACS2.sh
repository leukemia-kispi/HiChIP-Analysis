#!usr/bin/env bash
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

# Set output directories
OUTPUT_CHIP_ALIGN="$MAIN_DIR/4.ChIP_Alignment"
OUTPUT_CHIP_SUB="$MAIN_DIR/4.ChIP_Alignment/Outputs"
OUTPUT_MACS2="$MAIN_DIR/5.MACS2"
OUTPUT_MACS2_SORT="$MAIN_DIR/5.MACS2/SORT"
OUTPUT_MACS2_PERMISSIVE="$MAIN_DIR/5.MACS2/Permissive"
OUTPUT_MACS2_IDR="$MAIN_DIR/5.MACS2/IDR"

# Array containing replicate numbers and conditions found in filenames generated with trim-galore.
NUMBERS=("1" "2") # Replace with your actual replicate numbers
conditions=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me3") #For merging the replicates associated with each conditions and generating the final merged BAM files

#################################################################
### INITIALIZE MACS2 CONDA ENVIORONMENT #########################
#################################################################

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
for cond in "${conditions[@]}"; do
    # Set inpute files
    MACS2_INPUT_R1="BLF_ChIP_HAL01_${cond}_Rep1_cle_sort_dd.bam"
    MACS2_INPUT_R2="BLF_ChIP_HAL01_${cond}_Rep2_cle_sort_dd.bam"
    CONTROL_R1="BLF_ChIP_HAL01_control_Input_Rep1_cle_sort_dd.bam"
    CONTROL_R2="BLF_ChIP_HAL01_control_Input_Rep2_cle_sort_dd.bam"
    CONTROL_MERGED="BLF_ChIP_HAL01_control_Input_merged_cle_sort_dd.bam"

    # Set output file name for peak files. MACS2 adds its own suffix   
    MACS2_PEAK_R1="HAL01_${cond}_Rep1_cle_sort_dd_p0.05macs2"
    MACS2_PEAK_R2="HAL01_${cond}_Rep2_cle_sort_dd_p0.05macs2"
        
    #Call Peaks with permissive settings to be used for IDR.
    macs2 callpeak -t $OUTPUT_CHIP_SUB/$MACS2_INPUT_R1 -c $OUTPUT_CHIP_SUB/$CONTROL_R1 -f BAMPE  -g 2913022398 -p 0.05 -B --outdir $OUTPUT_MACS2_PERMISSIVE -n $MACS2_PEAK_R1
    macs2 callpeak -t $OUTPUT_CHIP_SUB/$MACS2_INPUT_R1 -c $OUTPUT_CHIP_SUB/$CONTROL_R1 -f BAMPE  -g 2913022398 -p 0.05 --broad -B --outdir $OUTPUT_MACS2_PERMISSIVE -n $MACS2_PEAK_R1
    macs2 callpeak -t $OUTPUT_CHIP_SUB/$MACS2_INPUT_R2 -c $OUTPUT_CHIP_SUB/$CONTROL_R2 -f BAMPE  -g 2913022398 -p 0.05 -B --outdir $OUTPUT_MACS2_PERMISSIVE -n $MACS2_PEAK_R2
    macs2 callpeak -t $OUTPUT_CHIP_SUB/$MACS2_INPUT_R2 -c $OUTPUT_CHIP_SUB/$CONTROL_R2 -f BAMPE  -g 2913022398 -p 0.05 --broad -B --outdir $OUTPUT_MACS2_PERMISSIVE -n $MACS2_PEAK_R2

    #Sort Peak files
    sort -k8,8nr $OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R1}_peaks.narrowPeak > $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.narrowPeak
    sort -k8,8nr $OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R1}_peaks.broadPeak > $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.broadPeak
    sort -k8,8nr $OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R2}_peaks.narrowPeak > $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.narrowPeak
    sort -k8,8nr $OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK_R2}_peaks.broadPeak > $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.broadPeak
done

for cond in "${conditions[@]}"; do
    # Set inpute files
    MACS2_INPUT="BLF_ChIP_HAL01_${cond}_merged_cle_sort_dd.bam"
    CONTROL_MERGED="BLF_ChIP_HAL01_control_Input_merged_cle_sort_dd.bam"

    # Set output file name for peak files. MACS2 adds its own suffix   
    MACS2_PEAK="HAL01_${cond}_merged_cle_sort_dd_p0.05macs2"

    #Enter Merged Folder
    macs2 callpeak -t  $OUTPUT_CHIP_ALIGN/$MACS2_INPUT -c  $OUTPUT_CHIP_ALIGN/$CONTROL_MERGED -f BAMPE -g 2913022398  -p 0.000000001 -B --outdir $OUTPUT_MACS2_PERMISSIVE -n $MACS2_PEAK
    macs2 callpeak -t  $OUTPUT_CHIP_ALIGN/$MACS2_INPUT -c  $OUTPUT_CHIP_ALIGN/$CONTROL_MERGED -f BAMPE -g 2913022398  -p 0.000000001 --broad -B --outdir $OUTPUT_MACS2_PERMISSIVE -n $MACS2_PEAK

    #Sort Peak files
    sort -k8,8nr $OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK}_peaks.narrowPeak > $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK}_peaks.narrowPeak
    sort -k8,8nr $OUTPUT_MACS2_PERMISSIVE/${MACS2_PEAK}_peaks.broadPeak > $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK}_peaks.broadPeak

done

####################################
### IDR ANALYSIS: OPTIONAL #########
####################################

for cond in "${conditions[@]}"; do
    # Set input file names
    MACS2_PEAK_R1="HAL01_${cond}_Rep1_cle_sort_dd_p0.05macs2"
    MACS2_PEAK_R2="HAL01_${cond}_Rep2_cle_sort_dd_p0.05macs2"
    MACS2_PEAK="HAL01_${cond}_merged_cle_sort_dd_p0.05macs2"

    # Set output file name for IDR files
    IDR_OUTPUT="HAL01_${cond}"

    for peaktype in "narrowPeak" "broadPeak"; do
        # IDR without Oracle Peak file
        idr --samples $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.${peaktype} $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.${peaktype} --input-file-type $peaktype --rank p.value --output-file $OUTPUT_MACS2_IDR/${IDR_OUTPUT}_cle_broad-idr --plot --log-output-file $OUTPUT_MACS2_IDR/${IDR_OUTPUT}_cle_${peaktype}.idr.log
	                                                                     
        # IDR with Oracle Peak file
        idr --samples $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R1}_peaks.${peaktype} $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK_R2}_peaks.${peaktype} --peak-list $OUTPUT_MACS2_SORT/Sort_${MACS2_PEAK}_peaks.${peaktype} --input-file-type ${peaktype} --rank p.value --output-file $OUTPUT_MACS2_IDR/Oracle_${IDR_OUTPUT}_cle_${peaktype}-idr --plot --log-output-file $OUTPUT_MACS2_IDR/Oracle_${IDR_OUTPUT}_cle_${peaktype}.idr.log
    done
done