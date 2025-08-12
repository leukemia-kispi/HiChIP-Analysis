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

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="$MAIN_DIR/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="$MAIN_DIR/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="$MAIN_DIR/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for read files, *fastq.gz and ID mapping file
FASTQ_DIR="$MAIN_DIR/1.RawData"
MAPPING_FILE="$MAIN_DIR/1.RawData/Histone_ChIP_IDs.txt"
# Set output directories
OUTPUT_DIR_TRIM="$MAIN_DIR/3.TRIM/ChIP"
OUTPUT_CHIP_ALIGN="$MAIN_DIR/4.ChIP_Alignment"
OUTPUT_CHIP_SUB="$MAIN_DIR/4.ChIP_Alignment/Outputs"
BIGWIG_Coverage="$MAIN_DIR/7.Deeptool_Matrix/Coverage/"

#############################
### MAPPING NEW IDs #########
#############################

#Assuming the fastqc files were for HAL-01 ChIP-seq for histone marks, downloaded from the European Nucleotide Archive deposited under accession number ERP109232. Use the Histone_ChIP_IDs.txt file to update the annotations.
# Accession IDs = acc
# Paired_Read = read_num
# True sample name = newname
# Assumes MAPPING_FILES is a .txt wile with mapping_file format (acc newname):
#     ERR2618839  ChIP_HAL01_H3K27ac_Rep1
#     ERR2618840  ChIP_HAL01_H3K27ac_Rep2

echo
echo "=== DRY RUN: Planned renames ==="
while read -r acc newname; do
    # Skip empty lines or lines starting with #
    [[ -z "$acc" || "$acc" == \#* ]] && continue

    for read_num in 1 2; do
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
    while read -r acc newname; do
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

##########################################################
### INITIALIZE DOVETAILHICHIP CONDA ENVIORONMENT #########
##########################################################

#Environment should have following installed: bwa-mem, samtools, pairtools, fastqc, bedtools, multiqc, trim-galore, deeptools

# Initialize Conda
eval "$(conda shell.bash hook)"

# Activate Conda Environment named DovetailHiChIP to use trim-galore
#CONDA_ENV="DovetailHiChIP"
#if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
#    conda activate $CONDA_ENV
#fi

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

# Array containing replicate numbers and conditions found in filenames generated with trim-galore.
NUMBERS=("1" "2") # Replace with your actual replicate numbers
conditions=("control_Input" "H3K27ac" "H3K27me3" "H3K4me1" "H3K4me3") #For merging the replicates associated with each conditions and generating the final merged BAM files

######################
### TRIMMING #########
######################

# Flag to check if trimming needs to be performed initially set to false
perform_trimming=false

# Loop through each pair of FASTQ files if working with paired-end read files
for num in "${NUMBERS[@]}"; do
    # Check if trimmed files already exist for the current replicate
    if [ ! -f "$OUTPUT_DIR_TRIM/*rep${num}_R1_val_1.fq.gz" ] || [ ! -f "$OUTPUT_DIR_TRIM/*rep${num}_R2_val_2.fq.gz" ]; then
        perform_trimming=true
        break  # No need to check other replicates once one is found missing
    fi
done

# Perform trimming only if the flag is set to true
if [ "$perform_trimming" = true ]; then
    for num in "${NUMBERS[@]}"; do
        # Set path to input FASTQ files using wildcard pattern
        READ1="$FASTQ_DIR/*rep${num}_R1.fastq.gz"
        READ2="$FASTQ_DIR/*rep${num}_R2.fastq.gz"
       
        echo "Trimming pair $READ1 and $Read2"
        
        # Trim samples and generate new FastQC files for all replicates
        trim_galore --fastqc --phred33 --length 30 --output_dir $OUTPUT_DIR_TRIM -j 4 --paired $READ1 $READ2

        echo "Trimming for sample $num completed."
    done
else
    echo "Trimming not needed as output files already exist."
fi

#######################
### ALIGNMENT #########
#######################

#Loop through each pair of trimmed FASTQ files to complete alignment
for cond in "${conditions[@]}"; do
    for num in "${NUMBERS[@]}"; do
        # Set path to trimmed FASTQ files usign wildcard pattern
        TRIMMED_R1_FILES="$OUTPUT_DIR_TRIM/*rep${num}_R1_val_1.fq.gz"
        TRIMMED_R2_FILES="$OUTPUT_DIR_TRIM/*rep${num}_R2_val_2.fq.gz"
        # Set output file name for BAM files.
        MAPPED_BAM="ChIP_${cond}_rep${num}_cle_sort.bam"
  
        #Create sorted BAM files with grep to remove alignments to alternative contigs, unlocalized sequence, or unplaced sequence.#####################
        bwa mem -5 -T25 -t32 $REF_FASTA $TRIMMED_R1_FILES $TRIMMED_R2_FILES | samtools view -hS | grep -v chrUn | grep -v random | grep -v _alt | samtools view -bS -@32 | samtools sort -@32 -o $OUTPUT_CHIP_ALIGN/$MAPPED_BAM 

        # Set output file name for BAM files. BLF referse to black list filtered file
        MAPPED_BLF_BAM="BLF_ChIP_${cond}_rep${num}_cle_sort.bam"

        #Remove BlackList region & index
        bedtools intersect -v -abam $OUTPUT_CHIP_ALIGN/$MAPPED_BAM -b $BLACKLIST > $OUTPUT_CHIP_ALIGN/$MAPPED_BLF_BAM
        samtools index $OUTPUT_CHIP_ALIGN/$MAPPED_BLF_BAM
    
        echo " Alignment, blacklist removal and indexing done for sample $MAPPED_BLF_BAM"
    done
done

##########################################################
### INITIALIZE PICARD CONDA ENVIORONMENT #################
##########################################################

# Activate Conda Environment named Picard
CONDA_ENV="Picard"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

###################################
### DEDUPLICATION: PICARD #########
###################################

for num in "${NUMBERS[@]}"; do
    # Set output file name for dedup files.
    MAPPED_BLF_BAM_DUPFLAG="BLF_*_rep${num}_cle_sort_dupsflag.bam"
    MAPPED_BLF_TXT_DUP="BLF_*_rep${num}_cle_sort_dups.txt"
   
    #Picard MarkDuplicates and Remove
    java -jar /home/ubuntu/miniconda3/envs/Picard/share/picard-2.25.7-0/picard.jar MarkDuplicates -I $OUTPUT_CHIP_ALIGN/$MAPPED_BLF_BAM -O $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DUPFLAG -M $OUTPUT_CHIP_SUB/$MAPPED_BLF_TXT_DUP --REMOVE_DUPLICATES false

done

#################################################################
### INITIALIZE DOVETAILHICHIP CONDA ENVIORONMENT ################
#################################################################

# Activate Conda Environment named DovetailHiChIP
CONDA_ENV="DovetailHiChIP"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

###########################################
### DEDUPLICATION: DOVETAILHICHIP #########
###########################################

for num in "${NUMBERS[@]}"; do
    # Set output file name for dedup files.
    MAPPED_BLF_BAM_DD="BLF_*_rep${num}_cle_sort_dd.bam"

    # Generate final deduplicated BAM files (remove flagged duplicates) and indexing
    samtools view -b -F 1024 $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DUPFLAG > $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD
    samtools index $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD
done

conditions=("control_Input" "H3K27ac" "H3K27me3" "H3K4me1" "H3K4me3") #For merging the replicates associated with each conditions and generating the final merged BAM files

# Set output file name for dedup merged files.
for cond in "${conditions[@]}"; do
    # Set file name for black list filterd, cleaned, sorted and deduplicated replicate BAM files
    MAPPED_BLF_BAM_DD_REP1="BLF_*_${cond}_rep1_cle_sort_dd.bam" 
    MAPPED_BLF_BAM_DD_REP2="BLF_*_${cond}_rep2_cle_sort_dd.bam" 

    # Set output file name for merged  BAM files
    MAPPED_MERGED_BLF_BAM_DD="BLF_*_${cond}_merged_cle_sort_dd.bam"

    # MERGE Replicate File and index
    samtools merge -f $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD_REP1 $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD_REP2   ####have to update file annotations 
    samtools index $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD
done


##############################
### GENOMIC COVERAGE #########
##############################

BIGWIG_BLF_DD="BLF_*_rep${num}_cle_sort_dd.bw"
BIGWIG_MERGED_BLF_DD="BLF_*_merged_cle_sort_dd.bw"

for num in "${NUMBERS[@]}"; do
# Generate Enrichment coverage files for IGV visualization for each replicate
bamCoverage -b $OUTPUT_CHIP_SUB/$MAPPED_BLF_BAM_DD -o $BIGWIG_Coverage/$BIGWIG_BLF_DD --effectiveGenomeSize 2913022398 -bl $BLACKLIST--normalizeUsing RPKM -p max -bs 10 \
--extendReads --ignoreForNormalization M

for cond in "${conditions[@]}"; do
# Generate Enrichment coverage files for IGV visualization for merged files
bamCoverage -b  $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD -o $BIGWIG_Coverage/$BIGWIG_MERGED_BLF_DD --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 \
--extendReads --ignoreForNormalization M

################################
### MACS2 PEAK CALLING #########
################################

# Activate Conda Environment named macs2
CONDA_ENV="macs2"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

#Call Peaks with permissive settings to be used for IDR.
macs2 callpeak -t BLF_ChIP_HAL01_H3K27ac_Rep1_cle_sort_dd.bam -c BLF_ChIP_HAL01_control_Input_Rep1_cle_sort_dd.bam -f BAMPE -g hs -p 0.05 -B --outdir /mnt/RawChIP_Histon/MACS2 -n BLF_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2


macs2 callpeak -t BLF_ChIP_HAL01_H3K27ac_Rep2_cle_sort_dd.bam -c BLF_ChIP_HAL01_control_Input_Rep2_cle_sort_dd.bam -f BAMPE -g hs -p 0.05 -B --outdir /mnt/RawChIP_Histon/MACS2 -n BLF_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05macs2


#Enter Merged Folder
macs2 callpeak -t  $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD -c  $OUTPUT_CHIP_ALIGN/$MAPPED_MERGED_BLF_BAM_DD_INPUT -f BAMPE -g hs -p 0.000000001 -B --outdir /mnt/RawChIP_Histon/MACS2 -n HAL01_H3K27ac_merged_cle_sort_dd_p9macs2

macs2 callpeak -t $OUTPUT_HICHIP_SUB/$PRIMARY_ALN --keep-dup 10 --min-length 300 -p 0.000000001 -g 2913022398 --bw 300 --mfold 5 50 -n $OUTPUT_MACS2/$MACS2_JoinedRep_Oracle


#Sort Peak files

sort -k8,8nr BLF_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak > /mnt/RawChIP_Histon/IDR/Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak

####################################
### IDR ANALYSIS: OPTIONAL #########
####################################

#IDR without Oracle Peak file
idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05macs2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file HAL-01_H3K27ac_cle-idr --plot --log-output-file HAL-01_H3K27ac_cle.idr.log

idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05bmacs2_peaks.broadPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05bmacs2_peaks.broadPeak  --input-file-type broadPeak --rank p.value --output-file HAL-01_H3K27ac_cle_broad-idr --plot --log-output-file HAL-01_H3K27ac_cle_broad.idr.log
	                                                                     

#IDR with Oracle Peak file
idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05macs2_peaks.narrowPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05macs2_peaks.narrowPeak --peak-list Sort_HAL01_H3K27ac_merged_cle_sort_dd_p9macs2_peaks.narrowPeak  --input-file-type narrowPeak --rank p.value --output-file Oracle_HAL-01_H3K27ac_cle-idr --plot --log-output-file Oracle_HAL-01_H3K27ac_cle.idr.log

idr --samples Sort_HAL01_H3K27ac_Rep1_cle_sort_dd_p0.05bmacs2_peaks.broadPeak Sort_HAL01_H3K27ac_Rep2_cle_sort_dd_p0.05bmacs2_peaks.broadPeak --peak-list Sort_HAL01_H3K27ac_merged_cle_sort_dd_p9bmacs2_peaks.broadPeak  --input-file-type broadPeak --rank p.value --output-file Oracle_HAL-01_H3K27ac_cle_broad-idr --plot --log-output-file Oracle_HAL-01_H3K27ac_cle_broad.idr.log