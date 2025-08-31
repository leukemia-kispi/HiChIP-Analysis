#!/usr/bin/env bash
set -e
shopt -s nullglob # make globbing return empty array if no match

#Using this script assumes the script DirectoryArchitecture&CondaEnv.sh was run beforhand to create the DirectoryArchitecture and generate conda environments with needed tools.

#Set path to reference genome index and blacklist. Genome index has to be generated first if not done. 
REF_FASTA="/mnt/0.GenomeAssembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" 
REF_GENOME="/mnt/0.GenomeAssembly/GRCh38_no_alt_ref.genome"
BLACKLIST="/mnt/0.BlackList/hg38-blacklist.v2.bed"
# Set Path for read before and after trimming, *fg.gz files
FASTQ_DIR="/mnt/1.RawData"
# Set output directories
OUTPUT_DIR_TRIM="/mnt/3.TRIM"
OUTPUT_HICHIP_ALIGN="/mnt/4.HiChIP_Alignment"
OUTPUT_HICHIP_SUB="/mnt/4.HiChIP_Alignment/Outputs"
# Thread usage
cores=32
#Thread usage for pairtools dedup and split processes
cores2=16
# Set Path to temporary directory
TEMP="/mnt/tmp"

# Initialize Conda
eval "$(conda shell.bash hook)"

# Activate Conda Environment named DovetailHiChIP
CONDA_ENV="DovetailHiChIP"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

# Array containing replicate numbers found in filenames generated with TrimGalore.
NUMBERS=("1" "2") # Replace with your actual replicate numbers

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
        READ1="$FASTQ_DIR/*_rep${num}_R1.fastq.gz"
        READ2="$FASTQ_DIR/*_rep${num}_R2.fastq.gz"

        # Trim samples and generate new FastQC files for all replicates
        trim_galore --fastqc --phred33 --length 50 --output_dir $OUTPUT_DIR_TRIM -j 4 --paired $READ1 $READ2

        echo "Trimming for sample $num completed."
    done
else
    echo "Trimming not needed as output files already exist."
fi

# Alignment Output directory
cd $OUTPUT_HICHIP_ALIGN

# Loop through each pair of FASTQ files if working with paired-end read files for SingleRep HiChIP alignment
for num in "${NUMBERS[@]}"; do
    # Set path to input FASTQ files using wildcard pattern
    HICHIP_R1="$OUTPUT_DIR_TRIM/*rep${num}_R1_val_1.fq.gz"
    HICHIP_R2="$OUTPUT_DIR_TRIM/*rep${num}_R2_val_2.fq.gz"
    MAPPED_PAIRS="Rep${num}_TCF3_HLF_hg38_nodd_mapped.pairs"
    MAPPED_BAM="Rep${num}_TCF3_HLF_hg38_nodd_mapped.PT.bam"
    MAPPED_BLF_BAM="BLF_Rep${num}_TCF3_HLF_hg38_nodd_mapped.PT.bam"

    bwa mem -5SP -T0 -t$cores $REF_FASTA $HICHIP_R1 $HICHIP_R2 | \
    pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores2 --nproc-out $cores2 --chroms-path $REF_GENOME | \
    pairtools sort --tmpdir=$TEMP --nproc $cores | \
    pairtools split --nproc-in $cores2 --nproc-out $cores2 --output-pairs $MAPPED_PAIRS --output-sam -|\
    samtools view -bS -@$cores | \
    samtools sort -@$cores -o $MAPPED_BAM;samtools index $MAPPED_BAM

    echo "HiCHIP Aligmnent Complete for rep${num}"

    cd /home/ubuntu

    #Remove black list
    bedtools intersect -v -abam $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -b $BLACKLIST > $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM
    samtools index $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM

    #QC compare ChIP-seq TCF3-HLF_FLAG
    bash /home/ubuntu/HiChiP/enrichment_stats.sh -g $REF_FASTA -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -p /home/ubuntu/HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -t $cores2 -x $OUTPUT_HICHIP_SUB/HiChIP_rep${num}_TCF3HLF_bw175

    #QC Plot ChIP-seq TCF3-HLF_FLAG
    python3 /home/ubuntu/HiChiP/plot_chip_enrichment_bed.py -bam $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -peaks /home/ubuntu/HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -output $OUTPUT_HICHIP_SUB/HiChIP_rep${num}_TCF3HLF_ChIP_FLAG_bw175_enrichment.png

    echo "HiCHIP Aligmnent QC Complete for rep${num}"
    
    #Enrichment for IGV
    bamCoverage -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BAM -o $OUTPUT_HICHIP_SUB/BLF_rep${num}_TCF3HLF_hg38_nodd_mapped.bw --effectiveGenomeSize 2913022398 -bl $BLACKLIST --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

    echo "Generated Bigwig file Complete for rep${num}"

    #ContacMaps
    java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/HiChiP/juicer_tools_1.22.01.jar pre --threads $cores $OUTPUT_HICHIP_ALIGN/$MAPPED_PAIRS $OUTPUT_HICHIP_SUB/rep${num}_TCF3HLF_HAL01_hg38_nodd_contact_map.hic $REF_GENOME

    echo "Generated .hic file Complete rep${num}"
done

####################
### FILTER PAIRS ###
####################

# The "paritools select" command in the DovtailHiChIP conda environment will filter
# and include unique and rescue mapped read pairs.

for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do

        # Input files
        MAPPED_PAIRS="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_mapped.pairs"

        # Output files
        MAPPED_PAIRS_FILTERED="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_mapped.filtered.pairs"

        pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' "$MAPPED_PAIRS" -o "$MAPPED_PAIRS_FILTERED"
    done
done

################################
### HICPRO VALID PAIRS FILES ### 
################################

# The paired files generated during aligment needs to be converted
# into valid HiC-Pro files to proceed with FitHiChIP.

for cell in "${CellLine[@]}"; do    
    for cond in "${conditions[@]}"; do

        # Input files
        MAPPED_PAIRS_FILTERED="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_mapped.filtered.pairs"

        # Output files
        MAPPED_PAIRS_HICPRO="$OUTPUT_HICHIP_ALIGN/Merge_${cell}_${cond}_hg38_nodd_hicpro_mapped.filtered.pairs.gz"

        grep -v '#' "$MAPPED_PAIRS"| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > "$MAPPED_PAIRS_HICPRO"
        
        echo "Converted $MAPPED_PAIRS into HiC-Pro format"
   done
done