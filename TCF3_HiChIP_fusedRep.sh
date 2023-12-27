#Set path to reference genome index. Genome index has to be generated first if not done. 
GENOME_INDEX="/mnt/0.STAR_hg38_v44_index"
# Set output directory
OUTPUT_DIR="/mnt/4.ALIGNMENT"
# Set Path for read *fg.gz files
FASTQ_DIR="/mnt/3.TRIM"

#TrimGalore
conda activate trimgalore

#For all replicates
trim_galore --fastqc --phred33 --length 50 --output_dir /mnt/RawHiChIP_trimgalor/TCF3-HLF/ -j 4 --paired /mnt/RawHiChIP/TCF3-HLF/20211112.A-o26422_1_1-HiChIP_HAL-01_aTCF3_rep1_R1.fastq.gz /mnt/RawHiChIP/TCF3-HLF/20211112.A-o26422_1_1-HiChIP_HAL-01_aTCF3_rep1_R2.fastq.gz 

trim_galore --fastqc --phred33 --length 50 --output_dir /mnt/RawHiChIP_trimgalor/TCF3-HLF/ -j 4 --paired /mnt/RawHiChIP/TCF3-HLF/20211112.A-o26422_1_2-HiChIP_HAL-01_aTCF3_rep2_R1.fastq.gz /mnt/RawHiChIP/TCF3-HLF/20211112.A-o26422_1_2-HiChIP_HAL-01_aTCF3_rep2_R2.fastq.gz

#Fuse Fasta files

cat 20211112.A-o26422_1_1-HiChIP_HAL-01_aTCF3_rep1_R2_val_1.fq.gz 20211112.A-o26422_1_2-HiChIP_HAL-01_aTCF3_rep2_R2_val_1.fq.gz > JoinedFastq_R1.fq.gz
cat 20211112.A-o26422_1_1-HiChIP_HAL-01_aTCF3_rep1_R2_val_2.fq.gz 20211112.A-o26422_1_2-HiChIP_HAL-01_aTCF3_rep2_R2_val_2.fq.gz > JoinedFastq_R2.fq.gz

# Alignment
bwa mem -5SP -T0 -t32 /mnt/GRCh38_assembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /mnt/RawHiChIP_trimgalor/TCF3-HLF/JoinedFastq_R1.fq.gz /mnt/RawHiChIP_trimgalor/TCF3-HLF/JoinedFastq_R2.fq.gz | pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 16 --nproc-out 16 --chroms-path /mnt/GRCh38_assembly/GRCh38_no_alt_ref.genome | pairtools sort --tmpdir=/mnt/temp/ --nproc 32 | pairtools split --nproc-in 16 --nproc-out 16 --output-pairs JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs --output-sam -| samtools view -bS -@32 | samtools sort -@32 -o JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam;samtools index JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam

#QC compare ChIP-seq TCF3-HLF_FLAG
bash enrichment_stats.sh -g /mnt/GRCh38_assembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -b /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -p /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_cle-idr.bed -t 16 -x HiChIP-TCF3-HLF
bash enrichment_stats.sh -g /mnt/GRCh38_assembly/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -b /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -p /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -t 16 -x HiChIP-TCF3-HLF_bw175

#QC Plot ChIP-seq TCF3-HLF_FLAG
python3 plot_chip_enrichment_bed.py -bam /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -peaks /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_cle-idr.bed -output HiChIP_TCF3_HLF_ChIP_FLAG_enrichment.png
python3 plot_chip_enrichment_bed.py -bam /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -peaks /mnt/ChIP-Seq/IDR_cle/Oracle2_HAL-01_TCF3-HLF_FLAG_bw175_cle-idr.bed -output HiChIP_TCF3_HLF_ChIP_FLAG_bw175_enrichment.png

#Enrichment for IGV
conda activate deeptool
bamCoverage -b JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -o JoinedRep_TCF3_HLF_hg38_nodd_mapped.bw -bl /mnt/BlastListPeaks/hg38-blacklist.v2.bed --normalizeUsing RPKM -p max -bs 10
bamCoverage -b JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -o BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.bw --effectiveGenomeSize 2913022398 -bl /mnt/BlastListPeaks/hg38-blacklist.v2.bed --normalizeUsing RPKM -p max -bs 10 --extendReads --ignoreForNormalization M

#ContacMaps
conda activate Juicebox
java -Xmx48000m  -Djava.awt.headless=true -jar /home/ubuntu/My_HiC-Pro/juicer_tools_1.22.01.jar pre --threads 32 JoinedRep_TCF3_HLF_hg38_nodd_mapped.pairs JoinedRep_TCF3-HLF_HAL01_hg38_nodd_contact_map.hic /mnt/GRCh38_assembly/GRCh38_no_alt_ref.genome

#call 1D peaks with MACS2
conda activate MACS2

#Remove black list
bedtools intersect -v -abam JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -b /mnt/BlastListPeaks/hg38-blacklist.v2.bed > BLF_JoinedRep_TCF3-HLF_HAL01_hg38_nodd_mapped.PT.bam

#MACS2 Peak calling on BlackList filtered BAM.
samtools view -b BLF_JoinedRep_TCF3_HLF_hg38_nodd_mapped.PT.bam -h -F 0x900 | bedtools bamtobed -i stdin > BLF_JoinedRep_TCF3_HLF_hg38_nodd_primary.aln.bed

macs2 callpeak -t BLF_JoinedRep_TCF3_HLF_hg38_nodd_primary.aln.bed -p 0.000000001 -n BLF_JoinedRep_TCF3_HLF_hg38_nodd_p9.macs2
macs2 callpeak -t BLF_JoinedRep_TCF3-HLF_HAL01_hg38_nodd_primary.aln.bed -p 0.000000001 -g 2913022398 -n BLF_JoinedRep_TCF3_HLF_gs_hg38_nodd_p9.macs2

#Oracle File for IDR and IDR
macs2 callpeak -t BLF_JoinedRep_TCF3-HLF_HAL01_hg38_nodd_primary.aln.bed --keep-dup 10 --min-length 300 -p 0.000000001 --bw 300 --mfold 5 50 -n BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2
sort -k8,8nr  BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2_peaks.narrowPeak > /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Oracle_MACS2/Sort_BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2_peaks.narrowPeak 
sort -k8,8nr  BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak > /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak
sort -k8,8nr  BLF_TCF3-HLF_HAL01_rep2_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak > /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep2_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak

idr --samples /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Permissive_MACS2/Sort_BLF_TCF3-HLF_HAL01_rep2_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak --peak-list /mnt/Dovetail_Pipeline/TCF3-HLF_HAL01_hg38_JoinedRep/Oracle_MACS2/Sort_BLF_JoinedRep_TCF3_HLF_hg38_nodd_kdp9Oracle.macs2_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file Oracle_HAL-01_TCF3_JoinedRep_cle-idr --plot --log-output-file Oracle_HAL-01_TCF3_JoinedRep_cle.idr.log

#FitHiChIP Loop calling
conda activate /home/ubuntu/My_HiCPro_ENV/

export PATH="/home/ubuntu/My_HiC-Pro/HiC-Pro_3.0.0/bin/":$PATH

BLF_TCF3-HLF_HAL01_rep1_hg38_nodd_kdp5bw300ml300Sort.macs2_peaks.narrowPeak
#FitHiChIP Loop Calling
#Filter pairs
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU") or (pair_type=="uu") or (pair_type=="Uu")  or (pair_type=="uU")' JoinedRep_TCF3_HLF_hg38_nodd_hicpro_mapped.pairs -o JoinedRep_TCF3-HLF_hg38_nodd_mapped.filtered.pairs

#HiCPro Valid Pairs Files
grep -v '#' JoinedRep_TCF3-HLF_hg38_nodd_mapped.filtered.pairs| awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' | gzip -c > JoinedRep_TCF3-HLF_hg38_nodd_hicpro_mapped.filtered.pairs.gz

cd /home/ubuntu/My_HiC-Pro/FitHiChIP/
bash ./FitHiChIP_Docker.sh -C /home/ubuntu/My_HiC-Pro/FitHiChIP/configfile_CB_TCF3_JoinedRep_5kb_Merge_50kb_3M