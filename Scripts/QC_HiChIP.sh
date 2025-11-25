

# Input files

# Output files


#QC compare ChIP-seq TCF3-HLF_FLAG
bash /home/$USER/HiChiP/enrichment_stats.sh -g $REF_FASTA -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -p /home/$USER/HiChIP_Analysis/Extras/ChIP_Seq_HAL01_TCF3HLF_FLAG.bed -t $cores2 -x $OUTPUT_HICHIP_SUB/HiChIPvsChIP_enrichment

#QC Plot ChIP-seq TCF3-HLF_FLAG
python3 /home/$USER/HiChiP/plot_chip_enrichment_bed.py -bam $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -peaks /home/$USER/HiChIP_Analysis/Extras/ChIP_Seq_HAL01_TCF3HLF_FLAG.bed -output $OUTPUT_HICHIP_SUB/HiChIPvsChIP_enrichment.png