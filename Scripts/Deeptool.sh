#!/bin/bash

#Set path to blacklist.
BLACKLIST="/mnt/DeepTool/Target_Region/Q1cutoff_3197_BlackList.bed"
# Set output directories
OUTPUT_DEEPTOOL="/mnt/7.Deeptool_Matrix/Outputs"
OUTPUT_DEEPTOOL_GRAPHS="/mnt/7.Deeptool_Matrix/Graphs"
OUTPUT_DEEPTOOL_SORTED="/mnt/7.Deeptool_Matrix/Sorted_Lists"
#Inputs
BW_HICHIP="/mnt/7.Deeptool_Matrix/Coverage/BLF_JoinedRep_TCF3HLF_hg38_nodd_mapped.bw"
BW_H3K27AC="/mnt/7.Deeptool_Matrix/Coverage/ChIP_HAL01_H3K27ac_merged_cle_dd_extended_RPKM.bw"
BW_H3K4M1="/mnt/7.Deeptool_Matrix/Coverage/ChIP_HAL01_H3K4me1_merged_cle_dd_extended_RPKM.bw "
BW_H3K4ME3="/mnt/7.Deeptool_Matrix/Coverage/ChIP_HAL01_H3K4me3_merged_cle_dd_extended_RPKM.bw"
BW_H3K27ME3="/mnt/7.Deeptool_Matrix/Coverage/ChIP_HAL01_H3K27me3_merged_cle_dd_extended_RPKM.bw"
BW_ATACSEQ="/mnt/7.Deeptool_Matrix/Coverage/GSM5663910_05_HAL01_ATAC-Seq.hg38.bw"
BED_INPUT_WHOLE="/mnt/7.Deeptool_Matrix/Target_Regions/BLF_JoinedRep_TCF3_HLF_gs_hg38_nood_kdp9_500bp_summit.bed"
BED_INPUT_MOTIF="/mnt/7.Deeptool_Matrix/Target_Regions/Anchor_UniqueHLF_motifs_3197.bed"
BED_INPUT_CLUSTER_1_2="/mnt/7.Deeptool_Matrix/Target_Regions/Cluster_1-2.bed"
BED_INPUT_CLUSTER_3="/mnt/7.Deeptool_Matrix/Target_Regions/Cluster_3.bed"
BED_INPUT_CLUSTER_4="/mnt/7.Deeptool_Matrix/Target_Regions/Cluster_4.bed"

# Initialize Conda
eval "$(conda shell.bash hook)"

# Activate Conda Environment named DovetailHiChIP
CONDA_ENV="DovetailHiChIP"
if [[ "$(conda info --base)" != "$(conda info --base --json | jq -r .conda_prefix)" ]]; then
    conda activate $CONDA_ENV
fi

#Whole HiChIP Peak List Generated from MACS2
computeMatrix reference-point -S $BW_HICHIP \ 
$BW_H3K27AC \
$BW_H3K4M1 \
$BW_H3K4ME3 \
$BW_H3K27ME3 \
$BW_ATACSEQ \
-R $BED_INPUT_WHOLE \
--referencePoint center -b 5000 -a 5000 -bs 10 -p "max" -out $OUTPUT_DEEPTOOL/TCF3-HLF_JoinedRep_nood_p9_500bp_summit_bs10_marks.mat \
--missingDataAsZero --skipZeros

#HiChIP 3197 HLF motif Peaks. Using MEME-suite filtered out to only inlcude the 3197 peaks with HLF motif
computeMatrix reference-point -S $BW_HICHIP \ 
$BW_H3K27AC \
$BW_H3K4M1 \
$BW_H3K4ME3 \
$BW_H3K27ME3 \
$BW_ATACSEQ \
-R $BED_INPUT_MOTIF \
--referencePoint center -b 5000 -a 5000 -bs 10 -p "max" -out $OUTPUT_DEEPTOOL/TCF3HLF_HLF_motifs_3197_marks.mat \
--missingDataAsZero --skipZeros

#HiChIP 3197 HLF motif Peaks with 1Quartile cutoff. Using R to see data in TCF3HLF_HLF_motifs_3197_marks.mat and remove the 1st Quartile of data (low TCF3::HLF signal) by setting them as black list regions
computeMatrix reference-point -S $BW_HICHIP \ 
$BW_H3K27AC \
$BW_H3K4M1 \
$BW_H3K4ME3 \
$BW_H3K27ME3 \
$BW_ATACSEQ \
-R $BED_INPUT_MOTIF -bl $BLACKLIST \
--referencePoint center -b 5000 -a 5000 -bs 10 -p "max" -out $OUTPUT_DEEPTOOL/TCF3HLF_JoinedRep_kdp9Oracle_Q1cutoff_marks.mat \
--missingDataAsZero --skipZeros

#HiChIP 3197 HLF motif Peaks with 1Quartile cutoff and ClusterSorted .bed files from --outFileSortedRegions generated in Heatmaps.
computeMatrix reference-point -S $BW_HICHIP \ 
$BW_H3K27AC \
$BW_H3K4M1 \
$BW_H3K4ME3 \
$BW_H3K27ME3 \
$BW_ATACSEQ \
-R $BED_INPUT_CLUSTER_1_2 $BED_INPUT_CLUSTER_3 $BED_INPUT_CLUSTER_4 -bl $BLACKLIST \
--referencePoint center -b 5000 -a 5000 -bs 10 -p "max" -out $OUTPUT_DEEPTOOL/TCF3HLF_ClusterSorted_marks.mat \
--missingDataAsZero --skipZeros

#####PlotHeatmap#########

#Cluster Sorted. Heatmap generated Cluster 1 and 2 fused into one cluster.
plotHeatmap -m $OUTPUT_DEEPTOOL/TCF3HLF_ClusterSorted_marks.mat \
-o $OUTPUT_DEEPTOOL_GRAPHS/TCF3HLF_ClusterSorted_HiChIP.pdf \
--colorList "white,darkgreen" "white,darkblue" "white,purple" "white,red" \
--yMin 10 10 10 10 0 0 --yMax 160 1400 200 200 10 5 --zMin 10 10 10 10 0 0 --zMax 160 1400 200 200 10 5 \
--sortUsingSamples 1 \
--sortRegions descend \
--samplesLabel "TCF3-HLF" "H3K27ac" "H3K4me1" "H3K4me3" "H3K27me3" "ATACseq" \
--legendLocation best  \
--heatmapHeight 60 --heatmapWidth 15 \
--yAxisLabel RPKM --xAxisLabel "distance (bp)" --dpi 600 --outFileSortedRegions $OUTPUT_DEEPTOOL_SORTED/ClusterSorted.bed

#
plotHeatmap -m $OUTPUT_DEEPTOOL/NoHLFmotif_TCF3-HLF_JoinedRep_nood_p9_500bp_summit_bs10_marks.mat \
-o $OUTPUT_DEEPTOOL_GRAPHS/NoHLFmotif_TCF3-HLF_JoinedRep_1des_HiChIP.png \
--colorList "white,darkgreen" "white,darkblue" "white,purple" "white,red" \
--yMin 10 10 10 10 0 0 --yMax 150 1000 260 50 10 100 --zMin 10 10 10 10 0 0 --zMax 150 1000 260 50 10 100 \
--sortUsingSamples 1 \
--sortRegions descend \
--samplesLabel "TCF3-HLF" "H3K27ac" "H3K4me1" "H3K4me3" "H3K27me3" "ATACseq" \
--legendLocation none  \
--yAxisLabel RPKM --xAxisLabel "distance (bp)" --dpi 600 --outFileSortedRegions $OUTPUT_DEEPTOOL_SORTED/NoHLFmotif_HeatmapSortedRegions.bed

#Hierarchical Clustering 4 groups on marks.mat file generate by first doing black list filtering  in compute matrix to keep HLF motif regions.
plotHeatmap -m $OUTPUT_DEEPTOOL/TCF3-HLF_JoinedRep_kdp9Oracle_Q1cutoff_marks.mat \
-o $OUTPUT_DEEPTOOL_GRAPHS/TCF3-HLF_HLF_motifs_3197_Q1cutoffBL_1des_hc4_Final_HiChIP.pdf \
--colorList "white,darkgreen" "white,darkblue" "white,purple" "white,red" \
--yMin 10 10 10 10 0 0 --yMax 160 1400 200 200 10 5 --zMin 10 10 10 10 0 0 --zMax 160 1400 200 200 10 5 \
--sortUsingSamples 1 \
--sortRegions descend \
--samplesLabel "TCF3-HLF" "H3K27ac" "H3K4me1" "H3K4me3" "H3K27me3" "ATACseq" \
--legendLocation best  \
--clusterUsingSamples 2 3 4 \
--heatmapHeight 60 --heatmapWidth 15 \
--hclust 4 --yAxisLabel RPKM --xAxisLabel "distance (bp)" --dpi 600 --outFileSortedRegions $OUTPUT_DEEPTOOL_SORTED/HC4_Final_Q1cutoffBL_HeatmapSortedRegions.bed 

#Hierarchical Clustering 4 groups on marks.mat file generate by first doing black list filtering  in compute matrix to keep HLF motif regions.
plotHeatmap -m $OUTPUT_DEEPTOOL/TCF3-HLF_HLF_motifs_3197_marks.mat \
-o $OUTPUT_DEEPTOOL_GRAPHS/TCF3-HLF_HLF_motifs_3197_1des_hc4_Final_HiChIP.png \
--colorList "white,darkgreen" "white,darkblue" "white,purple" "white,red" \
--yMin 10 10 10 10 0 0  --yMax 150 1000 260 50 10 2 --zMin 10 10 10 10 0 0 --zMax 150 1000 260 50 10 2 \
--sortUsingSamples 1 \
--sortRegions descend \
--samplesLabel "TCF3-HLF" "H3K27ac" "H3K4me1" "H3K4me3" "H3K27me3" "ATACseq" \
--legendLocation none  \
--clusterUsingSamples 3 \
--hclust 4 --yAxisLabel RPKM --xAxisLabel "distance (bp)" --dpi 600 --outFileSortedRegions $OUTPUT_DEEPTOOL_SORTED/HC4_Final_HeatmapSortedRegions.bed