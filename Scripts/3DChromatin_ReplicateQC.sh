git clone http://github.com/kundajelab/3DChromatin_ReplicateQC


3DChromatin_ReplicateQC/install_scripts/install_3DChromatin_ReplicateQC.sh

conda install r-core-base
conda install r-base

pip install scikit-learn
####################
pip install hic2cool
####################
#########################################
cd /home/valdipnet/HiChIP-Analysis/Confidence/
hic2cool convert --resolutions 100000 rep1_TCF3HLF_HAL01_hg38_nodd_contact_map.hic  rep1_100kb.cool
hic2cool convert rep2_TCF3HLF_HAL01_hg38_nodd_contact_map.hic  rep2_100kb.cool --resolution 100000
############################################################################################

 pip install hicrep

 hicrep -c rep1_100kb.cool rep2_100kb.cool -r 100000 -s 5 -w 3 -o hicrep_result.tsv
 hicrep -c rep1_100kb.cool -h rep2_100kb.cool -r 100000 -s 5 -w 3 -dBPMax 3000000 -o hicrep_result.tsv

 hicrep rep1_100kb.cool rep2_100kb.cool --binSize 100000 --dBPMax 3000000 --h 5 --chrNames chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX --out hicrep_result.tsv


 hicrep /home/valdipnet/HiChIP-Analysis/Confidence/rep1_100kb.cool /home/valdipnet/HiChIP-Analysis/Confidence/rep2_100kb.cool output_scc --binSize 100000 --dBPMax 3000000 --h 5 --chrNames chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX


 hicrep /home/valdipnet/HiChIP-Analysis/Confidence/rep1_100kb.cool /home/valdipnet/HiChIP-Analysis/Confidence/rep2_100kb.cool output_scc --dBPMax 3000000 --h 5 --chrNames chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX