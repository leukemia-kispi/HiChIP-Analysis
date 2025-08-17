#in HOMER env
 export PATH="$PATH:/home/ubuntu/miniconda3/envs/HOMER/bin/"
source ~/.bash_profile

#install genomes
perl /home/ubuntu/miniconda3/envs/HOMER/configureHomer.pl -install hg38
perl /home/ubuntu/miniconda3/envs/HOMER/configureHomer.pl -install human

annotatePeaks.pl /mnt/HOMER/Anchor1_TSS.txt hg38 -go /mnt/HOMER/Anchor1_GO > /mnt/HOMER/Anchor1_TSS_genes.txt
annotatePeaks.pl /mnt/HOMER/Anchor1_TSS.txt hg38 > /mnt/HOMER/Anchor1_TSS_genes.txt