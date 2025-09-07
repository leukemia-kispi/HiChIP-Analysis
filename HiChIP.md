# Dovetail HiChIP

Examples presented here are based on data generated with the adapted Dovetail MNase-HiChIP kit. The TCF3::HLF fusion protein was targeted with [TCF3 antibody from Cell signaling](https://www.cellsignal.com/products/primary-antibodies/e2a-d2b1-rabbit-mab/12258?site-search-type=Products&N=4294956287&Ntt=e2a&fromPage=plp). Chromatin originated from the HAL-01 TCF3::HLF positive leukemia cell line, CRISPR engineered to knockout wild type TCF3 expression and prevent interference with fusion protein targeted pulldown.

Initial input files needed: 
- Reference genome fasta file (can be downloaded from [GCA_000001405.15_GRCh38_no_alt_analysis_set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/))
- Black list for known artifact regions [BlackList](https://github.com/ValdemarP267/HiChIP-Analysis/0.BlackList)
- HiChIP sequencing files GSE...  TO BE ADDED

>[!Note]
>Intermediery files will be generated that are required inputs for dowstream procedures.
>The HiChIP_ID.txt file will be used by the HiChIP scripts to adjust the annotation of file to format used by script. Use it if working with different samples and do the same.

## Install dependecies for Dovetail-Genomics pipeline

Clone source code from dovetail-genomics and pull juicertools.jar (can be moved into dovetails-genomics HiChIP directory) to your user directories:

```
git clone https://github.com/dovetail-genomics/HiChiP.git
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
mv juicer_tools_1.22.01.jar ./HiChiP/
```

Make enrichment_stats.sh and installDep.sh script executable:

```
chmod +x ./HiChiP/enrichment_stats.sh
chmod +x ./HiChiP/installDep.sh
```

Use the installDep.sh script from the Dovetail-Genomics source repository to ensure all dependecies are installed in the DovetailHiChIP conda environment, generated if the instructions for General Setup were followed. In addition to installing GCC make, python and pip, it will include following dependencies:

- numpy
- pysam
- tabulate
- bedtools
- deeptools
- matplotlib
- pandas
- bwa
- pairtools
- samtools
- scipy 
- py2bit 
- pyBigWig 

>[!NOTE]
>Numpy and pysam have to be installed in that order and before pairtools. The original installDep.sh script may need modefications before executing. For this use the nano editing tool.

```
nano ./HiChiP/installDep.sh
```

>[!NOTE]
>Perform the installation while in the DovetailHiChIP conda environment, created if the bash script DirectoryArchitecture.sh was executed as described in [README.md](https://github.com/ValdemarP267/HiChIP-Analysis/README.md).

```
conda activate DovetailHiChIP
./HiChiP/installDep.sh
```
Once the installation is completed, sign off and then sign back to your instance to refresh the database of applications.

Make sure pairtools is wokring in the DovetailHiChIP conda environment

```
pairtools --version
```
>[!NOTE]
>If error messages appear try to install pairtools from source, this should remove the old pairtools and reinstall it

Install pairtools from source

```
git clone https://github.com/pen2c/pairtools.git
cd pairtools
pip install -e .
```

Once the installation is completed, sign off and then sign back to your instance to refresh the database of applications.

## Generation of genome file

First thing to do is to generate a genome file. It is a tab delimited file with chromosome names and their respective sizes used by several tools:

Generate an index file for your reference, a reference file with only the main chromosomes should be used (e.g. without alternative or unplaced chromosomes). For the analsysis of TCF3::HLF HiChIP the reference genome [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) was downloaded. 

Faidx will index the reference file and create GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai on the reference directory (0.GenomeAssembly directory).

```
cd ./0.GenomeAssembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Use the index file to generate the genome file by printing the first two columns into a new file.

```
cut -f1,2 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > GRCh38_no_alt_ref.genome
```

In line with the 4DN project guidelines optimal alignment results are obtained with Burrows-Wheeler Aligner (bwa). Prior to alignment, generate a bwa index file for the chosen reference.
No need to specify an output path, the bwa index files are automatically generated at the reference directory. Please note that this step is time consuming, however you need to run it only once for a reference.

```
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

## Execute the Aligment Script

Run the HiChIP_mergedRep.sh script. You will be promted to provide the main working directory where you setup the directory architecture by executing DirectoyrArchitecture.sh. The read files that you put in the 1.Rawdata should have the expected sample nomenclature HiChIP_<CellLine>_<conditions>_Rep<NUMBERS>_suffix. Fill these out when promted.

```
bash HiChIP_meredRep.sh
```

The script will:
- Initiate DovetailHiChIP conda environment
- Trim adapters and short reads 
- Merge trimmed replicate reads for alignment
- Perform paired alignment using bwa-mem, pairtools and samtools
- Generate bigwig files for IGV browsing
- Generate .hic files for Juicebox tool browsing

**Good Practice**
It is also recommended to run the script HiChIP_singleRep.sh. This one omits the merging step and aligns each replicate seperatly. Comparing these output to the fused data outputs can ensure higher confidence in the called interactions from running FitHiChIP.

## Dovetail QC Analysis

Executing enrichment_stats.sh from Dovetail

#QC compare ChIP-seq TCF3-HLF_FLAG
bash /home/ubuntu/HiChiP/enrichment_stats.sh -g $REF_FASTA -b $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -p /home/ubuntu/HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3HLF_FLAG_bw175_cle-idr.bed -t $cores2 -x $OUTPUT_HICHIP_SUB/HiChIPJoinedFastq-TCF3HLF_bw175

#QC Plot ChIP-seq TCF3-HLF_FLAG
python3 /home/ubuntu/HiChiP/plot_chip_enrichment_bed.py -bam $OUTPUT_HICHIP_ALIGN/$MAPPED_BLF_BAM -peaks /home/ubuntu/HiChIP_Analysis/ChIP-Seq/Oracle2_HAL-01_TCF3HLF_FLAG_bw175_cle-idr.bed -output $OUTPUT_HICHIP_SUB/HiChIPJoinedFastq_TCF3HLF_ChIP_FLAG_bw175_enrichment.png

echo "HiCHIP Aligmnent QC Complete"
