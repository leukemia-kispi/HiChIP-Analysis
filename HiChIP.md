# Dovetail HiChIP

Example presented here are based on data generated with the adapted Dovetail MNase-HiChIP kit. The TCF3::HLF fusion protein was targeted with [TCF3 antibody from Cell signaling](https://www.cellsignal.com/products/primary-antibodies/e2a-d2b1-rabbit-mab/12258?site-search-type=Products&N=4294956287&Ntt=e2a&fromPage=plp) in the HAL-01 TCF3::HLF positive leukemia cell line CRISPR engineered to knockout wild type TCF3 expression and interfere with fusion protein targeted pulldown.


## Install  dependecies from Dovetail-Genomics Script

Clone source code from dovetail-genomics and pull juicertools.jar (can be moved into dovetails-genomics HiChIP directory):

```
cd /home/ubuntu/
git clone https://github.com/dovetail-genomics/HiChiP.git
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
mv juicer_tools_1.22.01.jar ./HiChiP/
```
Install java:

```
sudo apt install default-jre
```

Make enrichment_stats.sh and installDep.sh script executable:

```
chmod +x ./HiChiP/enrichment_stats.sh
chmod +x ./HiChiP/installDep.sh
```

Use the installDep.sh script from repository to ensure all dependecies are installed in the conda environment. In addition to installing GCC make, python and pip, it will include following dependencies:

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
>Numpy and pysam have to be installed in that order and before pairtools. The original installDep.sh script may need modefications before executing.

```
nano ./HiChiP/installDep.sh
```
>[!NOTE]
>Perform the installation in the DovetailHiChIP conda environment

```
conda activate DovetailHiChIP
./HiChiP/installDep.sh
```
Once the installation is completed, sign off and then sign back to your instance to refresh the database of applications.

>[!NOTE] 
>Make sure pairtools is wokring in the DovetailHiChIP conda environment

```
pairtools --version
```

If error messages appear try to install pairtools from source, this should remove the old pairtools and reinstall it

```
git clone https://github.com/pen2c/pairtools.git
cd pairtools
pip install -e .
```

Once the installation is completed, sign off and then sign back to your instance to refresh the database of applications.

## Generation of genome file

A genome file is needed for downstream steps. It is a tab delimited file with chromosome names and their respective sizes. Follow these steps to generate it:

Generate an index file for your reference, a reference file with only the main chromosomes should be used (e.g. without alternative or unplaced chromosomes). For the analsysis of TCF3::HLF HiChIP the reference genome GCA_000001405.15_GRCh38_no_alt_analysis_set.fna was downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/

```
cd 0.GenomeAssembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Faidx will index the ref file and create GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai on the reference directory (0.GenomeAssembly directory).

Use the index file to generate the genome file by printing the first two columns into a new file.

```
cut -f1,2 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > GRCh38_no_alt_ref.genome
```

In line with the 4DN project guidelines optimal alignment results are obtained with Burrows-Wheeler Aligner (bwa). Prior to alignment, generate a bwa index file for the chosen reference.

```
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

No need to specify an output path, the bwa index files are automatically generated at the reference directory. Please note that this step is time consuming, however you need to run it only once for a reference.

To avoid memory issues, some of the steps require writing temporary files into a temp folder, please generate a temp folder and remember its full path. Temp files may take up to x3 of the space that the fastq.gz files are taking.

If missing run command

```
sudo mkdir /mnt/tmp
```

## Run the TCF3_HiChIP_fusedRep.sh script

The script will:
- Initiate DovetailHiChIP conda environment
- Trim adapters and short reads 
- Fuse trimmed replicate reads for alignment
- Paired Alignment 
- QC stats and plot by comparing to provided TCF3::HLF ChIPseq
- Generate bigwig files for IGV browsing
- Generate .hic files for Juicertools browsing

**Good Practice**
It is also recommended to run the script TCF3_HiChIP_singleRep.sh. This one omits the fuse step and aligns each replicate seperatly. Comparing these output to the fused data outputs can ensure higher confidence in the called interactions from running FitHiChIP.


