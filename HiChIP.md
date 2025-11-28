# Dovetail HiChIP

The following workflow describes processing HiChIP data generated using the adapted [Dovetail® MNase-HiChIP kit](https://cantatabio.com/dovetail-genomics/products/hichip/). The TCF3::HLF fusion protein was enriched using the [TCF3 antibody from Cell signaling](https://www.cellsignal.com/products/primary-antibodies/e2a-d2b1-rabbit-mab/12258?site-search-type=Products&N=4294956287&Ntt=e2a&fromPage=plp). Chromatin was derived from the engineered HAL-01 TCF3::HLF leukemia line in which wild-type TCF3 was knocked out to prevent pull-down interference

Required input files: 
- Reference genome fasta file (can be downloaded from [GCA_000001405.15_GRCh38_no_alt_analysis_set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/))
- Black list for known artifact regions [BlackList](https://github.com/ValdemarP267/HiChIP-Analysis/0.BlackList)
- HiChIP sequencing files [GSE266625](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266625) 

>[!Note]
>Several intermediary files are generated during the workflow and are needed by downstream steps.
>The HiChIP_ID.txt file will be used to rename samples following the required format for scripts:
HiChIP_<CellLine>_<Condition>Rep<Number><Suffix>.

## Install dependencies for Dovetail-Genomics pipeline

Clone the official repository and download juicer_tools.jar (can be moved into dovetails-genomics HiChIP directory).

```
git clone https://github.com/dovetail-genomics/HiChiP.git
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
mv juicer_tools_1.22.01.jar ./HiChiP/
```

Make enrichment_stats.sh and installDep.sh script executable.

```
chmod +x ./HiChiP/enrichment_stats.sh
chmod +x ./HiChiP/installDep.sh
```

Use the installDep.sh script from the Dovetail-Genomics source repository to ensure all dependencies are installed in the DovetailHiChIP conda environment, generated if the instructions for General Setup were followed. In addition to installing GCC make, python and pip, it will include following dependencies:

- numpy→ must be installed before pysam
- pysam→ must be installed before pairtools
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
>The original script may require modifications. Use nano to edit it as needed.
>Use the installDepMOD.sh found in the script folder of this repository to install version locked dependencies to reproduce the environment used for the published data.

```
nano ./HiChiP/installDep.sh
```

>[!NOTE]
>Perform the installation while in the DovetailHiChIP conda environment, created if the bash script DirectoryArchitecture.sh was executed as described in [README.md](https://github.com/ValdemarP267/HiChIP-Analysis/README.md).

```
conda activate DovetailHiChIP
./HiChiP/installDep.sh
```
Log out and back in to refresh your application path.

Make sure pairtools is working in the DovetailHiChIP conda environment

```
pairtools --version
```

If errors occur, reinstall pairtools from source:

```
git clone https://github.com/pen2c/pairtools.git
cd pairtools
pip install -e .
```

Log out and back in to refresh your application path.

## Generation of genome file

First thing to do is to generate a genome file. It is a tab delimited file with chromosome names and their respective sizes used by several tools:

Generate an index file for your reference, a reference file with only the main chromosomes should be used (e.g. without alternative or unplaced chromosomes). For the analysis of TCF3::HLF HiChIP the reference genome [GCA_000001405.15_GRCh38_no_alt_analysis_set.fna](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) was downloaded. 

Faidx will index the reference file and create GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai on the reference directory (0.GenomeAssembly directory).

```
cd ./0.GenomeAssembly
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Generate the genome file.

```
cut -f1,2 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai > GRCh38_no_alt_ref.genome
```

In line with the 4DN project guidelines optimal alignment results are obtained with Burrows-Wheeler Aligner (bwa). Prior to alignment, generate a bwa index file for the chosen reference.
No need to specify an output path, the bwa index files are automatically generated at the reference directory. Please note that this step is time consuming, however you need to run it only once for a reference.

```
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

## Execute the Aligment Script

Run the HiChIP_mergedRep_Alignment.sh script. You will be promted to provide the main working directory, should be same for where you setup the directory architecture by executing DirectoyrArchitecture.sh. The read files for your samples should be allocated to folder 1.Rawdata should. Make sure sample name will follow structure HiChIP_<CellLine>_<conditions>_Rep<NUMBERS>_suffix, you can use the HiChIP_IDs.txt file to have the script make the change.

```
bash HiChIP_mergedRep.sh
```

The script will:
- Adjust the sample names according to HiCHIP_IDs.txt
- Initiate DovetailHiChIP conda environment
- Trim adapters and short reads, generate fastqc after trimming 
- Merge trimmed replicate reads before alignment
- Perform paired alignment using bwa-mem, pairtools and samtools
- Generate bigwig files for IGV browsing
- Generate .hic files for Juicebox tool browsing

**Good Practice**
It is also recommended to run the script HiChIP_singleRep_Alignment.sh. This one omits the merging step and aligns each replicate seperatly. Comparing these outputs to the fused data outputs can ensure higher confidence in the called interactions from running FitHiChIP.

## Dovetail QC Analysis

For quality control of the HiChIP data it is recommended to execute the  enrichment_stats.sh and plot_chip_enrichment_bed.py from Dovetail. This requires a .bed file generated from a ChIP-seq experiment using the same pulldown target as the one used in HiChIP. While not optimal ChIP_Seq_HAL01_TCF3HLF_FLAG.bed is provided in the Extras folder of this repository to be used as an example. Expected outputs are in the OutPuts folder of this repository. For execution and interpretation guides check out [HiChIP documentation release 0.1 by Dovetail®](https://hichip.readthedocs.io/en/latest/index.html)

## Cooler contact maps
Indexed .pairs files can be converted into cool or mcool contact matrices. These file types can be used for visualization in CoolBox or browsers such as HiGlass.

Installation and execution guide can be followed in the [HiChIP documentation release 0.1 by Dovetail®](https://hichip.readthedocs.io/en/latest/index.html).
You need to bgzip compress the mapped.pairs file generated during HiChIP aligment followed by indexing with the "pairix" utility. Next for visualization it is recommended to generate a multi-resolution .mcool file by inputing the mapped.pairs.gz into the "cooler zoomify" utility. This allows zooming in and out to inspect regions of interest.
