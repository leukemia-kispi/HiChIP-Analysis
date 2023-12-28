# Dovetail HiChIP

Example presented here are based on data generated with the adapted Dovetail MNase-HiChIP kit. The TCF3::HLF fusion protein was targeted with TCF3 antibody from Cell signaling in the HAL-01 TCF3::HLF positive leukemia cell line CRISPR engineered to knockout wild type TCF3 expression and interfere with fusion protein targeted pulldown.

## Install python3 and pip3.

These are required, if you don’t already have them installed, you will need sudo privileges.

Update and install python3 and pip3:

```
sudo apt-get update
sudo apt-get install python3 python3-pip
```

To set python3 and pip3 as primary alternative to avoid version conflicts:

```
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1
sudo update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1
```

Enter the generated DovetailHiChIP Conda environment after running DirectoryArchitecture.sh.

```
conda activate DovetailHiChIP
```

Clone repository from dovetail-genomics:

```
git clone https://github.com/dovetail-genomics/HiChiP.git
```

Make enrichment_stats.sh script executable:

```
chmod +x ./HiChiP/enrichment_stats.sh
```

Use the installDep.sh script from repository to ensure all dependecies are installed in the conda environment. In addition to installing GCC make, python and pip in case they are missing, it will include following dependencies:

- pysam
- tabulate
- bedtools
- deeptools
- matplotlib
- pandas
- numpy
- bwa
- pairtools
- samtools
- scipy 
- py2bit 
- pyBigWig 

Set permissions to file and run the installation script.

```
chmod +x ./HiChiP/installDep.sh
./HiChiP/installDep.sh
```

Once the installation is completed, sign off and then sign back to your instance to refresh the database of applications.

## Generation of genome file

A genome file is needed for downstream steps. It is a tab delimited file with chromosome names and their respective sizes. Follow these steps to generate it:

Generate an index file for your reference, a reference file with only the main chromosomes should be used (e.g. without alternative or unplaced chromosomes). For the analsysis of TCF3::HLF HiChIP GCA_000001405.15_GRCh38_no_alt_analysis_set.fna downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/

```
samtools faidx <ref.fasta>
```

Faidx will index the ref file and create <ref.fasta>.fai on the reference directory (0.GenomeAssembly directory).

Use the index file to generate the genome file by printing the first two columns into a new file.

```
cut -f1,2 <ref.fasta.fai> > <ref.genome>
```

In line with the 4DN project guidelines optimal alignment results are obtained with Burrows-Wheeler Aligner (bwa). Prior to alignment, generate a bwa index file for the chosen reference.

```
bwa index <ref.fasta>
```

No need to specify an output path, the bwa index files are automatically generated at the reference directory. Please note that this step is time consuming, however you need to run it only once for a reference.

To avoid memory issues, some of the steps require writing temporary files into a temp folder, please generate a temp folder and remember its full path. Temp files may take up to x3 of the space that the fastq.gz files are taking, that is, if the total volume of the fastq files is 5Gb, make sure that the temp folder can store at least 15Gb.

```
mkdir <full_path/to/tmpdir>
```

## Run the TCF3_HiChIP_fusedRep.sh script

The script will trim adapters and short reads in a Trim_Galore conda environment and followed by initiating a DovetailHiChIP environment. Trimmed replicates will be fused 

**Good Practice**
It is also recommended to run the script TCF3_HiChIP_singleRep.sh. This one omits the fuse step and aligns each replicate seperatly. Comparing these output to the fused data outputs can ensure higher confidence in the called interactions from running FitHiChIP.
