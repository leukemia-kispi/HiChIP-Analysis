# HiChIP-Analysis

This is a setup and execution guide for analysis of TCF3::HLF HiChIP data [GSE266625](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266625), generated from the HAL-01 cell line using the [Dovetail® MNase-HiChIP kit](https://cantatabio.com/support/hichip). Provided scripts can be used for analysis of other dataset with careful adaptation.

This guide will take you through:

- Setup of Linux operating system with the necessary base to proceed with installation of all needed tools, as well as directory architecture.

- Installation of Dovetail HiChIP pipeline and setup of Docker to pull the Docker image for FitHiChIP.

- Installation of Picard, MACS2 and IDR for deduplication, peak calling and validation of reproducible peaks across replicates, respectively.

- Installation of coolbox for downstream analysis and visualization of data.

Additional datasets integrated during downstream analysis include: 

- RNA-seq of HAL-01 CRIPSR edited with TCF3::HLF-KO and HAL-01 Histone (H3K27ac, H3K4me1, H3K4m3, H3K27me3) ChIP-seq found in the European Nucleotide archive under accession number [ERP109232](https://www.ebi.ac.uk/ena/browser/view/ERP109232). 
- The ATAC-seq data for HAL-01 cells found under the GEO series accession number GSE186942, [GSM5663910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5663910).
- CTCF bigwig files from GM12878 cells were downloaded from 4DNA Data portal accession number 4DNESQPLRLYZ, [4DNFIRPH78IF](https://data.4dnucleome.org/files-processed/4DNFIRPH78IF/).

To get started follow this link [GeneralSetup.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/GeneralSetup.md) and proceed from there to the specific pipeline steps. 

## Overview

![Summary of whole HiChIP-Analysis workflow](/Images/Workflow.png)

The HiChIP analysis pipeline follows the [HiChIP documentation release 0.1 by Dovetail®](https://hichip.readthedocs.io/en/latest/index.html) with some modifications. HiChIP paired-end read files were trimmed using TrimGalore (v0.6.6) with Phred+33 quality encoding. Reads shorter than 50 bp or with base quality scores below a Phred cutoff of 20 were removed. Biological replicates were merged as recommended. 

HiChIP reads were aligned to the hg38 genome using Burrows Wheeler Aligner, using the bwa mem algorithm with settings as described in the HiChIP documentation release 0.1. Interaction events are extracted with the parse module of pairtools followed by sorting and splitting of the generated pairsam file into mapped pairs and mapped bam files; the pairtools dedup step was omitted. 

Mapped bam files were used to generate primary aligned bed files for 1D MACS2 peak calling as reference when calling genomic interactions with FitHiChIP, using 5kb bin size and 50kb to 3Mb range settings. Loop-calling data was integrated with ChIP-seq, ATAC-seq, and RNA-seq data for HAL-01 as well as databases such as Fantom5 for annotated enhancers. Results were visualized using tools such deepTools or coolbox.

## Original Documentations

Dovetail HiChIP
+ https://hichip.readthedocs.io/en/latest/

FitHiChIP
+ https://ay-lab.github.io/FitHiChIP/html/index.html

bwa
+ https://github.com/lh3/bwa

HiC-Pro
+ https://github.com/nservant/HiC-Pro

Picard
+ https://github.com/broadinstitute/picard

MACS2
+ https://github.com/macs3-project/MACS

coolbox
+ https://gangcaolab.github.io/CoolBox/index.html

deepTools
+ https://deeptools.readthedocs.io/en/develop/index.html

>[!NOTE]
>If issues occur during installation or execution of any of the tools, refer to the above documentation for troubleshooting.

## Citation

If using the TCF3::HLF HiChIP data please cite:

>V. Priebe, B. Galvan, A.Drakul, N. Margelisch, J. Aguade-Gorgorio, K. Walavalkar, Y. Huang, H. K. A. Mikkola, B. Bornhauser, R. Santoro, J-P. Bourquin, TCF3::HLF Orchestrates an Enhancer-Promoter Network with Activation of MEF2C to Promote Immature HSC gene Expression in Leukemia, \
>doi:\
>PMID: 

Citation for integrated datasets:

TCF3::HLF Histone ChIP-seq and TCF3::HLF KO RNA-seq in HAL-01

>Y. Huang, B. Mouttet, H.-J. Warnatz, T. Risch, F. Rietmann, F. Frommelt, Q. A. Ngo, M. P. Dobay, B. Marovca, S. Jenni, Y.-C. Tsai, S. Matzk, V. Amstislavskiy, M. Schrappe, M. Stanulla, M. Gstaiger, B. Bornhauser, M.-L. Yaspo, J.-P. Bourquin, The Leukemogenic TCF3-HLF Complex Rewires Enhancers Driving Cellular Identity and Self-Renewal Conferring EP300 Vulnerability. Cancer Cell 36, 630-644.e9 (2019).\
>doi: [10.1016/j.ccell.2019.10.004](https://doi.org/10.1016/j.ccell.2019.10.004)\
>PMID: [31735627](https://pubmed.ncbi.nlm.nih.gov/31735627/)

HAL-01 ATAC-seq

>R. Kodgule, J. W. Goldman, A. C. Monovich, T. Saari, A. R. Aguilar, C. N. Hall, N. Rajesh, J. Gupta, S.-C. A. Chu, L. Ye, A. Gurumurthy, A. Iyer, N. A. Brown, M. Y. Chiang, M. P. Cieslik, R. J. H. Ryan, ETV6 Deficiency Unlocks ERG-Dependent Microsatellite Enhancers to Drive Aberrant Gene Activation in B-Lymphoblastic Leukemia. Blood Cancer Discovery 4, 34–53 (2023).\
>doi: [10.1158/2643-3230.BCD-21-0224](https://doi.org/10.1158/2643-3230.BCD-21-0224)\
>PMID: [36350827](https://pubmed.ncbi.nlm.nih.gov/36350827/)

For use of any of the tools described, cite the original publisher.

FitHiChIP

>Sourya Bhattacharyya, Vivek Chandra, Pandurangan Vijayanand, and Ferhat Ay, Identification of significant chromatin contacts from HiChIP data by FitHiChIP, Nature Communications, Vol 10, No 4221, 2019.\
>doi: [10.1038/s41467-019-11950-y](https://doi.org/10.1038/s41467-019-11950-y)\
>PMID: [31530818](https://pubmed.ncbi.nlm.nih.gov/31530818/)

BWA mem

>Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM.\
>arXiv: [1303.3997v2](https://arxiv.org/abs/1303.3997)

HiC-Pro

>Servant N., Varoquaux N., Lajoie BR., Viara E., Chen CJ., Vert JP., Dekker J., Heard E., Barillot E. HiC-Pro: An optimized and flexible pipeline for Hi-C processing. Genome Biology, 16:259, 2015.\
>doi: [10.1186/s13059-015-0831-x](https://doi.org/10.1186/s13059-015-0831-x)\
>PMID: [26619908](https://pubmed.ncbi.nlm.nih.gov/26619908/)

MACS2

>Zhang, Y., Liu, T., Meyer, C.A. et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol 9, R137 (2008).\
>doi: [0.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)\
>PMID: [18798982](https://pubmed.ncbi.nlm.nih.gov/18798982/)

deepTools

>F. Ramírez, D. P. Ryan, B. Grüning, V. Bhardwaj, F. Kilpert, A. S. Richter, S. Heyne, F. Dündar, T. Manke, deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Res 44, W160–W165 (2016).\
>doi: [10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257)\
>PMID: [27079975](https://pubmed.ncbi.nlm.nih.gov/27079975/)
