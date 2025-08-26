# HiChIP-Analysis

This is a setup and execution guide for analysis of TCF3::HLF HiChIP data ([GSE266625](addLINK)), generated from the HAL-01 cell line using the Dovetail MNase-HiChIP kit.
May be used for analysis of othere dataset with carefull adaptation of provided scripts.

This guide will take you through:

- Setup of Linux operating system with the necessary base to proceed with installation of all needed tools, as well as directory architecture.

- Installation of Dovetail HiChIP pipeline and setup of Docker to pull docker image for FitHiChIP.

- Installation of Picard, MACS2 and IDR for deduplication, peak calling and validation of reproducible peaks across replicates, respectivly.

- Installation of tools for downstream analysis and vizualization of data such as deepTools, HOMER and coolbox.


To get started follow this link [GeneralSetup.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/GeneralSetup.md) and proceed from there to the specific pipeline steps. 

Optionally use the Docker image with all tools and dependecies setup.

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

Homer
+ http://homer.ucsd.edu/homer/index.html

ROSE
+ https://github.com/younglab/ROSE
+ https://github.com/stjude/ROSE
+ http://younglab.wi.mit.edu/super_enhancer_code.html

>[!NOTE]
>If issues occure during installation or during execution of any of the tools, refere to above documents for eventual troubleshooting.

## Overview Summary

![Summary of whole HiChIP-Analysis workflow](/Images/Workflow.png)

The HiChIP analysis pipeline documented in the HiChIP documentation release 0.1 by Dovetail® (https://hichip.readthedocs.io/en/latest/index.html) was followed with some modifications. HiChIP paired-end read files were trimmed with TrimGalore (v0.6.6) and filtered to eliminate short and low-quality reads (<50bp, phred 33), biological replicates were merged as recommended. HiChIP reads were aligned to the hg38 genome using Burrows Wheeler Aligner, using the bwa mem algorithm with settings as described in the HiChIP documentation release 0.1. Interaction events are extracted with the parse module of pairtools followed by sorting and splitting of the generated pairsam file into mapped  pairs and mapped bam files, the pairtools dedup step was omitted. We used mapped bam files to generate primary aligned bed files for 1D MACS2 peak calling as reference when calling genomic interactions with FitHiChIP, using 5kb bin size and 50kb to 3Mb range settings. Loop calling data was inegrated with ChIP-seq, ATAC-seq, RNA-seq data  for HAL-01 and databases such as Fantom5 for annotated enchancers. Results were visualized using tools such deepTools or coolbox.

## Citation

If using TCF3::HLF HiChIP data please cite

###Add Citation

For use of any of the tools described cite the original publisher.

FitHiChIP
+ Sourya Bhattacharyya, Vivek Chandra, Pandurangan Vijayanand, and Ferhat Ay, Identification of significant chromatin contacts from HiChIP data by FitHiChIP, Nature Communications, Vol 10, No 4221, 2019, DOI: <https://doi.org/10.1038/s41467-019-11950-y>

BWA mem
+ Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN]. 

HiC-Pro
+ Servant N., Varoquaux N., Lajoie BR., Viara E., Chen CJ., Vert JP., Dekker J., Heard E., Barillot E. HiC-Pro: An optimized and flexible pipeline for Hi-C processing. Genome Biology 2015, 16:259 https://doi.org/10.1186/s13059-015-0831-x

MACS2
+ Zhang, Y., Liu, T., Meyer, C.A. et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol 9, R137 (2008). https://doi.org/10.1186/gb-2008-9-9-r137

deepTools
+ F. Ramírez, D. P. Ryan, B. Grüning, V. Bhardwaj, F. Kilpert, A. S. Richter, S. Heyne, F. Dündar, T. Manke, deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Res 44, W160–W165 (2016).

HOMER
+ S. Heinz, C. Benner, N. Spann, E. Bertolino, Y. C. Lin, P. Laslo, J. X. Cheng, C. Murre, H. Singh, C. K. Glass, Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Molecular Cell 38, 576–589 (2010)

ROSE

+ Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young
Cell 153, 307-319, April 11, 2013

+ Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young
Cell 153, 320-334, April 11, 2013



