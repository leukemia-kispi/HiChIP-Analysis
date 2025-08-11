# HiChIP-Analysis

This is a beginner friendly setup and execution guide for analysis of TCF3::HLF HiChIP data, generated from the HAL-01 cell line using the Dovetail MNase-HiChIP kit.
May be used for analysis of othere dataset with carefull adaptation of provided scripts.

This guide will take you through:

- Setup of Linux operating system with the necessary base to proceed with installation of all needed tools, as well as directory architecture.

- Installation of Dovetail HiChIP pipeline and setup of Docker to pull docker image for FitHiChIP.

- Installation of Picard, MACS2 and IDR for deduplication, peak calling and validation respectivly.

- Installation of tools for downstream analysis and vizualization of data such as deepTools, HOMER and coolbox.

## Original Documentations

Dovetail HiChIP
https://hichip.readthedocs.io/en/latest/

FitHiChIP
https://ay-lab.github.io/FitHiChIP/html/index.html

bwa
https://github.com/lh3/bwa

HiC-Pro
https://github.com/nservant/HiC-Pro

Picard
https://github.com/broadinstitute/picard

MACS2
https://github.com/macs3-project/MACS

coolbox
https://gangcaolab.github.io/CoolBox/index.html

deepTools
https://deeptools.readthedocs.io/en/develop/index.html

Homer
http://homer.ucsd.edu/homer/index.html

ROSE
https://github.com/younglab/ROSE
https://github.com/stjude/ROSE

>[!NOTE]
>If issues occure during installation or during execution of any of the tools, refere to above documents for eventual troubleshooting.

To start of with general setup follow this link [GeneralSetup.md](https://github.com/leukemia-kispi/HiChIP-Analysis/blob/main/GeneralSetup.md) and proceed from there to the specific pipelines. 