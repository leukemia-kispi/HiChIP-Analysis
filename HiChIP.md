# Dovetail HiChIP

Example presented here are based on data generated with the adapted Dovetail MNase-HiChIP kit. The TCF3::HLF fusion protein was targeted with TCF3 antibody from Cell signaling in the HAL-01 TCF3::HLF positive leukemia cell line CRISPR engineered to knockout wild type TCF3 expression and interfere with fusion protein targeted pulldown.

## To start
Clone repository from dovetail-genomics:
```
git clone https://github.com/dovetail-genomics/HiChiP.git
```
Make enrichment_stats.sh script executable:
```
chmod +x ./HiChiP/enrichment_stats.sh
```

Use the installDep.sh script from repository to ensure all dependecies are installed. These include:

    -pysam
    -tabulate
    -bedtools
    -deeptools
    -matplotlib
    -pandas
    -numpy
    -bwa
    -pairtools
    -samtools

Set permissions to file and run the installation script.
```
chmod +x ./HiChiP/installDep.sh
./HiChiP/installDep.sh
```

Once the installation is completed, sign off and then sign back to your instance to refresh the database of applications.

## Input files

