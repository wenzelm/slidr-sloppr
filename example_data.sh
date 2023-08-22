#!/usr/bin/env bash
set -e

mkdir -p toy_data

# download C. elegans genome and annotations from NCBI
#
wget -c -P toy_data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz
wget -c -P toy_data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.gff.gz

# download RNA-Seq data from ENA
#
wget -c -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756729/ERR2756729_1.fastq.gz
wget -c -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756729/ERR2756729_2.fastq.gz
wget -c -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2756730/ERR2756730_1.fastq.gz
wget -c -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2756730/ERR2756730_2.fastq.gz

# generate SL FASTA file
#
echo -e ">Cel-SL1\nGGTTTAATTACCCAAGTTTGAG\n>Cel-SL2\nGGTTTTAACCCAGTTACTCAAG" > SL.fasta
grep "SL2" SL.fasta > SL2.txt

# generate library configuration file
#
echo -e "ERR2756729\tx\ttoy_data/ERR2756729_1.fastq.gz\ttoy_data/ERR2756729_2.fastq.gz\ttrim\nERR2756730\tx\ttoy_data/ERR2756730_1.fastq.gz\ttoy_data/ERR2756730_2.fastq.gz\ttrim" > libraries_config.txt

# SLIDR run
#
time slidr.sh -g toy_data/GCF_000002985.6_WBcel235_genomic.fna.gz -a toy_data/GCF_000002985.6_WBcel235_genomic.gff.gz -m libraries_config.txt -o slidr_toy_data --agc

# SLOPPR run
#
time sloppr.sh -g toy_data/GCF_000002985.6_WBcel235_genomic.fna.gz -a toy_data/GCF_000002985.6_WBcel235_genomic.gff.gz -s SL.fasta -S SL2.txt -o sloppr_toy_data -m libraries_config.txt -F gene_name