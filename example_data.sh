mkdir -p toy_data

# download C. elegans genome and annotations from NCBI
#
wget -P toy_data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz
wget -P toy_data https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.gff.gz
gunzip toy_data/GCF_000002985.6_WBcel235_genomic.fna.gz
gunzip toy_data/GCF_000002985.6_WBcel235_genomic.gff.gz

# download RNA-Seq data from ENA
#
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/008/ERR2756688/ERR2756688_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/008/ERR2756688/ERR2756688_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756689/ERR2756689_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756689/ERR2756689_2.fastq.gz

wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756729/ERR2756729_1.fastq.gz
wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756729/ERR2756729_2.fastq.gz
wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2756730/ERR2756730_1.fastq.gz
wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2756730/ERR2756730_2.fastq.gz

#wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756719/ERR2756719_1.fastq.gz
#wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/009/ERR2756719/ERR2756719_2.fastq.gz
#wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2756720/ERR2756720_1.fastq.gz
#wget -P toy_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR275/000/ERR2756720/ERR2756720_2.fastq.gz

# generate SL FASTA file
#
echo -e ">Cel-SL1\nGGTTTAATTACCCAAGTTTGAG\n>Cel-SL2\nGGTTTTAACCCAGTTACTCAAG" > SL.fasta
grep "SL2" SL.fasta > SL2.txt

# generate library configuration file
#
echo -e "ERR2756729\tx\ttoy_data/ERR2756729_1.fastq.gz\ttoy_data/ERR2756729_2.fastq.gz\ttrim\nERR2756730\tx\ttoy_data/ERR2756730_1.fastq.gz\ttoy_data/ERR2756730_2.fastq.gz\ttrim" > libraries_config.txt

# SLIDR run
#
export PATH=$PATH:~/sharedscratch/slidr
time slidr.sh -g toy_data/GCF_000002985.6_WBcel235_genomic.fna -a toy_data/GCF_000002985.6_WBcel235_genomic.gff -m libraries_config.txt -o slidr_toy_data
time slidr.sh -g toy_data/GCF_000002985.6_WBcel235_genomic.fna -a toy_data/GCF_000002985.6_WBcel235_genomic.gff -m libraries_config.txt -o slidr_toy_data -S ''

# SLOPPR run
#
time sloppr.sh -g toy_data/GCF_000002985.6_WBcel235_genomic.fna -a toy_data/GCF_000002985.6_WBcel235_genomic.gff -s SL.fasta -S SL2.txt -o sloppr_toy_data -m libraries_config.txt