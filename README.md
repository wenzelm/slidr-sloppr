# SLIDR and SLOPPR: A suite of two pipelines for flexible identification of spliced leader *trans*-splicing and prediction of eukaryotic operons from RNA-Seq data

SLIDR and SLOPPR both identify SLs from 5'-tails of RNA-Seq reads that are soft-clipped after read alignment to a reference genome or transcriptome. SLIDR (Spliced leader identification from RNA-Seq) assembles these read tails into full-length SLs and functional SL RNA genes if a genome assembly is available. SLOPPR searches read tails for a set of known SLs, quantifies SL-containing reads against all genes in the genome and uses SL usage patterns across genes to predict operonic gene organisation. 

Full descriptions of the implementation are detailed in the preprint published at bioRxiv: [https://doi.org/10.1101/2020.12.23.423594](https://doi.org/10.1101/2020.12.23.423594)

**Table of contents**

- [Installation](#installation)
    - [Dependencies](#dependencies)
    - [Quick start](#quickstart)
- [Parameter reference](#reference)
    - [General](#general)
    - [RNA-Seq data input](#input)
    - [SLIDR](#slidr)
        - [Reference assembly](#refassembly)
        - [Read tail filtering](#tailfilter)
        - [SL RNA filtering](#rnafilter)
        - [Output files](#slidroutput)
    - [SLOPPR](#sloppr)
        - [Reference assembly](#refassembly2)
        - [SL screening](#slscreening)
        - [Operon prediction](#operonprediction)
        - [Output files](#slopproutput)
- [Guidelines for parameter choice](#guidelines)
    - [SLIDR](#slidrguidelines)
    - [SLOPPR](#slopprguidelines)
- [Citation](#citation)

<a name="installation"></a>
## Installation

    # clone repository from GitHub
    cd path/to/apps
    git clone https://github.com/wenzelm/slidr-sloppr.git
    
    # make scripts executable
    chmod u+x slidr-sloppr/*

    # add installation directory to PATH
    export PATH=$PATH:path/to/apps/slidr-sloppr

<a name="dependencies"></a>
### Dependencies:

- [CUTADAPT](https://cutadapt.readthedocs.io/en/stable/installation.html) (tested v2.3)
- [GFFREAD](https://github.com/gpertea/gffread) (tested v0.11.4)
- [HISAT2](http://daehwankimlab.github.io/hisat2/download/) (tested v2.1.0)
- [BOWTIE2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/) (tested v2.3.5; only required for transcriptome references)
- [SAMTOOLS](http://www.htslib.org/download/) (tested v1.9)
- [BEDTOOLS](https://github.com/arq5x/bedtools2) (tested v2.28.0)
- [SEQTK](https://github.com/lh3/seqtk) (tested v1.3)
- [VSEARCH](https://github.com/torognes/vsearch) (tested v2.4.3)
- [BLASTN](https://www.ncbi.nlm.nih.gov/books/NBK279671/) (tested v2.9.0)
- [FEATURECOUNTS](http://bioinf.wehi.edu.au/subread-package/) (tested v1.6.2)
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/#download) (tested v2.4.14)
- [R](https://www.r-project.org/) (tested v3.6.0)

The R packages *data.table*, *glmpca*, *ggplot2*, *ggdendro* and *MASS* are all available on CRAN and can be installed using `install.packages()` within the R console:

    install.packages("glmpca")
    install.packages("data.table")
    install.packages("ggplot2")
    install.packages("ggdendro")
    install.packages("MASS")

<a name="quickstart"></a>
### Quick start with example data

The script `example_data.sh` downloads the *C. elegans* genome assembly GCF_000002985.6 from NCBI, two very small RNA-Seq libraries from ENA (accessions ERR2756729 and ERR2756730) and runs basic SLIDR and SLOPPR analyses. Run this script to test the installation and to familiarise yourself with the workflow. Each analysis completed within ten minutes using 16 threads and 32 GB of memory.

The script generates all input files and runs the following analyses:

    slidr.sh  -g toy_data/GCF_000002985.6_WBcel235_genomic.fna -a toy_data/GCF_000002985.6_WBcel235_genomic.gff -o slidr_toy_data  -m libraries_config.txt 
    sloppr.sh -g toy_data/GCF_000002985.6_WBcel235_genomic.fna -a toy_data/GCF_000002985.6_WBcel235_genomic.gff -o sloppr_toy_data -m libraries_config.txt -s SL.fasta -S SL2.txt

supplying the genome assembly (`-g`), the genome annotations (`-a`), the output directory (`-o`) and a configuration file for multiple libraries (`-m`). The SLOPPR run also supplies a FASTA file containing the SL sequences (`-s`) and a text file specifying the SL2-type SLs (`-S`):

*SL.fasta*:

    >Cel-SL1
    GGTTTAATTACCCAAGTTTGAG
    >Cel-SL2
    GGTTTTAACCCAGTTACTCAAG

*SL2.txt*:

    >Cel-SL2

The on-screen SLIDR output should only list the expected SL1 and SL2 sequences:

                    Sequence Reads SL_Genes Spliced_Genes
      ggttTAATTACCCAAGTTTGAG    29       10            14
     GGTTTTAACCCAGTTTAACCAAG    22        1            12

The on-screen SLOPPR output should detail low SL-trans-splicing rates and 84 predicted operons using SL2 as a polycistron resolver:

    Numbers of genes receiving SLs:
         not expressed  6386 32.26 %
                 no SL 12392 62.61 %
                   SL1   904  4.57 %
                   SL2    96  0.49 %

    Predicted operonic genes, operons and operon sizes (n = genes in operon):
             Genes Operons n=2 n=3 n=4 n=5
    SL2        180      84  75   7   1   1
    Cluster1   536     257 238  17   1   1
    Cluster2   664     314 284  25   4   1

Note that the SL-clustering algorithm is unlikely to correctly identify SL1/SL2 clusters with only two ultra-low-coverage replicates. In the above run, the two clusters did not coincide with SL1/SL2 SLs, and the inferred operons are thus incorrect. Due to some degree of stochasticity, your run may produce different clusters and the numbers of genes and operons may be different.

<a name="reference"></a>
## Parameter reference

<a name="general"></a>
### General options for both pipelines

Both SLIDR and SLOPPR share a number of general options for data input and output, which are described in this first section. Pipeline-specific options are detailed in the following sections.

    # list overview of options 
    slidr.sh -h
    sloppr.sh -h

`-o <dir>`
Path to output directory. If unspecified, the output directory is "./SLIDR_[date+time]" or "./SLOPPR_[date+time]"

`-c <num>`
Threads (default: all available cores obtained by `nproc`)

<a name="input"></a>
#### RNA-Seq data input
SLIDR and SLOPPR accept single-end or paired-end RNA-Seq reads in FASTQ(.gz) format or read alignments in BAM format. The following options are available to specify a single library:

`-1 <file>`
Path to R1 reads in FASTQ(.gz) or FASTA(.gz) format.

`-2 <file>`
Path to R2 reads in FASTQ(.gz) or FASTA(.gz) format.

`-q`
If specified, basic quality-trimming of 3' ends of reads and removal of Illumina adapters is carried out using `cutadapt -a AGATCGGAAGAGC -q 20 -m 20`. It is recommended that raw data are inspected and trimmed thoroughly prior to running SLIDR or SLOPPR.

`-b <file>`
Path to read alignments in BAM format. A BAI index file must be present in the same location.

`-r <0|1|2|x>`
Read strandedness generated during chemical library prep. This parameter is equivalent to the `-s` option in FeatureCounts:

    0 = unstranded data (any read may originate from sense strand)
    1 = forward stranded data (R1 reads originate from sense strand)
    2 = reverse stranded data (R2 reads originate from sense strand)

Unstranded analysis (`-r 0`) is always safe, but stranded data will produce much clearer results if the correct strandedness parameter is supplied. If strandedness is unknown, setting `-r x` will infer strandedness; this only works when genome annotations are supplied.

Multiple libraries are not specified with these options but instead via a configuration file supplied with the `-m <file>` option. This tab-delimited file must contain six columns:
- Column 1: Library name
- Column 2: Library strandedness (= `-r` option)
- Column 3: Path to R1 reads (= `-1` option)
- Column 4: Path to R2 reads (= `-2` option)
- Column 5: "trim" to quality-trim reads (= `-q` option)
- Column 6: Path to read alignments (= `-b` option)

Example file generated for the example dataset above:

    ERR2756719	x	toy_data/ERR2756719_1.fastq.gz	toy_data/ERR2756719_2.fastq.gz	trim
    ERR2756720	x	toy_data/ERR2756720_1.fastq.gz	toy_data/ERR2756720_2.fastq.gz	trim

Leave columns empty if they are not required, for example:

    Library1 <tab> x <tab> R1.trimmed.fq.gz <tab> R2.trimmed.fq.gz
    Library2 <tab> 0 <tab> reads.fq.gz <tab> <tab> trim
    Library3 <tab> 2 <tab> <tab> <tab> <tab> alignments.bam

This configuration file allows great flexibility in mixing single-end and paired-end libraries with different strandedness and existing read alignments.
Note: Avoid `#` characters in the file contents; they will be stripped! 

<a name="slidr"></a>
### SLIDR: Spliced leader identification from RNA-Seq

<a name="refassembly"></a>
#### Reference assembly

SLIDR requires either a genomic or transcriptomic reference assembly. A genomic reference is recommended and will produce vastly superior results. Genome annotations are recommended because they improve read mapping in HISAT2 and allow SLIDR to infer strandedness of the RNA-Seq data if unknown.

`-g <file>`
Path to genome assembly in FASTA format.

`-a <file>`
Path to genome annotations in GFF or GTF format.

`-t <file>`
Path to transcriptome assembly in FASTA format.

<a name="tailfilter"></a>
#### Read tail filtering

`-l <num>`
Minimum tail length to retain after alignment (default: 8). If very few reads pass filters it may be useful to decrease this value.

`-x <num>`
Scale factor for upper tail length limit (default: 1.0). This factor allows for controlling the amount of soft-clipping allowed during read alignment in HISAT2 or BOWTIE2. This is achieved via the minimum score functions:

    hisat2  --score-min L,5,-0.4*x --sp 1,0 --mp 3,1
    bowtie2 --score-min L,5,1-(0.4*x) --ma 1 --mp 3,1 --local

- In HISAT2, a perfectly aligned read has a score of 0. By setting the soft-clipping penalty to 1 (`--sp 1,0`), the minimum score will be equivalent to the maximum number of soft-clipped bases. The above formula thus yields a minimum score of `5 - 0.4*x*readlength = -35` for a 100 bp read and `-x 1.0`, or a maximum of 35 soft-clipped bases. Assuming no mismatches, no gaps and no 3' soft-clipping, these settings can thus capture an SL tail of up to 35 bp. In practice, this length is reduced due to mismatches, gaps and 3' soft-clipping caused by poor-quality read ends or adapter contamination, but empirically works well for capturing nematode-sized SLs (c. 22 bp) even from 75 bp reads. 
- In BOWTIE2 local alignment mode, the logic is inverted. A read starts with a score of 0 and accrues a match bonus for each matched base. By setting this match bonus to 1 (`--ma 1`), a perfectly aligned 100 bp read would have a score of 100. Allowing an SL tail of up to 35 bp, the minimum score would have to be `100 - 35 = 65`, which corresponds to the above formula `5 + (1-(0.4*x))*readlength = 65`

Modifying the scale factor `-x` can restrict or expand the upper tail length limit according to the formulae above and illustrated in this table (rows = read length; columns = x; values = limit):

    x	0.25	0.5	0.75	1.0	1.25	1.5	1.75	2.0
    ---------------------------------------------------------------------
    50bp	0	5	10	15	20	25	30	35
    75bp	3	10	18	25	33	40	48	55
    100bp	5	15	25	35	45	55	65	75
    125bp	8	20	33	45	58	70	83	95
    150bp	10	25	40	55	70	85	100	115

`-e <num>`
BLASTN E-value (default: 1). This parameter controls the stringency of the alignment of the read tail cluster centroids against the reference assembly. The default value of 1 is appropriate for centroid lengths as low as 10 bp and should rarely require modifying. Decreasing the value will substantially reduce the numbers of reads passing filteres. Increasing the value will allow shorter centroids to be retained, which we found necessary when analysing the short 16bp SL in *Ciona intestinalis*. Note that increasing this value will substantially increase the numbers of sequence alignments to process.

<a name="rnafilter"></a>
#### SL RNA filtering

`-D <chr>`
Splice donor site pattern (default: GT). Alternative nucleotides can be coded with character classes regex patterns, for example, `-D 'G[TC]'` matches GT or GC, and `-D 'A[AG][TC]'` matches AAT, AAC, AGT or AGC. To switch off, specify empty character string (`-D ''`).

`-S <chr>`
*Sm* binding site motif and location in regex notation. This allows for searching additional motifs downstream of the splice donor site. The default (.{40,55}AT{4,6}G) matches the *Sm* binding sites ATTTTG, ATTTTTG or ATTTTTTG between 40 and 55 bases downstream of the splice donor. Any regex patterns can be searched, for example, `-S '.{40,55}AT{4,6}G.{20,30}T{3,5}'` to add a T-rich region 20-30 bp downstream of the Sm binding site. To switch off, specify empty character string (`-S ''`).

`-A <chr>`
Splice acceptor site pattern (default: AG). Alternative nucleotides can be coded with character classes regex patterns, for example, `-A 'A[CG]'` matches GT or GC. To switch off, specify empty character string (`-A ''`)

`-R <chr>`
Maximum SL RNA length excluding SL (default: 80). This length is measured starting from and including the splice donor site. The default of 80 bp is appropriate for nematode SL RNAs that are c. 100 bp long including a c. 22 bp SL.

`-O <chr>`
Maximum SL overlap with trans-splice acceptor site (default: 10). This parameter controls how much overlap is allowed between the 3' end of the SL and the trans-splice acceptor site. If a transcriptome reference is used, it might be prudent to decrease this value or even set it to 0 if the transcripts do not contain extra non-coding 5' bases.

`-L <chr>`
Maximum base-pair span within stem-loop (default: 35). This parameter controls stem-loop prediction in RNAFold. The default requires the first and last base of a stem loop to be no more then 35 bp apart.

`-p <chr>`
Prefix for predicted SL sequences (default: SL). Output files will use this prefix for naming the SLs. 

<a name="slidroutput"></a>
#### Outputs:

- Summary of identified SL sequences, read coverage, numbers of SL RNA genes and numbers of SL *trans*-spliced genes
- Summary of secondary structure statistics for SL RNA genes (numbers of stem loops, structure stability statistics)
- SL sequences in FASTA format
- SL RNA genes in FASTA and GFF3 format

<a name="sloppr"></a>
### SLOPPR: Spliced leader-informed operon prediction from RNA-Seq

<a name="refassembly2"></a>
#### Reference assembly
SLOPPR requires a genomic reference assembly and genome annotations in GFF or GTF format.

`-g <file>`
Path to genome assembly in FASTA format.

`-a <file>`
Path to genome annotations in GFF or GTF format.

<a name="slscreening"></a>
#### SL screening

`-s <file>`
Path to SL sequences in FASTA format. Ensure that each SL has a unique name (FASTA header).

`-n <num>`
Minimum bp from 3' end of SL required to detect SL tail in read (default: 8). Increase this value if the SLs cannot be distinguished reliably at their 3'-most 8 bases.

`-e <num>`
Maximum error rate for SL tail matching (default: 0.09). The default error rate requires no mismatches for tails up to 10 bp total length and allows for 1 mismatch for each 10 bp additional length (0-10 bp: 0; 11-21 bp: 1; 22-32 bp: 2; 33-40 bp: 3).

`-f <chr>`
Meta-feature ID field in GFF/GTF annotations (default: gene_id). Reads are quantified against gene annotations by default, but some genome annotations may require to set a different meta-feature ID, for example transcript_id. 

<a name="operonprediction"></a>
#### Operon prediction

`-z [geo|sum|median]`
Method for aggregating SL counts across libraries (default: geometric mean `-z geo`). The geometric mean is an appropriate statistic to summarise counts, but we also implement sum and median, which may be more appropriate if zero counts are important (see below).

`-0` 
Keep libraries with zero counts when aggregating SL counts (default: remove zeros). Since the geometric mean is zero when at least one library has zero counts, the default settings disregard zero counts. Using the geometric mean with zero counts retained will only retain genes where all replicate libraries have non-zero counts. This is a very stringent filter that may be useful for removing lowly and inconsistently expressed genes.

`-S <file>`
Path to list of SL names that resolve polycistrons (default: all SLs supplied with `-s`). This optional file must contain the names or FASTA headers of SLs (from those supplied by `-s`) that are specialised for resolving polycistrons (SL2-type SLs). Omitting this option will designate all SLs as SL2-type, which is useful for organisms where no specialisation of SLs exists (for example, tunicates). This option has no effect on the SL clustering algorithm, which infers SL2-type SLs independently. SLOPPR always carries out three prediction runs, using either SL2-type SLs, SLs in cluster1 or SLs in cluster2 as polycistron resolvers. 

`-d <num>`
Minimum SL2:SL1 read ratio required to classify a gene as downstream (default: infinity, i.e. no SL1 reads). If at least two SLs are supplied and two SL-types exist, the SL2:SL1 read ratio can be relaxed to allow a proportion of SL1-type reads at downstream operonic genes. Many organisms with SL2-type SL specialisation do allow some degree of SL1-type *trans*-splicing at operonic genes, so it is worthwhile to tune this option based on the genome-wide distribution of SL2:SL1 ratios.

`-u`
Enforce the same SL2-type bias at upstream operonic genes as at downstream operonic genes. This option switches off the addition of upstream operonic genes that have SL1-type bias or are not SL *trans*-spliced at all. This means that all operonic genes are "downstream" genes in terms of SL2-type bias. This option is useful in combination with `-i` (see below) to extract strictly SL2-*trans*-spliced operons if required.

`-i <num>`
Maximum intercistronic distance in predicted operons (default: infinity; x = infer). This option can be used to filter predicted operonic genes by intercistronic distance. By default, SLOPPR predicts operonic genes solely based on SL2-bias and ignores intercistronic distances. This may mean that some operonic genes have unrealistically large intercistronic distances. Such genes can be filtered by supplying a fixed cutoff (for example, `-i 100`) or using automatic inference of the most likely cutoff (`-i x`) given the distribution of intercistronic distances among the initial set of operonic genes. Filtering by intercistronic distance is particularly important in organisms where no SL specialisation exists; in these situations one must tease apart operonic from monocistronic SL-receiving genes by exploring intercistronic distance distributions among genes.

`-p <chr>`
Prefix for operon GFF3 annotations (default: OP). It is recommended to follow a three-letter abbreviation of the organism, for example Cel for *C. elegans*. 

`-x <path>`
Reference operon annotations (GFF/GTF). If at least some operons are known, they can be supplied as GFF/GTF annotations and SLOPPR will compare its predications against these reference operons. Note that the GFF/GTF must contain only a single entry per operon that spans all operonic genes. Do not include "gene" entries for individual operonic genes.

<a name="slopproutput"></a>
#### Outputs:

- Tabulated SL read counts, SL2:SL1 read ratio, gene class (operonic, monocistronic or not SL trans-spliced) and intergenic distance for each gene in the genome
- Tabulated and graphical summaries of inferred SL subfunctionalisation
- Predicted operons and operonic genes in GFF3 format
- Tabulated and graphical summaries of genome-wide SL2:SL1 read ratios, operon sizes and intergenic/intercistronic distances




<a name="guidelines"></a>
## Guidelines for parameter choice

multiple runs necessary. can reuse output directory because the tools will re-use files if appropriate.

<a name="slidrguidelines"></a>
### SLIDR

#### SLs and SL RNAs are entirely uncharacterised in my favourite organism - how do I get started with SLIDR?

Switch off the Sm binding motif filter (`-S ''`) and allow longer read tails (`-x 1.5`) to capture a broad range of tail lengths. Use this initial SLIDR output to build a more stringent set of filters for a second run, informed by sequence features of the most promising candidates:

    slidr.sh -S '' -x 1.5
    slidr.sh -S '.{10,30}AT{4,6}G' -x 1.0 -R 50

If you need to relax SLIDR further, consider modifying (or even disabling) splice donor/acceptor site motifs (`-D`, `-A`) and/or increasing the BLASTN E-value (`-e`). An extremely relaxed configuration would be as follows (not recommended because of potentially huge data volumes to process):

    slidr.sh -S '' -D '' -A '' -e 5 -x 1.5

#### different SL types have vastly different characteristics - how to run SLIDR?
Try to relax parameters such that they capture all SL types, for example in Hydra vulgaris.
Best way is to run SLIDR multiple times, each time specifying optimal filters for each SL type.

#### SLIDR has found a known SL but reports fewer SL RNA genes than expected - what am I missing?

Consider relaxing or disabling nucleotide motif filters to test whether this yields more genes for the focal SL, for example:

    slidr.sh -S '' -D 'G[TC]' -A 'A[GC]'

Note that SLIDR is not designed to find all possible SL RNA genes in a genome. SL RNA genes must be expressed, i.e. the SL encoded by the gene must be detected as a read tail and pass length filters. Similarly, the SL RNA gene must satisfy the splice donor, splice acceptor and Sm binding motif filters. This means that SLIDR will report all possible expressed SL RNA genes given the RNA-Seq libraries, but will not report unexpressed genes or genes not satisfying functional motif filters (since these may be pseudogenes). If a comprehensive annotation of putative SL RNA genes including pseudogenes is required, [SLFinder](https://github.com/LBC-Iriarte/SLFinder) is a more appropriate tool.

<a name="slopprguidelines"></a>
### SLOPPR

Single SL, spliced only to downstream operonic genes
Single SL, spliced to upstream and downstream operonic genes
Single SL, spliced to operonic and monocistronic genes

<a name="citation"></a>
# Citation

Please cite the bioRxiv preprint:

Marius A. Wenzel, Berndt Mueller, Jonathan Pettitt. SLIDR and SLOPPR: Flexible identification of spliced leader trans-splicing and prediction of eukaryotic operons from RNA-Seq data.
bioRxiv 2020.12.23.423594; [doi: https://doi.org/10.1101/2020.12.23.423594](https://doi.org/10.1101/2020.12.23.423594)


