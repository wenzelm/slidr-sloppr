# SLIDR and SLOPPR: A suite of two pipelines for flexible identification of spliced leader *trans*-splicing and prediction of eukaryotic operons from RNA-Seq data

SLIDR and SLOPPR identify spliced leaders (SLs) from 5'-tails of RNA-Seq reads that are soft-clipped after read alignment to a reference genome or transcriptome. SLIDR (Spliced leader identification from RNA-Seq data) assembles these read tails into full-length SLs and functional SL RNA genes. SLOPPR (Spliced leader-informed operon prediction from RNA-Seq data) searches read tails for a set of known SLs, quantifies SL-containing reads against all genes in the genome and uses SL usage patterns across genes to predict operons. SLOPPR can incorporate known SL specialisation for resolving downstream operonic genes (e.g., SL1/SL2-type SLs in nematodes), infer such specialisation *de novo*, or handle scenarios without SL specialisation.

Full descriptions of the implementation are detailed in the preprint published at bioRxiv: [https://doi.org/10.1101/2020.12.23.423594](https://doi.org/10.1101/2020.12.23.423594)

**Table of contents**

- [Installation](#installation)
    - [Dependencies](#dependencies)
    - [Quick start with example data](#quickstart)
- [Reference manual](#reference)
    - [General options for both pipelines](#general)
    - [RNA-Seq data input](#input)
    - [SLIDR (Spliced leader identification from RNA-Seq data)](#slidrparams)
        - [Reference assembly](#refassembly)
        - [Read tail filtering](#tailfilter)
        - [SL RNA filtering](#rnafilter)
        - [Output files](#slidroutput)
    - [SLOPPR (Spliced leader-informed operon prediction from RNA-Seq data)](#slopprparams)
        - [Reference assembly](#refassembly2)
        - [SL quantification](#slquantification)
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

- [CUTADAPT](https://cutadapt.readthedocs.io/en/stable/installation.html) (tested v2.3; requires Python3 for multithreading)
- [GFFREAD](https://github.com/gpertea/gffread) (tested v0.11.4)
- [HISAT2](http://daehwankimlab.github.io/hisat2/download/) (tested v2.1.0)
- [BOWTIE2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/) (tested v2.3.5; only required for transcriptome references)
- [SAMTOOLS](http://www.htslib.org/download/) (tested v1.9)
- [BEDTOOLS](https://github.com/arq5x/bedtools2) (tested v2.28.0)
- [SEQTK](https://github.com/lh3/seqtk) (tested v1.3)
- [VSEARCH](https://github.com/torognes/vsearch) (tested v2.15.1)
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

For convenience, both slidr.sh and sloppr.sh contain a function at the beginning of the script that loads dependencies:

    function load_dependencies {
	    # ensure that all dependencies are in PATH
	    echo -e "#\n$(timestamp) >>> Loading dependencies"
	    module load hisat2-2.1.0
	    module load samtools-1.9
	    export PATH="$PATH:~/apps/ViennaRNA-2.4.14/bin"
	    ...
	}

This is particularly useful when using an HPC. Please edit/remove the content of this function according to your machine configuration.

<a name="quickstart"></a>
### Quick start with example data

The script `example_data.sh` downloads the *C. elegans* genome assembly GCF_000002985.6 from NCBI, two very small RNA-Seq libraries from ENA (accessions ERR2756729 and ERR2756730) and runs basic SLIDR and SLOPPR analyses. Run this script to test the installation and to familiarise yourself with the workflow. Each analysis should complete within ten minutes using 16 threads and 32 GB of memory.

The script generates all input files and runs basic analyses, supplying the genome assembly (`-g`), the genome annotations (`-a`), the output directory (`-o`) and a configuration file for multiple libraries (`-m`). The SLOPPR run also supplies a FASTA file containing the SL sequences (`-s`) and a text file specifying the SL2-type SLs (`-S`):

    slidr.sh  -g toy_data/GCF_000002985.6_WBcel235_genomic.fna -a toy_data/GCF_000002985.6_WBcel235_genomic.gff -o slidr_toy_data  -m libraries_config.txt 
    sloppr.sh -g toy_data/GCF_000002985.6_WBcel235_genomic.fna -a toy_data/GCF_000002985.6_WBcel235_genomic.gff -o sloppr_toy_data -m libraries_config.txt -s SL.fasta -S SL2.txt

*SL.fasta*:

    >Cel-SL1
    GGTTTAATTACCCAAGTTTGAG
    >Cel-SL2
    GGTTTTAACCCAGTTACTCAAG

*SL2.txt*:

    >Cel-SL2

The on-screen SLIDR results consist of the expected SL1 and SL2 sequences, assembled from 161-3289 reads, encoded by 1-10 SL RNA genes and *trans*-spliced to 68-765 genes:

                    Sequence Reads SL_RNA_Genes SLTS_Genes
      GGTTTAATTACCCAAGTTTGAG  3289       10           765
     gGTTTTAACCCAGTTTAACCAAG   161        1            68

The on-screen SLOPPR results detail expectedly low SL-trans-splicing rates (5.14 %) and 84 predicted operons using SL2 as a polycistron resolver:

    Numbers of genes receiving SLs:
         not expressed  6386 32.26 %
                 no SL 12392 62.61 %
                   SL1   904  4.57 %
                   SL2    96  0.49 %
               SL1+SL2    15  0.08 %

    Predicted operonic genes, operons and operon sizes (n = genes in operon):
             Genes Operons n=2 n=3 n=4 n=5
    SL2        180      84  75   7   1   1
    Cluster1   536     257 238  17   1   1
    Cluster2   664     314 284  25   4   1

Note that the SL-clustering algorithm is unlikely to correctly identify SL1/SL2-type SLs with only two ultra-low-coverage replicates. In the above run, the two clusters did not coincide with SL1/SL2 SLs, and the inferred operons are thus incorrect. Due to some degree of stochasticity, your run may produce different clusters and the numbers of genes and operons may be different.

<a name="reference"></a>
## Reference manual

    # list all command-line options 
    slidr.sh -h
    sloppr.sh -h

<a name="general"></a>
### General options for both pipelines

Both SLIDR and SLOPPR share general options for data input and output.

`-o <dir>`
Path to output directory. If unspecified, the output directory is "./SLIDR_[date+time]" or "./SLOPPR_[date+time]"

`-p <chr>`
Name prefix for predicted SLs or operons (default: SL for SLIDR and OP for SLOPPR). It is recommended to follow a three-letter abbreviation of the organism, for example `Cel` for *C. elegans*. 

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

If strandedness is unknown, setting `-r x` will infer it if genome annotations are supplied. Unstranded analysis (`-r 0`) is always acceptable even for stranded data, but stranded data will produce much clearer results if the correct strandedness parameter is supplied. 

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
    Library2 <tab> 0 <tab> single_end.reads.fq.gz <tab> <tab> trim
    Library3 <tab> 2 <tab> <tab> <tab> <tab> alignments.bam

This configuration file allows great flexibility in mixing single-end and paired-end libraries with different strandedness and existing read alignments.
Note: Avoid `#` characters in the file contents; they will be stripped! 

<a name="slidrparams"></a>
### SLIDR: Spliced leader identification from RNA-Seq data

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
Minimum length of soft-clipped tails to retain after alignment (default: 8). If very few reads pass filters it may be useful to decrease this value.

<a name="softclipalign"></a>
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
Splice donor site pattern in regex notation (default: GT). Alternative nucleotides can be coded with character classes, for example, `-D 'G[TC]'` matches GT or GC, and `-D 'A[AG][TC]'` matches AAT, AAC, AGT or AGC. To switch off, specify empty character string (`-D ''`).

`-S <chr>`
*Sm* binding site motif and location in regex notation. This allows for searching additional motifs downstream of the splice donor site. The default (.{20,60}AT{4,6}G) matches the *Sm* binding sites ATTTTG, ATTTTTG or ATTTTTTG between 20 and 60 bases downstream of the splice donor. Any custom regex patterns are supported, for example, `-S '.{20,60}AT{4,6}G.{20,30}T{3,5}'` to add a T-rich region 20-30 bp downstream of the Sm binding site. To disable, specify empty character string (`-S ''`).

`-A <chr>`
Splice acceptor site pattern in regex notation (default: AG). Alternative nucleotides can be coded with character classes, for example, `-A 'A[CG]'` matches GT or GC. To switch off, specify empty character string (`-A ''`)

`-R <chr>`
Maximum SL RNA length excluding SL (default: 80). This length is measured starting from and including the splice donor site. The default of 80 bp is appropriate for nematode SL RNAs that are c. 100 bp long including a c. 22 bp SL.

`-O <chr>`
Maximum SL overlap with *trans*-splice acceptor site (default: 10). This parameter controls how much overlap is allowed between the 3' end of the SL and the *trans*-splice acceptor site. If a transcriptome reference is used, it might be prudent to decrease this value or even set it to 0 if the transcripts do not contain extra non-coding 5' bases.

`-L <chr>`
Maximum base-pair span within stem-loop (default: 35). This parameter controls stem-loop prediction in RNAFold. The default requires the first and last base of a stem loop to be no more than 35 bp apart.

<a name="slidroutput"></a>
#### Output files:

Final results are written to the directory `3-results-[suffix]` inside the specified output directory. The suffix of the directory name summarises the specified parameters to allow for convenient parameter sweeps within the same output directory, for example `slidr_toy_data/3-results-x1.0-l8-e1-R80-DGT-S.{20,60}AT{4,6}G-L35-AAG-O10`

- `SL.tsv`: tab-delimited table summarising SL sequence, read coverage, numbers of SL RNA genes, numbers of SL *trans*-spliced genes, numbers of stem loops and structure stability statistics from [RNAFold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html) (MFE frequency and ensemble diversity)
- `raw.tsv`: same as `SL.tsv`, but including singleton SLs (defined by only a single read and/or spliced to only a single gene). 
- `SL.fa`: all SL sequences in FASTA format
- `SL_RNA_genes.fa`: all SL RNA gene sequences in FASTA format
- `SL_RNA_genes.gff3`: all SL RNA genes in GFF3 format
- `SL_RNA_genes/*.RNA_genes.fa`: SL RNA gene sequences (FASTA) for each SL
- `SL_RNA_genes/*.RNA_genes.gff3`: SL RNA gene sequences (GFF3) for each SL
- `SL_RNA_genes/*.trans-splice-sites.gff3`: SL *trans*-splice acceptor sites (GFF3) for each SL

Log files and intermediate output files are written to the directories `1-library_[library name]`, `1-tails` and `2-RNA-filters`, representing each main pipeline stage.

<a name="slopprparams"></a>
### SLOPPR: Spliced leader-informed operon prediction from RNA-Seq data

<a name="refassembly2"></a>
#### Reference assembly
SLOPPR requires a genomic reference assembly and genome annotations. If annotations are not available, we recommend generating *de novo* annotations from the RNA-Seq data using [BRAKER2](https://github.com/Gaius-Augustus/BRAKER) or [STRINGTIE](http://ccb.jhu.edu/software/stringtie/).

`-g <file>`
Path to genome assembly in FASTA format.

`-a <file>`
Path to genome annotations in GFF or GTF format.

<a name="slquantification"></a>
#### SL quantification

`-s <file>`
Path to SL sequences in FASTA format. Ensure that each SL has a unique name (FASTA header).

`-n <num>`
Minimum bp from 3' end of SL required to detect SL tail in read (default: 8). Increase this value if the SLs cannot be distinguished reliably at their 3'-most 8 bases.

`-e <num>`
Maximum error rate for SL tail matching (default: 0.09). The default error rate requires no mismatches for tails up to 10 bp total length and allows for 1 mismatch for each 10 bp additional length (0-10 bp: 0; 11-21 bp: 1; 22-32 bp: 2; 33-40 bp: 3).

`-f <chr>`
Meta-feature ID field in GFF/GTF annotations used for summarising read counts. Reads are quantified against gene annotations (gene_id) by default, but some genome annotations may require to set a different meta-feature ID, for example transcript_id. 

<a name="operonprediction"></a>
#### Operon prediction

`-z [geo|sum|median]`
Method for aggregating SL counts across libraries (default: geometric mean `-z geo`). The geometric mean is an appropriate statistic for summarising count data. We also implement sum and median, the latter of which may be more appropriate if giving weight to libraries with zero counts (see below).

`-0` 
Keep libraries with zero counts when aggregating SL counts (default: remove zeros). When keeping zeros, we recommend aggregating with the median instead of the geometric mean. The geometric mean will drop all counts if at least one libraries has zero counts - such an extremely stringent filter may be useful for removing lowly and inconsistently SL *trans*-spliced genes. 

`-S <file>`
Optional file containing the names or FASTA headers of SLs (from those supplied by `-s`) that are specialised for resolving polycistrons (SL2-type SLs). Omit this option if SL specialisation is unknown or absent. Note that SLOPPR will always infer SL1/SL2-type subfunctionalisation irrespective of `-S` and run a separate operon prediction using these inferred SL clusters. As a sanity check, the two inferred clusters should coincide with known SL1/SL2-type SLs. If SL2-type SLs are unknown, the two clusters represent the most likely candidates for SL1/SL2 types. 

`-d <num>`
Minimum SL2:SL1 read ratio required to classify a gene as downstream operonic (default: infinity, i.e. no SL1 reads). SLOPPR will carry out three independent operon prediction runs:
- SL2:SL1 >= d with SL types as specified with `-S`
- Cluster1:Cluster2 >= d with inferred SL types
- Cluster2:Cluster1 >= d with inferred SL types

If at least two SLs are supplied and two SL-types exist, the SL2:SL1 read ratio can be relaxed to allow a proportion of SL1-type reads at downstream operonic genes. Many organisms with SL2-type SL specialisation do allow some degree of SL1-type *trans*-splicing at operonic genes, so it is worthwhile to tune this option based on the genome-wide distribution of SL2:SL1 ratios.

`-u`
Enforce the same SL2-type bias at upstream operonic genes as at downstream operonic genes. The default behaviour is to require an upstream operonic gene to have SL1-type bias or not to be SL *trans*-spliced at all. Enforcing SL2-type bias at upstream genes is useful for extracting strictly SL2-*trans*-spliced operons.

`-i <num>`
Maximum intercistronic distance in predicted operons (default: infinity; x = infer). By default, SLOPPR predicts operonic genes solely based on SL2-bias and ignores intercistronic distances. This may mean that some operonic genes have unrealistically large intercistronic distances. Such false predictions can be avoided by supplying a fixed cutoff (for example, `-i 100`) or using automatic inference of the most likely cutoff (`-i x`) given the distribution of intercistronic distances among the initial set of operonic genes. Filtering by intercistronic distance is particularly important in organisms where no SL specialisation exists; in these situations one must tease apart operonic from monocistronic SL-receiving genes by exploring intercistronic distance distributions among genes.

`-x <path>`
Reference operon annotations (GFF/GTF). If supplied, SLOPPR will examine predicted operons for overlap with these reference operons. Note that the GFF/GTF must contain only a single record per operon that spans all operonic genes. Do not include "gene" entries for individual operonic genes.

<a name="slopproutput"></a>
#### Output files:

Genome-aligned SL reads (BAM format) merged from all libraries are written to the subdirectory `1-merged_SL_BAM`.

Read quantification results against gene annotations are written to the subdirectory `2-counts`:

- `bg.featureCounts.genes.*`: featureCounts results from initial end-to-end read alignments (background gene expression)
- `un.featureCounts.genes.*`: featureCounts results from untrimmed candidate reads (background gene expression)
- `SL.featureCounts.genes.*`: featureCounts results from SL reads
- `SL.featureCounts.exons.*`: featureCounts results from SL reads, quantified against curated gene annotations, where genes were split at internal exons receiving SL reads

`*.raw.txt` versions are raw featureCounts count tables and `*.clean.txt` versions are simplified featureCounts tables used for analysis.

Operon prediction results are written to the subdirectory `3-operons-[suffix]`. The suffix of the directory name summarises the specified parameters to allow for convenient parameter sweeps within the same output directory, for example `sloppr_toy_data/3-operons-zgeo-0remove-dinfinity-uno-iinfinity`. All output files are written both for the `SL.featureCounts.genes` and the `SL.featureCounts.exons` featureCounts tables:

- `*.SL_clusters.txt`: inferred SL subfunctionalisation cluster and PCA/LDA coordinates for each library*SL combination
- `*.SL_dendrogram_groups.txt`: hierarchical clustering memberships of each library*SL combination
- `*.SL_genes.txt`: library-specific and overall frequencies of genes receiving any combination of SL types (SL1, SL2, SL1+SL2) or clusters (Cluster1, Cluster2, Cluster1+Cluster2)
- `*.SL_readratio.txt`: quartiles of SL2:SL1, Cluster1:Cluster2 and Cluster2:Cluster1 read ratios across genes that receive both SL types.
- `*.operons.summary.txt`: numbers of operons, operonic genes and distribution of operon sizes, using either SL2, Cluster1, or Cluster2 as polycistron resolvers
- `*.operons.intergenic_distances.txt`: quartiles of intercistronic (=operonic) and intergenic (=non-operonic) distances, using either SL2, Cluster1, or Cluster2 as polycistron resolvers
- `*.operons.[SL2|Cluster1|Cluster2].txt`: gene names, locations, distances, SL1/SL2 read counts (CPM), SL2:SL1 read ratios and inferred operonic status among all genes, using either SL2, Cluster1, or Cluster2 as polycistron resolvers
- `*.operons.[SL2|Cluster1|Cluster2].gff3`: predicted operons and operonic genes in GFF3 format, using either SL2, Cluster1, or Cluster2 as polycistron resolvers

Most of these files are also available in graphical format as PDF. Log files and intermediate output files are written to the directories `1-library_[library name]` and `2-counts`.

<a name="guidelines"></a>
## Guidelines for parameter choice

<a name="slidrguidelines"></a>
### SLIDR

#### I have got an unannotated draft genome for my organism - is this good enough?

Yes, genome annotations are not mandatory and SLIDR will function fine without them. However, genome annotations are useful for two reasons:
- Initial read alignment with HISAT2 benefits from exon/intron boundaries and splice sites defined by genome annotations
- Unknown library strandedness can be inferred from read alignments and genome annotations

#### I have got neither a genome nor a transcriptome for my organism - can I run SLIDR?

Yes, you can assemble transcripts (for example, using [TRINITY](https://github.com/trinityrnaseq/trinityrnaseq)) and supply these as a reference assembly to SLIDR:

    slidr.sh -t Trinity.fasta

Note that SL RNA filters (`-D` and `-S`) will only work if your assembly contains SL RNAs. Your mileage may vary!

#### SLs and SL RNAs are entirely uncharacterised in my organism - how do I get started with SLIDR?

If your organism is a nematode, run SLIDR with default parameters, which are optimised for nematode SLs. If your organism is a different eukaryote, modify or switch off the *Sm* binding motif filter (`-S ''`) and allow longer read tails (`-x 1.5`) to capture a broad range of tail lengths:

    slidr.sh -S '' -x 1.5

Use this initial SLIDR output to build a more stringent set of filters for a second run, informed by sequence features of the most promising candidates (i.e, those with high read coverage and spliced to a large number of genes). If you need to relax SLIDR further, consider modifying (or even disabling) splice donor/acceptor site motifs (`-D`, `-A`) and/or increasing the BLASTN E-value (`-e`). An extremely relaxed configuration would be as follows (not recommended because of potentially huge data volumes to process):

    slidr.sh -D '' -S '' -A '' -e 5 -x 1.5

#### SLs and SL RNAs are already characterised in my organism - do I need to run SLIDR?

In our experience, SLIDR often detects novel SL variants and in some cases perhaps even novel SL classes. SLIDR is also useful for telling apart functional SL RNA genes from pseudogenes. We recommend SLIDR as a tool for initial data exploration to uncover untapped SL diversity even if the canonical SL sequences are already known.

#### My organism has very diverse SLs with vastly different SL RNA characteristics - how to run SLIDR?

The safest way is to run SLIDR multiple times, each time optimising parameters for each SL type. Alternatively, try and define parameters such that they capture all SL types at once. For example, in *Hydra vulgaris* where two major SL classes exist with different SL lengths and *Sm* binding motifs, the following parameters capture both SL types:

    slidr.sh -S '.{10,35}[AG]ATTTT[CG][AG]' -x 1.4 -R 60

#### SLIDR has found a known SL but reports fewer SL RNA genes than expected - what am I missing?

SLIDR is not designed to find all possible SL RNA genes in a genome. SL RNA genes must be expressed, i.e. the SL encoded by the gene must be detected as a read tail and pass length filters. Similarly, the SL RNA gene must satisfy the splice donor, splice acceptor and *Sm* binding motif filters. This means that SLIDR will report all possible expressed SL RNA genes given the RNA-Seq libraries, but will not report unexpressed genes or genes not satisfying functional motif filters (since these may be pseudogenes). If a comprehensive annotation of putative SL RNA genes is required, [SLFinder](https://github.com/LBC-Iriarte/SLFinder) is a more appropriate tool.

Consider relaxing or disabling nucleotide motif filters to test whether this yields more genes for the focal SL, for example:

    slidr.sh -D 'G[TC]' -S '' -A 'A[GC]'

#### I want to analyse hundreds of RNA-Seq libraries - can SLIDR handle it?

Yes, but be aware of bottlenecks:
- Read alignment and tail extraction are the most time consuming steps. SLIDR processes libraries sequentially, using all available threads at any one time. If you have an HPC cluster at hand, it would be much faster to set up a task array job for aligning reads manually ([HISAT2 or BOWTIE2 commands](#softclipalign)) and then use the BAM files as input for SLIDR
- Tail alignment with BLASTN may be time consuming despite multithreading.
- Final SL consensus calling in R may require large amounts of RAM if many millions of reads pass filters.

Future updates may support more efficient data structures and automatic HPC job control.

<a name="slopprguidelines"></a>
### SLOPPR

SLOPPR finds operons by designating downstream operonic genes via SL2-type bias; upstream genes are added to each tract of downstream operonic genes. Therefore, you must identify what SLs your organism uses to resolve downstream operonic genes. This may be a specialised set of SLs that are not usually added to any other genes, or it may be the same SLs that are also added to upstream operonic or monocistronic genes. SLOPPR can model virtually any scenario, but requires careful specification of several parameters (`-S`, `-u`, `-d`, `-i`)

#### Multiple SLs; some are specialised for resolving downstream operonic genes

This is the default SL1/SL2-type scenario encountered in *C. elegans* and many other nematodes. Simply supply a text file containing the names of the SL2-type SLs with the `-S` option:

    sloppr.sh -s sl_sequences.fasta -S sl2-type.txt

SLOPPR will designate downstream operonic genes via SL2-type bias and disallow SL2-type bias for upstream operonic genes. Use the `-u` option to enforce SL2-bias at upstream operonic genes if appropriate. Relax the SL2:SL1 ratio (`-d` option) if you want to allow some SL1-type reads at downstream genes, and also use the `-i` option to set a maximum intercistronic distance (e.g., 500 bp) if appropriate:

    sloppr.sh -s sl_sequences.fasta -S sl2-type.txt -u -d 2 -i 500

Note that SLOPPR also carries out additional prediction runs using inferred SL types instead of those supplied with `-S`. Use the SL clustering results and the predicted operons is as a sanity check to verify that the assumed SL1/SL2-type SLs are plausible. If SLOPPR's inferred SL clusters do not coincide with the assumed SL1/SL2 types, consider that you may have mis-specified your SLs or your SL types may not be as well characterised as expected. 

#### Multiple SLs; all are specialised for resolving downstream operonic genes

Omit the `-S` option to automatically designate all SLs as "SL2-type":

    sloppr.sh -s sl_sequences.fasta

Since SL2:SL1 ratios are undefined in this case, every gene receiving any SL will be designated as downstream operonic. Each operon will also receive an upstream gene that is not SL *trans*-spliced. By definition, no monocistronic SL-receiving genes are designated.

Use the inferred SL clustering as a sanity check to verify absence of SL subfunctionalisation.

#### Multiple SLs; all are specialised for resolving upstream and downstream operonic genes

As above, but use the `-u` option to enforce SL2-type bias at upstream operonic genes:

    sloppr.sh -s sl_sequences.fasta -u

#### Multiple SLs; specialisation unknown 

Omit the `-S` option and carefully inspect the SL clustering results for plausible SL1/SL2-type subfunctionalisation:

    sloppr.sh -s sl_sequences.fasta

#### Multiple SLs; specialisation absent

In this scenario, both operonic and monocistronic genes receive the same SLs. Instead of using SL2:SL1 read ratios, we must tease out operonic genes via short intercistronic distances from monocistronic genes with large intergenic distances. 

First, omit the `-S` option and observe vastly inflated intercistronic distances (because every gene receiving any SL will be designated as downstream operonic):

    sloppr.sh -s sl_sequences.fasta

Then, run the automatic intercistronic distance filter, or specify a cutoff directly (e.g., 500 bp):

    sloppr.sh -s sl_sequences.fasta -i x
    sloppr.sh -s sl_sequences.fasta -i 500

These filters will cause SL-receiving genes with large intergenic distances to be designated as "monocistronic".

#### Single SL; specialised for resolving downstream operonic genes

Omit the `-S` option and ignore all SL clustering results:

    sloppr.sh -s sl_sequences.fasta

#### Single SL; specialised for resolving upstream and downstream operonic genes

As above, but use the `-u` option to enforce SLs at upstream operonic genes:

    sloppr.sh -s sl_sequences.fasta -u

#### Single SL; specialisation absent

As above, but use intercistronic distance filtering to designate monocistronic genes:

    sloppr.sh -s sl_sequences.fasta -i x
    sloppr.sh -s sl_sequences.fasta -i 500

#### I want to analyse hundreds of RNA-Seq libraries - can SLOPPR handle it?

Yes, but be aware of bottlenecks:
- Read alignment and SL quantification are the most time consuming steps. SLOPPR processes libraries sequentially, using all available threads at any one time. If you have an HPC cluster at hand, it would be much faster to set up a task array job for aligning reads manually (`hisat2 --no-softclip --no-discordant`) and then use the BAM files as input for SLOPPR.
- SL clustering in R may be time consuming for hundreds of samples

Future updates may support more efficient data structures and automatic HPC job control.

<a name="citation"></a>
# Citation

Please cite the bioRxiv preprint:

Marius A. Wenzel, Berndt Mueller, Jonathan Pettitt. SLIDR and SLOPPR: Flexible identification of spliced leader trans-splicing and prediction of eukaryotic operons from RNA-Seq data.
bioRxiv 2020.12.23.423594; [doi: https://doi.org/10.1101/2020.12.23.423594](https://doi.org/10.1101/2020.12.23.423594)


