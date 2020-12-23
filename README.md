# SLIDR and SLOPPR: Flexible identification of spliced leader *trans*-splicing and prediction of eukaryotic operons from RNA-Seq data

## General options for both pipelines

Both SLIDR and SLOPPR share a number of general options for data input and output, which are described in this first section. Pipeline-specific options are detailed in the following sections.

    bash slidr.sh -h
    bash sloppr.sh -h

`-o <dir>`
Path to output directory. If unspecified, the output directory is "./SLIDR_[date+time]" or "./SLOPPR_[date+time]"

`-c <num>`
Threads (default: all available cores obtained by `nproc`)

### RNA-Seq data input
SLIDR and SLOPPR accept single-end or paired-end RNA-Seq reads in FASTQ(.gz) or read alignments in BAM format. The following options are available to specify a single library:

`-1 <file>`
Path to R1 reads in FASTQ(.gz) or FASTA(.gz) format.

`-2 <file>`
Path to R2 reads in FASTQ(.gz) or FASTA(.gz) format.

`-q`
Basic quality-trimming of 3' ends of reads and removal of Illumina adapters using `cutadapt -a AGATCGGAAGAGC -q 20 -m 20`. It is recommended that raw data are inspected and trimmed thoroughly prior to running SLIDR or SLOPPR.

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

Leave columns empty if they are not required, for example:

    Library1 <tab> x <tab> R1.trimmed.fq.gz <tab> R2.trimmed.fq.gz
    Library2 <tab> 0 <tab> reads.fq.gz <tab> <tab> trim
    Library3 <tab> 2 <tab> <tab> <tab> <tab> alignments.bam

This configuration file allows great flexibility in mixing single-end and paired-end libraries with different strandedness and existing read alignments if desired.
Note: Avoid `#` characters in the file contents; they will be stripped! 

## SLIDR: Spliced leader identification from RNA-Seq

### Reference assembly

SLIDR requires either a genomic or transcriptomic reference assembly. A genomic reference is recommended and will produce vastly superior results. Genome annotations are recommended because they improve read mapping in HISAT2 and allow SLIDR to infer strandedness of the RNA-Seq data if unknown.

`-g <file>`
Path to genome assembly in FASTA format.

`-a <file>`
Path to genome annotations in GFF or GTF format.

`-t <file>`
Path to transcriptome assembly in FASTA format.

### Read tail filtering

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

### SL RNA filtering

`-D <chr>`
Splice donor site pattern (default: GT). Alternative nucleotides can be coded with character classes regex patterns, for example, `-D G[TC]` matches GT or GC, and `-D A[AG][TC]` matches AAT, AAC, AGT or AGC. To switch off, specify empty character string (`-D ''`).

`-S <chr>`
Sm binding site pattern in regex notation. This allows for searching additional motifs downstream of the splice donor site. The default (.{40,55}AT{4,6}G) matches the Sm binding sites ATTTTG, ATTTTTG or ATTTTTTG between 40 and 55 bases downstream of the splice donor. Any regex patterns can be searched, for example, `-S .{40,55}AT{4,6}G.{20,30}T{3,5}` to add a T-rich region 20-30 bp downstream of the Sm binding site. To switch off, specify empty character string (`-S ''`).

`-A <chr>`
Splice acceptor site pattern (default: AG). Alternative nucleotides can be coded with character classes regex patterns, for example, `-A A[CG]` matches GT or GC. To switch off, specify empty character string (`-A ''`)

`-R <chr>`
Maximum SL RNA length excluding SL (default: 80). This length is measured starting from the splice donor site. The default of 80 bp is appropriate for nematode SL RNAs that are c. 100 bp long including a c. 22 bp SL.

`-O <chr>`
Maximum SL overlap with trans-splice acceptor site (default: 10). This parameter controls how much overlap is allowed between the 3' end of the SL and the trans-splice acceptor site. If a transcriptome reference is used, it might be prudent to decrease this value or even set it to 0 if the transcripts do not contain extra non-coding 5' bases.

`-L <chr>`
Maximum base-pair span within stem-loop (default: 35). This parameter controls stem-loop prediction in RNAFold. The default requires the first and last base of a stem loop to be no more then 35 bp apart.

`-p <chr>`
Prefix for predicted SL sequences (default: SL). Output files will use this prefix for naming the SLs. 

## SLOPPR: Spliced leader-informed operon prediction from RNA-Seq

### Reference assembly
SLOPPR requires a genomic reference assembly and genome annotations in GFF or GTF format.

`-g <file>`
Path to genome assembly in FASTA format.

`-a <file>`
Path to genome annotations in GFF or GTF format.

### SL screening

`-s <file>`
Path to SL sequences in FASTA format. Ensure that each SL has a unique name (FASTA header).

`-n <num>`
Minimum bp from 3' end of SL required to detect SL tail in read (default: 8). Increase this value if the SLs cannot be distinguished reliably at their 3'-most 8 bases.

`-e <num>`
Maximum error rate for SL tail matching (default: 0.09). The default error rate requires no mismatches for tails up to 10 bp total length and allows for 1 mismatch for each 10 bp additional length (0-10 bp: 0; 11-21 bp: 1; 22-32 bp: 2; 33-40 bp: 3).

`-f <chr>`
Meta-feature ID field in GFF/GTF annotations (default: gene_id). Reads are quantified against gene annotations by default, but some genome annotations may require to set a different meta-feature ID, for example transcript_id. 

### Operon prediction

`-z [geo|sum|median]`
Method for aggregating SL counts across libraries (default: geometric mean `-z geo`). The geometric mean is an appropriate statistic to summarise counts, but we also implement sum and median, which may be more appropriate if zero counts are important (see below).

`-0` 
Keep libraries with zero counts when aggregating SL counts (default: remove zeros). Since the geometric mean is zero when at least one library has zero counts, the default settings disregard zero counts. Using the geometric mean with zero counts retained will only retain genes where all replicate libraries have non-zero counts. This is a very stringent filter that may be useful for removing lowly and inconsistently expressed genes.

`-S <file>`
Path to list of SL names that resolve polycistrons default: all SLs supplied with `-s`). This optional file must contain the names or FASTA headers of SLs (from those supplied by `-s`) that are specialised for resolving polycistrons (SL2-type SLs). Omitting this option will designate all SLs as SL2-type, which is useful for organisms where no specialisation of SLs exists (for example, tunicates). This option has no effect on the SL clustering algorithm, which infers SL2-type SLs independently. SLOPPR always carries out three prediction runs, using either SL2-type SLs, SLs in cluster1 or SLs in cluster2 as polycistron resolvers. 

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
