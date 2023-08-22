#!/usr/bin/env bash
set -e

function timestamp { 
	date +"[%Y-%m-%d %H:%M:%S]"
}

function load_dependencies {
	# ensure that all dependencies are in PATH
	# MODIFY/DELETE THESE CALLS ACCORDING TO YOUR MACHINE
	echo -e "#\n$(timestamp) >>> Loading dependencies"
	export PATH="$PATH:~/sharedscratch/apps/ViennaRNA-2.4.14/bin"
	export PATH="$PATH:~/sharedscratch/apps"
	export PATH="$PATH:~/sharedscratch/apps/hisat2-2.2.1/"
	module load samtools
	module load blast-plus
	module load r
	module load cutadapt
	module load bowtie2
}

title="#\n# ***********************************************************\n# * SLIDR - Spliced Leader IDentification from RNA-Seq data *\n# * Version 1.2                                             *\n# ***********************************************************\n#"

function printhelp {
	echo -e "$title"
	echo -e "# Usage: slidr.sh [-o <>] [-p <>] [-c <>] -g <> | -t <> [-a <>]"
	echo -e "#   \t\t  {-1 <> [-2 <>] | -b <> | -m <>} [-r 0|1|2|x] [-q]"
	echo -e "#   \t\t  [-l <>] [-x <>] [-e <>] [-D <>] [-S <>] [-A <>] [-R <>] [-O <>] [-L <>]"
	echo -e "#\n# General options:"
	echo -e '#    -h            Print this help page'
	echo -e '#    -o <dir>      Path to output directory (default: "SLIDR_[date+time]")'
	echo -e '#    -p <chr>      Name prefix for predicted SL sequences (default: SL)'
	echo -e '#    -c <num>      CPU threads (default: 8)'
	echo -e '#    --tmp <dir>   TEMP directory (default: $TMPDIR)'
	echo -e '#    --hpc         Run on HPC with automatic job submission'
	echo -e '#    --hpcp <chr>  HPC submit command for parallel jobs (one job per library)'
	echo -e '#                  (default: "sbatch -c <threads> --mem <threads*4G> --parsable")'
	echo -e '#    --hpcd <chr>  HPC submit command for final dependent job'
	echo -e '#                  (default: "sbatch -c <threads> --mem <threads*8G> -d afterok:")'
	echo -e '#    --hpcs <chr>  Separator for job IDs in job dependency option (default: ":")'
	echo -e "#\n# Reference assembly:"
	echo -e '#    -g <file>     Path to genome assembly (FASTA[.gz])'
	echo -e '#    -a <file>     Path to genome annotations (GFF/GTF[.gz])'
	echo -e '#    -t <file>     Path to transcriptome assembly (FASTA[.gz])'
	echo -e "#\n# Single RNA-Seq library:"
	echo -e '#    -1 <file>     Path to R1 reads (FASTQ[.gz])'
	echo -e '#    -2 <file>     Path to R2 reads (FASTQ[.gz])'
	echo -e '#    -b <file>     Path to read alignments (BAM)'
	echo -e '#    -r <num>      Read strandedness (0=unstranded; 1=stranded; 2=reverse-stranded; x=infer; default: x)'
	echo -e "#    -q            Quality-trim 3' ends of reads and remove Illumina adapters"
	echo -e "#\n# Multiple RNA-Seq libraries:"
	echo -e '#    -m <file>     Path to tab-delimited file describing library configurations:'
	echo -e '#                  Column 1: Library name'
	echo -e '#                  Column 2: Read strandedness (= -r option)'
	echo -e '#                  Column 3: Path to R1 reads (= -1 option)'
	echo -e '#                  Column 4: Path to R2 reads (= -2 option)'
	echo -e '#                  Column 5: "trim" to quality-trim reads (= -q option)'
	echo -e '#                  Column 6: Path to read alignments (= -b option)'
	echo -e "#\n# Read tail clustering:"
	echo -e "#    -l <num>      Minimum tail length (default: 8)"
	echo -e "#    -x <num>      Scale factor for upper tail length limit (default: 1.0)"
	echo -e "#    --agc         Use abundance-based greedy clustering (default: use DGC)"
	echo -e "#\n# SL RNA filtering:"
	echo -e "#    -e <num>      BLASTN E-value (default: 1)"
	echo -e "#    -D <chr>      Splice donor site regex (default: 'GT')"
	echo -e "#    -A <chr>      Splice acceptor site regex (default: 'AG')"
	echo -e "#    -S <chr>      Sm binding site regex (default: '.{20,60}AT{4,6}G')"
	echo -e "#    -R <num>      Maximum SL RNA length excluding SL (default: 80)"
	echo -e "#    -O <num>      Maximum SL overlap with trans-splice acceptor site (default: 10)"
	echo -e "#    -L <num>      Maximum base-pair span within stem-loop (default: 35)"
	echo -e "#"
}

#
# PARSE COMMAND LINE OPTIONS AND SANITY CHECK INPUTS
#

# fundamental sanity check
if [[ $# -eq 0 ]] ; then
	printhelp
	exit 0
fi

# default parameters
outdir="SLIDR_$(date +"%Y-%m-%d_%H-%M-%S")"
threads=8 # $(nproc)
stranded="x"
clean="no"
tcull="no"
tlength=8
tscale=1.0
agc="DGC"
evalue=1
donor="GT"
sm=".{20,60}AT{4,6}G"
acceptor="AG"
rlength=80
overlap=10
sloop=35
slprefix="SL"
# HPC control
hpc="no"
plocal=""
hpc_submit=""
hpc_depend=""
hpc_sep=":"
jobdeps=""

# parse command-line options
cmdline="slidr.sh"
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
	
		-h)		printhelp
		exit 0;;
		-g)		genome=$2
				cmdline+=" -g $2"
		shift; shift;;
		-a)		ann=$2
				cmdline+=" -a $2"
		shift; shift;;
		-t)		transcriptome=$2
				cmdline+=" -t $2"
		shift; shift;;
		-T)		tcull="yes"
				cmdline+=" -T"
		shift;;
	    -1)		R1=$2
				cmdline+=" -1 $2"
		shift; shift;;
	    -2)		R2=$2
				cmdline+=" -2 $2"
		shift; shift;;
		-b)		bam=$2
				cmdline+=" -b $2"
		shift; shift;;
	    -r)		stranded=$2
				cmdline+=" -r $2"
		shift; shift;;
		-q)		clean="trim"
				cmdline+=" -q"
		shift;;
		-m)		design=$2
				cmdline+=" -m $2"
		shift; shift;;
		-l)		tlength=$2
				cmdline+=" -l $2"
		shift; shift;;
		-x)		tscale=$2
				cmdline+=" -x $2"
		shift; shift;;
		--agc)	agc="AGC"
				vsopt="--sizeorder"
				cmdline+=" --agc"
		shift;;
		-e)		evalue=$2
				cmdline+=" -e $2"
		shift; shift;;
		-D)		donor=$2
				cmdline+=" -D '$2'"
		shift; shift;;
		-S)		sm=$2
				cmdline+=" -S '$2'"
		shift; shift;;
		-A)		acceptor=$2
				cmdline+=" -A '$2'"
		shift; shift;;
		-R)		rlength=$2
				cmdline+=" -R $2"
		shift; shift;;
		-O)		overlap=$2
				cmdline+=" -O $2"
		shift; shift;;
		-L)		sloop=$2
				cmdline+=" -L $2"
		shift; shift;;
		-o)		outdir=$2
				cmdline+=" -o '$2'"
		shift; shift;;
		-p)		slprefix=$2
				cmdline+=" -p $2"
		shift; shift;;
	    -c)		threads=$2
				cmdline+=" -c $2"
		shift; shift;;
		--tmp)	export TMPDIR=$2
				cmdline+=" --tmp $2"
		shift; shift;;
		--hpc)	hpc="yes"
				cmdline+=" --hpc"
		shift;;
		--local)	plocal="yes"		# override hpc submission for final part of pipeline
		shift;;
		--hpcp)
				hpc="yes"
				hpc_submit=$2
				cmdline+=" --hpcp '$2'"
		shift; shift;;
		--hpcd)
				hpc="yes"
				hpc_depend=$2
				cmdline+=" --hpcd '$2'"
		shift; shift;;
		--hpcs)
				hpc="yes"
				hpc_sep=$2
				cmdline+=" --hpcs '$2'"
		shift; shift;;
	    -)	            # unknown option
		shift; shift;;
	esac
done

# HPC control

# ensure that automatic output directory is passed on
if [[ ! "$cmdline" =~ " -o " ]]; then 
		cmdline+=" -o '$outdir'"
fi

# if HPC is requested but submit command is omitted, use default SLURM command
if [ "$hpc" == "yes" ] && [ "$hpc_submit" == "" ]; then
	hpc_submit="sbatch -c $threads --mem $(($threads*4))G --parsable"
fi
# same for job dependency command
# if HPC is requested but submit command is omitted, use default SLURM command
if [ "$hpc" == "yes" ] && [ "$hpc_depend" == "" ]; then
	hpc_depend="sbatch -c $threads --mem $(($threads*8))G -d afterok:"
fi

# print input summary on screen and log file
function print_summary {
	echo -e "$title"
	echo -e "# General:"
	echo -e "#   > Output directory\t $outdir"
	echo -e "#   > SL name prefix\t $slprefix"
	echo -e "#   > Threads\t\t $threads"
	echo -e "#   > TEMP directory\t $TMPDIR"
	if [ "$hpc" == "yes" ]; then
		echo -e "#   > HPC parallel\t $hpc_submit"
		echo -e "#   > HPC sequential\t $hpc_depend"
	fi
	echo -e "#\n# Reference:"
	echo -e "#   > Genome\t\t $genome"
	echo -e "#   > Annotations\t $ann"
	echo -e "#   > Transcriptome\t $transcriptome"
	echo -e "#\n# Read tail clustering:"
	echo -e "#   > Minimum length\t $tlength"
	echo -e "#   > Maximum length\t ${tscale}X"
	echo -e "#   > Clustering mode\t $agc"
	echo -e "#\n# SL RNA filtering:"
	echo -e "#   > BLASTN E-value\t $evalue"
	echo -e "#   > Splice donor\t $donor"
	echo -e "#   > Splice acceptor\t $acceptor"
	echo -e "#   > Sm binding site\t $sm"
	echo -e "#   > SL RNA length\t $rlength"
	echo -e "#   > Acceptor overlap\t $overlap"
	echo -e "#   > Stem-loop span\t $sloop"
	echo -e "#\n# RNA-Seq data:"
	if [ ! "$design" == "" ] ; then
		echo -e  "#   > Libraries\t\t $design"
	else
		echo -e  "#   > R1 reads\t\t $R1"
		echo -e  "#   > R2 reads\t\t $R2"
		echo -e  "#   > Quality trimming\t $clean"
		echo -e  "#   > BAM alignments\t $bam"
		echo -e  "#   > Strandedness\t $stranded"
	fi
}
print_summary

# sanity check input
function sanity_check {
	# we need to resolve "~" in file paths supplied via text file
	bam="${bam/\~/$HOME}"
	R1="${R1/\~/$HOME}"
	R2="${R2/\~/$HOME}"
	
	# if transcriptome is supplied, splice acceptor site must be empty
	if [ ! "$transcriptome" == "" ]; then
		if [ ! "$acceptor" == "" ]; then
			echo "WARNING: Splice acceptor sites cannot be identified in a transcriptome reference. Removing splice acceptor motif."
		fi
		if [ "$overlap" -gt 0 ]; then
			echo "WARNING: It is recommended to set trans-splice acceptor overlap to 0 when using a transcriptome reference."
		fi
		acceptor=""
	fi
	
	if ([ "$genome" == "" ] && [ "$transcriptome" == "" ]) || \
		([ ! "$genome" == "" ] && [ ! "$transcriptome" == "" ]) ; then
		echo "ERROR: You must specify either a genome or a transcriptome! Exiting."
		exit 1
	elif [ ! -f $genome ]; then
		echo "ERROR: The genome $genome does not exist! Exiting."
		exit 1
	elif [ ! -f $transcriptome ]; then
		echo "ERROR: The transcriptome $transcriptome does not exist! Exiting."
		exit 1
	fi
	if [ ! "$ann" == "" ] && [ ! "$transcriptome" == "" ]; then
		echo "WARNING: The annotation file will be ignored when a transcriptome is supplied."
	fi
	if [ ! -f $ann ]; then
		echo "ERROR: The annotation file $ann does not exist! Exiting."
		exit 1
	fi
	if [ "$bam" == "" ] && [ "$R1" == "" ] && [ "$design" == "" ]; then
		echo "ERROR: You must specify either a BAM, R1 reads or library configuration file! Exiting."
		exit 1
	fi
	if [ ! "$bam" == "" ] && [ ! -f $bam ]; then
		echo "ERROR: The BAM file $bam does not exist! Exiting."
		exit 1
	fi
	if [ ! "$R1" == "" ] && [ ! -f $R1 ]; then
		echo "ERROR: R1 FASTQ file $R1 does not exist! Exiting."
		exit 1
	fi
	if [ ! "$R2" == "" ] && [ ! -f $R2 ]; then
		echo "ERROR: R2 FASTQ file $R2 does not exist! Exiting."
		exit 1
	fi
	if [ ! "$design" == "" ] && [ ! -f $design ]; then
		echo "ERROR: The library configuration file $design does not exist! Exiting."
		exit 1
	fi
	if [ ! "$bam" == "" ] && ([ ! "$R1" == "" ] || [ ! "$R2" == "" ]); then
		echo "WARNING: You have specified both a BAM file and reads. The reads will be ignored."
	fi
}
# first sweep of checks
sanity_check

# if library design was supplied, sanity check that, too:
if [ ! "$design" == "" ] && ([ ! "$bam" == "" ] || [ ! "$R1" == "" ] || [ ! "$R2" == "" ]); then
	echo "You have specified a library configuration file. Specified BAM or reads files will be ignored."
fi
if [ ! "$design" == "" ]; then
	while IFS='#' read lib stranded R1 R2 clean bam
	do
		echo "#   * $lib (-r $stranded) -1 $R1 -2 $R2 $clean $bam"
		sanity_check
	done < <(tr -d '#' < $design | tr '\t' '#' )
fi

#
# END OF SANITY CHECKS
#
load_dependencies
mkdir -p "$outdir"
echo "$cmdline" > "$outdir/0-command_summary.txt"
print_summary >> "$outdir/0-command_summary.txt"

#
# 1) PREPARE GENOME/TRANSCRIPTOME INDEX
#
echo "$(timestamp) >>> STAGE 1: SL read identification"

if [ ! "$genome" == "" ]; then
	ref="$genome"
elif [ ! "$transcriptome" == "" ]; then
	ref="$transcriptome"
fi

# copy genome/transcriptome
# unzip if required
if [ ! -f "$outdir/reference.fa" ]; then
	if [ "$(file $ref | grep -i gzip)" ]; then
		gunzip -c "$ref" > "$outdir/reference.fa"
	else
		cp "$ref" "$outdir/reference.fa"
	fi
fi
ref="$outdir/reference.fa"

# copy annotations
# unzip if required
# convert GFF to GTF if required
if [ ! "$ann" == "" ]; then
	gtf="$outdir/annotations.gtf"
	if [ ! -f "$gtf" ]; then
		echo "$(timestamp) Processing genome annotations ..."
		gropts="--force-exons --gene2exon --keep-genes -T"
		fform="$(file $ann)"
		if [[ "$fform" =~ "gzip" ]] && [[ "$fform" =~ ".gtf" ]]; then
			# gzipped GTF
			gunzip -c "$ann" > "$gtf"
		elif [[ "$fform" =~ "gzip" ]] && [[ ! "$fform" =~ ".gtf" ]]; then
			# gzipped GFF
			gunzip -c "$ann" | gffread $gropts -o "$gtf"
		elif [[ ! "$fform" =~ "gzip" ]] && [[ "$fform" =~ ".gtf" ]]; then
			# unzipped GTF
			cp "$ann" "$gtf"
		else
			# unzipped GFF
			gffread "$ann" $gropts -o "$gtf"
		fi
	fi
	ann="$gtf"
	bedtools sort -i "$ann" | awk '$3~"transcript|gene|mRNA"' > "$ann.genes"
fi

# make hisat2 index for genome
if [ ! "$genome" == "" ] && [ ! -d "$outdir/hisat2_index" ]; then
	echo "$(timestamp) Generating HISAT2 index ..."
	mkdir -p "$outdir/hisat2_index"
	# if we have genome annotations, use them
	if [ ! "$ann" == "" ]; then
		# extract splice sites for hisat2
		hisat2_extract_splice_sites.py "$ann" > "$outdir/hisat2_index/hisat2_splicesites.txt"
		hisat2_extract_exons.py "$ann" > "$outdir/hisat2_index/hisat2_exons.txt"
		if [ -s "$outdir/hisat2_index/hisat2_splicesites.txt" ] && [ -s "$outdir/hisat2_index/hisat2_exons.txt" ]; then
			hisatgffopts+="--ss $outdir/hisat2_index/hisat2_splicesites.txt --exon $outdir/hisat2_index/hisat2_exons.txt"
		fi
	fi
	hisat2-build -p $threads $hisatgffopts "$ref" "$outdir/hisat2_index/genome" \
		> "$outdir/hisat2_index/log_hisat2-build.txt" 2>&1
fi


#
# IMPLEMENT TRANSCRIPTOME CULLING HERE
#

# make bowtie2 index for transcriptome
if [ ! "$transcriptome" == "" ] && [ ! -d "$outdir/bowtie2_index" ]; then
	echo "$(timestamp) Generating BOWTIE2 index ..."
	mkdir -p "$outdir/bowtie2_index"
	bowtie2-build --threads $threads \
		"$ref" "$outdir/bowtie2_index/transcriptome" \
		> "$outdir/bowtie2_index/log_bowtie2-build.txt" 2>&1
fi


# we need to extract tails for each library
# if a library-design file is specified, loop through libraries
# for each library, check what needs done
# new: check whether library.OK file exists. If not, submit job to check what needs done

function library_pipeline {
	echo "$(timestamp) >>> Processing library $lib"
	lib="1-library_$lib"
	# does library.OK file exist? If yes, library has already been processed
	if [ -f "$outdir/$lib/library-x$tscale.OK" ]; then 
		echo "$(timestamp) >>> Nothing to do"
	else
		# process via HPC or local
		if [ "$hpc" == "yes" ]; then
			# HPC
			jobID=$(eval "$hpc_submit slidr_library.sh -var=\"$threads\" -var=\"$outdir\" -var=\"$lib\" -var=\"$stranded\" -var=\"$tscale\" -var=\"$R1\" -var=\"$R2\" -var=\"$clean\" -var=\"$bam\"")
			echo "$(timestamp) >>> Submitted HPC job $jobID"
			jobdeps+="$jobID$hpc_sep"
		else
			# local
			slidr_library.sh -var="$threads" -var="$outdir" -var="$lib" -var="$stranded" -var="$tscale" -var="$R1" -var="$R2" -var="$clean" -var="$bam"
		fi
	fi
}

if [ ! "$design" == "" ]; then
	while IFS='#' read lib stranded R1 R2 clean bam
	do
		library_pipeline
	done < <(tr -d '#' < "$design" | tr '\t' '#')
	else
		# just do it once
		lib="data"
		library_pipeline
fi

### from here on, everything will be done on all libraries together
# if HPC was not used, all libraries will be done and we can continue straight away
# if HPC was used, it's complicated
# if at least one library was submitted, we'll have to wait until that job is finished
# if no libraries were submitted, we are ready 

if [ "$hpc" == "yes" ] && [ ! "$plocal" == "yes" ]; then
	# submit self and exit
	# check if dependencies are empty; in that case we don't need them!
	# the dependency command will then just get a dummy job ID (1 for SLURM)
	jobdeps=${jobdeps%$hpc_sep}
	if [ "$jobdeps" == "" ]; then
		jobdeps="1"
	else
		echo "$(timestamp) >>> Pipeline will continue once HPC jobs $jobdeps have finished"
	fi
	eval "${hpc_depend}${jobdeps} $cmdline --local"
	echo "$(timestamp) >>> Submitted current run to HPC. Exiting."
	exit 0
fi

# continue rest of pipeline as normal

echo "$(timestamp) >>> STAGE 2: SL tail clustering"

#
# 4) TAIL CLUSTERING
#

mkdir -p $outdir/1-tails

# collect all tails
if [ ! -f $outdir/1-tails/tails-x$tscale.fa.gz ]; then
	echo "$(timestamp) Collecting tails ..."
	cat $outdir/1-library_*/alignments-x$tscale.bam.tails*p.fa.gz \
		> $outdir/1-tails/tails-x$tscale.fa.gz
fi

# dereplicate 
# remove sequences with Ns
tailprefix=$outdir/1-tails/tails-x$tscale-l$tlength
if [ ! -f $tailprefix.derep.fa.gz ]; then
	echo "$(timestamp) Dereplicating tails ..."
	gunzip -c $outdir/1-tails/tails-x$tscale.fa.gz \
	| awk '!/^>/ && $0~"N"{$0=""}{print}' \
	| vsearch --derep_fulllength - \
		--output $tailprefix.derep.fa \
		--uc $tailprefix.derep.out.txt \
		--sizeout \
		--minseqlength $tlength \
		--threads $threads \
		> $tailprefix.derep.log-vsearch.txt 2>&1

	# clusterinfo; two columns: representative, read
	awk '$1=="H"{print $10, $9}$1=="S"{print $9, $9}' OFS="\t" $tailprefix.derep.out.txt \
		| sort -S90% --batch-size=100 -t "$(printf '\t')" -k 1,1 \
		| gzip > $tailprefix.derep.clusterinfo.txt.gz
	gzip $tailprefix.derep.fa
	gzip $tailprefix.derep.out.txt
	
	ns=$(grep -o '[0-9]* seqs' $tailprefix.derep.log-vsearch.txt | grep -o '[0-9]*')
	nt=$(grep -o '[0-9]* unique' $tailprefix.derep.log-vsearch.txt | grep -o '[0-9]*')
	echo "$(timestamp) Dereplicated $ns tails (>= $tlength bp) to $nt unique tails"	
fi

# cluster
clusterprefix=$tailprefix.cluster-$agc 
if [ ! -f $clusterprefix.centroids.fa.gz ]; then
	echo "$(timestamp) Clustering tails ($agc) ..."
	# --wordlength 8 --minwordmatches 1000 ensures that all words must match
	# --sizeorder --maxaccepts 0 --maxrejects 0 ensures that AGC finds the most abundant sequence at 100% identity
	gunzip -c $tailprefix.derep.fa.gz \
	| vsearch --cluster_fast - --qmask none \
		--uc $clusterprefix.out --centroids $clusterprefix.centroids.fa \
		--consout $clusterprefix.consensus.fa --msaout $clusterprefix.msa.fa \
		--clusterout_sort --sizein --sizeout \
		--minseqlength $tlength --threads $threads \
		--id 1.0 --iddef 0 --rightjust \
		--wordlength 8 --minwordmatches 1000 \
		$vsopt --maxaccepts 0 --maxrejects 0 \
		> $clusterprefix.log-vsearch.txt 2>&1
	
	# clusterinfo; two columns: centroid, representative
	# sort by representative to aid joining
	awk '$1=="H"{print $10, $9}$1=="S"{print $9, $9}' OFS="\t" $clusterprefix.out \
		| sed 's/;size=[0-9]\+//g' \
		| sort -S90% --batch-size=100 -t "$(printf '\t')" -k 2,2 \
		| gzip > $clusterprefix.clusterinfo.txt.gz
	gzip $clusterprefix.out
	gzip $clusterprefix.centroids.fa
	gzip $clusterprefix.msa.fa
	
	echo "$(timestamp) Generated $(zgrep -c '^>' $clusterprefix.centroids.fa.gz) tail clusters"
fi

# summarise read information from derep and cluster
# two columns: read, centroid
if [ ! -f $clusterprefix.readinfo.txt.gz ]; then
	echo "$(timestamp) Extracting cluster information ..."
	join -e '-' -t "$(printf '\t')" -1 1 -2 2 \
		<(gunzip -c $tailprefix.derep.clusterinfo.txt.gz) \
		<(gunzip -c $clusterprefix.clusterinfo.txt.gz) \
		| cut -f 2,3 | gzip > $clusterprefix.readinfo.txt.gz
fi

#
# 5) TAIL FILTERING
#
echo "$(timestamp) >>> STAGE 3: SL RNA identification"
mkdir -p $outdir/2-RNA_filters
# map to the reference
# ensure that the 3' end of the tail is aligned (5' end can contain noise, so not important) $8==$14
# have considered keeping only the longest match(es) for each tail (not done currently): 
# awk '$8==$14{ if($1==q){if($4==l)print} else {q=$1;l=$4; print} }'

echo "$(timestamp) Aligning tail cluster centroids to reference ..."
if [ ! -d $outdir/blast_db ]; then
	mkdir -p $outdir/blast_db
	makeblastdb -dbtype nucl -in $ref -out $outdir/blast_db/reference > $outdir/blast_db/log_makeblastdb.txt 2>&1
fi

blastprefix=$outdir/2-RNA_filters/centroids-x$tscale-l$tlength-$agc.blastn-e$evalue
if [ ! -f $blastprefix.out.gz ]; then
	gunzip -c $clusterprefix.centroids.fa.gz \
		| dustmasker -in - -level 10 -window 10 -outfmt fasta \
		| blastn -db $outdir/blast_db/reference -query - \
			-word_size 8 -num_threads $threads -perc_identity 100 -lcase_masking -evalue $evalue \
			-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq qlen' \
			2> $blastprefix.log_blastn.txt \
			| awk '$8==$14' > $blastprefix.out
	gzip $blastprefix.out

	nbl=$(gunzip -c $blastprefix.out.gz | wc -l)
	echo "$(timestamp) Generated $nbl alignment records"
fi

#
# 6) SPLICE DONOR AND SM BINDING SITES
#
# find splice donor and Sm binding site
# if we only have a transcriptome, this could still work if we have full SL RNAs
# the code is robust to empty donor or Sm sites; it will simply return empty overlaps and the full $rlength bp.
echo "$(timestamp) Identifying splice donor and Sm binding sites (${donor}${sm}) ..."

smout="$blastprefix.slRNA-R$rlength-D$donor-S$sm.txt.gz"
if [ ! -f "$smout" ]; then
	samtools faidx $outdir/reference.fa
	gunzip -c $blastprefix.out.gz \
	 | awk '$9<$10{print $2, $10, $10+1, $13"::"$1, ".", "+"} 
			$9>$10{print $2, ($10-2>0?$10-2:0), ($10-1>1?$10-1:1), $13"::"$1, ".", "-"}' OFS="\t" \
	 | bedtools slop -g $outdir/reference.fa.fai -i - -s -l 0 -r $(($overlap+$rlength)) \
	 | bedtools getfasta -name+ -tab -s -fi $outdir/reference.fa -bed - \
	 | awk --re-interval -F'\t|::' -v p="${donor}${sm}" -v l="$rlength" -v x="$overlap" \
		'{  s=match(toupper($4), toupper(p)); o=toupper(substr($4, 1, s-1)); if(o==""){o=""}; gsub(":", "\t", $3); 
			if(s>0 && s<=(x+1)){print $2, $3, $1, o, $1""o""substr(toupper($4),RSTART,RLENGTH+l)}}' OFS="\t" \
		| sed 's/;size=[0-9]\+/\t&/g' | gzip > "$smout"
		# have the ;size= information in a separate column to facilitate read ID merging later
	
	nsm=$(gunzip -c "$smout" | wc -l)
	nbl=$(gunzip -c $blastprefix.out.gz | wc -l)
	echo "$(timestamp) Identified $nsm donor splice sites and Sm binding sites ($((100*$nsm/$nbl))% of alignments)"
fi

#
# 7) SPLICE ACCEPTOR SITES
#
# we can find acceptor sites when we have a genome.
# this won't work with a transcriptome, but we still need to collect read alignment locations.
# filter by acceptor site motif.
tailbams="$outdir/2-RNA_filters/alignments-x$tscale"
spliceout="$tailbams.spliceacceptor-$acceptor-O$overlap.txt.gz"
if [ ! -f "$spliceout" ]; then
	echo "$(timestamp) Extracting splice acceptor sites ($acceptor) ..."
	
	# collect read alignments
	{ if [ ! -f "$tailbams.tails+5p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails+5p.bam $outdir/1-library_*/*x$tscale.bam.tails+5p.sam.gz
	fi; } &
	{ if [ ! -f "$tailbams.tails-3p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails-3p.bam $outdir/1-library_*/*x$tscale.bam.tails-3p.sam.gz
	fi; } &
	{ if [ ! -f "$tailbams.tails-5p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails-5p.bam $outdir/1-library_*/*x$tscale.bam.tails-5p.sam.gz
	fi; } &
	{ if [ ! -f "$tailbams.tails+3p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails+3p.bam $outdir/1-library_*/*x$tscale.bam.tails+3p.sam.gz
	fi; } & wait
		
	# length of acceptor site (AG: 2bp)
	# if regex is used, try and work out how many characters using bracket pairs
	# e.g., G[TA][AC] would be three characters
	acclen=$(($(echo "$acceptor" | sed 's/\[[^][]*\]/./g' | wc -m)-1))
		
	# collect acceptor regions and filter by acceptor site
	# only use reads that were clustered (shorter reads than $tlength were dropped)
	#readids=$(zless $outdir/1-tails/tails.derep.clusterinfo.txt.gz | cut -f 2)
	#clinfo="$tailprefix.derep.clusterinfo.txt.gz"
	#gunzip -c $clinfo | cut -f 2 > $tailbams.reads.tmp
	
	{ bedtools bamtobed -i $tailbams.tails+5p.bam \
		| awk '{$3=$2+1; print}' OFS="\t" > $tailbams.tails+5p.bed; } &
	{ bedtools bamtobed -i $tailbams.tails-3p.bam \
		| awk '{$3=$2+1; $6="+"; print}' OFS="\t" > $tailbams.tails-3p.bed; } &
	{ bedtools bamtobed -i $tailbams.tails-5p.bam \
		| awk '{$2=$3-1; print}' OFS="\t" > $tailbams.tails-5p.bed; } &
	{ bedtools bamtobed -i $tailbams.tails+3p.bam \
		| awk '{$2=$3-1; $6="-"; print}' OFS="\t" > $tailbams.tails+3p.bed; } & wait
	
	cat $tailbams.tails*.bed \
		| bedtools slop -g $outdir/reference.fa.fai -i - -s -l $acclen -r $overlap \
		| bedtools getfasta -fi $outdir/reference.fa -bed - -s -tab -name+ \
		| awk -F'\t|::' -v a="$acceptor" 'toupper($3)~toupper(a){
				gsub(":", "\t", $2); 
				gsub("/[12]$", "", $1); print $1, $2, toupper($3)
			}' OFS="\t" \
		| gzip > "$spliceout"
	
	rm $tailbams.tails*.bed
	#rm $tailbams.reads.tmp
	
	nas=$(gunzip -c $spliceout | wc -l)
	echo "$(timestamp) Extracted $nas splice acceptor sites"
fi

#
# 8) RNAfold on unique donor sequences
#
rnafout="$smout.RNAfold-L$sloop.txt.gz"
if [ ! -f "$rnafout" ]; then
	echo "$(timestamp) Predicting RNA secondary structures ..."
	# we want the Sm site to be unpaired, so need to add a constraint string (...xxxx) below the sequence
	donlen=$(($(echo "$donor" | sed 's/\[[^][]*\]/./g' | wc -m)-1))
	smlen=$(($(echo "$sm" | sed -e 's/\[[^][]*\]/./g' -e 's/{[^}{]*}/./g' | wc -m)-2))
	gunzip -c "$smout" | cut -f 7 | sort | uniq \
		| awk -v d="$donlen" -v s="$smlen" -v p="${donor}${sm}" '{match(toupper($0), p); 
			print ">"$0"\n"$0; 
			c=1; if(s<0){s=0};
			while(c<(RSTART+RLENGTH-s)) { printf "."; c++ }
			while(c<(RSTART+RLENGTH)) { printf "x"; c++ }
			while(c<=length($0)) { printf "."; c++ }
			print "\n"}' \
			| RNAfold -j$threads -C --noPS --noDP --noconv -p --MEA -d0 --maxBPspan=$sloop --noLP \
			| awk '{if($0~/^>/){gsub(">", "", $0); d=$0}
					else if($0~/mfe/){gsub(";", "", $7); m=$7; e=$10; print d,l,m,e }
					else if($0~/)$/){l=$1}}' OFS="\t" \
			| gzip > "$rnafout"
fi
ud=$(gunzip -c "$rnafout" | wc -l)
echo "$(timestamp) Processed $ud unique RNAs"

#
# 9) reconstruct 3' end of SL tails and cluster
#
echo "$(timestamp) >>> STAGE 4: Consensus SL construction"

# new: cluster directly from slRNA output. No need to go through all the joining trouble.
if [ ! -f "$smout.centroids.fa.gz" ]; then
	echo "$(timestamp) Re-clustering SL candidates ..."
	gunzip -c "$smout" | awk -F'\t' '{print ">"$1"::"$5"::"$6$2"::nl::"$5$6}'  \
		| sort -S90% -t "$(printf '\t')" | uniq | sed 's/::nl::/\n/g' | gzip > "$smout.SLtails.fa.gz"
	gunzip -c "$smout.SLtails.fa.gz" \
	| vsearch --cluster_fast - --qmask none \
		--uc "$smout.clusters.out" --centroids "$smout.centroids.fa" \
		--msaout "$smout.msa.fa" \
		--id 1.0 --iddef 0 --rightjust \
		--clusterout_sort --sizein --sizeout \
		--minseqlength $tlength --threads $threads \
		--wordlength 8 --minwordmatches 1000 \
		--sizeorder --maxaccepts 0 --maxrejects 0 \
		> "$smout.log_vsearch-cluster.txt" 2>&1
		
	awk -F'::|\t' '$1~"S|H"{print $2, $9, $10, $11}' OFS="\t" "$smout.clusters.out" \
		| sed 's/;size=[0-9]\+//g' | gzip > "$smout.clusterinfo.txt.gz"
	gzip "$smout.clusters.out"
	gzip "$smout.centroids.fa"
	gzip "$smout.msa.fa"
fi
nc=$(grep -o "Clusters: [0-9]*" $smout.log_vsearch-cluster.txt | grep -o "[0-9]*")
echo "$(timestamp) Generated $nc clusters"

#
# 10) final filtering and consensus construction in R
#
echo "$(timestamp) Resolving consensus SLs ..."
resultsdir="$outdir/3-results-x$tscale-l$tlength-$agc-e$evalue-R$rlength-D$donor-S$sm-L$sloop-A$acceptor-O$overlap"
mkdir -p "$resultsdir"
# deposit file of input file names
# read info
echo "$clusterprefix.readinfo.txt.gz" > "$resultsdir/filters.fofn"
# acceptor info
echo "$spliceout" >> "$resultsdir/filters.fofn"
# donor info
echo "$smout" >> "$resultsdir/filters.fofn"
# RNAfold info
echo "$rnafout" >> "$resultsdir/filters.fofn"
# SL tail clusterinfo
echo "$smout.clusterinfo.txt.gz" >> "$resultsdir/filters.fofn"

slidr_consensus.R "$resultsdir" " $acceptor" $slprefix $threads

#
# 11) Overlap predicted SL RNA genes and SLTS acceptor sites with reference annotations (if available)
#
if [ -f "$ann" ]; then
	echo "$(timestamp) Finding SL RNA genes and trans-spliced genes in reference ..."

	for sl in $(find "$resultsdir" -name "*.RNA_genes.gff3" | sort)
	do
		acc=${sl/RNA_genes.gff3/acceptor_sites.gff3}
		sln=$(basename ${sl/.RNA_genes.gff3/})
		
		bedtools intersect -s -wb \
			-a <(bedtools sort -i "$sl") -b "$ann.genes" \
			> "$sl.reference.txt"
	
		bedtools intersect -s -wb \
			-a <(bedtools sort -i "$acc" | bedtools slop -r 100 -l 0 -s -i - -g "$ref.fai") -b "$ann.genes" \
			> "$acc.reference.txt"
	
		nslrnap=$(grep -c -v '^#' "$sl")
		nsltsap=$(grep -c -v '^#' "$acc")
		nslrna=$(wc -l < "$sl.reference.txt")
		nsltsa=$(wc -l < "$acc.reference.txt")
		
		echo "$(timestamp) Found $nslrna out of $nslrnap ($((100*$nslrna/$nslrnap))%) predicted $sln RNA genes ..."	
		echo "$(timestamp) Found $nsltsa $sln trans-spliced genes ..."	
		#echo "$(timestamp) Found $nsltsa out of $nsltsap ($((100*$nsltsa/$nsltsap))%) predicted $sln trans-spliced genes/transcripts ..."	
	done
fi

echo "$(timestamp) Finished!"

exit 0