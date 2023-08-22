#!/usr/bin/env bash
set -e

function timestamp { 
	date +"[%Y-%m-%d %H:%M:%S]"
}

function load_dependencies {
	# ensure that all dependencies are in PATH
	# MODIFY/DELETE THESE CALLS ACCORDING TO YOUR MACHINE
	echo -e "#\n$(timestamp) >>> Loading dependencies"
	export PATH="$PATH:~/sharedscratch/apps"
	module load hisat2
	module load samtools
	module load bedtools2
	module load seqtk
	module load r
	module load cutadapt
	module load subread
	module load bowtie2
}

title="#\n# ************************************************************************\n# * SLOPPR - Spliced-Leader-informed Operon Prediction from RNA-Seq data *\n# * Version 1.2                                                          *\n# ************************************************************************\n#"

function printhelp {
	echo -e "$title"
	echo -e "# Usage: sloppr.sh [-o <>] [-c <>] -g <> -a <>"
	echo -e "# \t\t   {-1 <> [-2 <>] | -b <> | -m <>} [-r 0|1|2|x] [-q]"
	echo -e "# \t\t   -s <> [-n <>] [-e <>] [f <>]"
	echo -e "# \t\t   [-z geo|sum|median] [-0] [-S <>] [-d <>] [-u] [-i <>] [-p <>] [-x <>]"
	echo -e "#\n# General options:"
	echo -e '#    -h            Print this help page'
	echo -e '#    -o <dir>      Path to output directory (default: "SLOPPR_[date+time]")'
	echo -e "#    -p <chr>      Name prefix for predicted operons (default: OP)"
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
	#echo -e '#    -t <file>     Path to transcriptome assembly (FASTA[.gz])'
	echo -e "#\n# Single RNA-Seq library:"
	echo -e '#    -1 <file>     Path to R1 reads (FASTQ[.gz])'
	echo -e '#    -2 <file>     Path to R2 reads (FASTQ[.gz])'
	echo -e '#    -b <file>     Path to read alignments (BAM)'
	echo -e '#    -r <num>      Read strandedness (0=unstranded; 1=stranded; 2=reverse-stranded; x=infer; default: x)'
	echo -e "#    -q            Quality-trim 3' ends of reads and remove adapters"
	echo -e "#\n# Multiple RNA-Seq libraries:"
	echo -e '#    -m <file>     Path to tab-delimited file describing library configurations:'
	echo -e '#                  Column 1: Library name'
	echo -e '#                  Column 2: Library strandedness (= -r option)'
	echo -e '#                  Column 3: Path to R1 reads (= -1 option)'
	echo -e '#                  Column 4: Path to R2 reads (= -2 option)'
	echo -e '#                  Column 5: "trim" to quality-trim reads (= -q option)'
	echo -e '#                  Column 6: Path to read alignments (= -b option)'
	echo -e "#\n# SL quantification:"
	echo -e '#    -s <file>     Path to SL sequences (FASTA)'
	echo -e "#    -n <num>      Minimum 3' SL overlap with 5' read end (default: 8)"
	echo -e "#    -e <num>      Maximum error rate for nucleotide matching (default: 0.09)"
	echo -e "#    -f <chr>      Feature ID for read quantification (default: exon)"
	echo -e "#    -F <chr>      Meta-feature ID for read quantification (default: gene_id)"
	echo -e "#\n# Operon prediction:"
	echo -e "#    -z <chr>      [geo|sum|median] Method for aggregating SL counts across libraries (default: geometric mean [geo])"
	echo -e "#    -0            Keep libraries with zero counts when aggregating SL counts (default: remove zeros)"
	echo -e "#    -S <file>     Path to list of SL names (from -s) that resolve polycistrons (SL2-type; default: all SLs)"
	echo -e "#    -d <num>      Minimum SL2:SL1 read ratio required for downstream operonic genes (default: infinity, i.e. no SL1 reads)"
	echo -e "#    -u            Enforce SL2-bias (from -d) at upstream operonic genes"
	echo -e "#    -i <num>      Maximum intercistronic distance (default: infinity; x = infer)"
	echo -e "#    -x <file>     Reference operon annotations (GFF/GTF)"
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
outdir="SLOPPR_$(date +"%Y-%m-%d_%H-%M-%S")"
threads=8 # $(nproc)
stranded="x"
slength=8
err=0.09
clean="no"
opp="OP"
agg="geo"
zero="remove"
slrr="infinity"
upstream="no"
dist="infinity"
metafeatureid="gene_id"
featureid="exon"
# HPC control
hpc="no"
plocal=""
hpc_submit=""
hpc_depend=""
hpc_sep=":"
jobdeps=""

# parse command-line options
cmdline="sloppr.sh"
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
		-s)		sls=$2
				cmdline+=" -s $2"
		shift; shift;;
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
		-n)		slength=$2
				cmdline+=" -n $2"
		shift; shift ;;
		-e)		err=$2
				cmdline+=" -e $2"
		shift; shift ;;
		-F)		metafeatureid=$2
				cmdline+=" -F $2"
		shift; shift ;;
		-f)		featureid=$2
				cmdline+=" -f $2"
		shift; shift ;;
		-S)		sl2s=$2
				cmdline+=" -S $2"
		shift; shift;;
		-d)		slrr=$2
				cmdline+=" -d $2"
		shift; shift;;
		-u)		upstream="yes"
				cmdline+=" -u"
		shift;;
		-i)		dist=$2
				cmdline+=" -i $2"
		shift; shift;;
		-p)		opp=$2
				cmdline+=" -p $2"
		shift; shift;;
		-z)		agg=$2
				cmdline+=" -z $2"
		shift; shift;;
		-0)		zero="keep"
				cmdline+=" -0"
		shift;;
		-x)		refops=$2
				cmdline+=" -x $2"
		shift; shift;;
		-o)		outdir=$2
				cmdline+=" -o $2"
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
	echo -e "#   > Output directory\t\t $outdir"
	echo -e "#   > Threads\t\t\t $threads"
	echo -e "#   > Operon prefix\t\t $opp"
	echo -e "#   > Reference operons\t\t $refops"
	echo -e "#   > TEMP directory\t\t $TMPDIR"
	if [ "$hpc" == "yes" ]; then
		echo -e "#   > HPC parallel\t\t $hpc_submit"
		echo -e "#   > HPC sequential\t\t $hpc_depend"
	fi
	echo -e "#\n# Reference:"
	echo -e "#   > Genome\t\t\t $genome"
	echo -e "#   > Annotations\t\t $ann"
	echo -e "#   > Transcriptome\t $transcriptome"
	echo -e "#\n# SL quantification:"
	echo -e "#   > SL sequences\t\t $sls"
	echo -e "#   > Minimum SL tail\t\t $slength"
	echo -e "#   > Error rate\t\t $err"
	echo -e "#   > Feature ID\t\t $featureid"
	echo -e "#   > Meta-feature ID\t\t $metafeatureid"
	echo -e "#\n# Operon prediction:"
	echo -e "#   > SL-count aggregation\t $agg"
	echo -e "#   > Zero SL counts\t\t $zero"
	echo -e "#   > SL2-type SLs\t\t $sl2s"
	echo -e "#   > Downstream SL2:SL1 ratio\t >= $slrr"
	echo -e "#   > Upstream SL2-bias \t $upstream"
	echo -e "#   > Intercistronic distance\t <= $dist"
	echo -e "#\n# RNA-Seq data"
	if [ ! "$design" == "" ] ; then
		echo -e  "#   > Libraries\t\t $design"
	else
		echo -e  "#   > R1 reads\t\t\t $R1"
		echo -e  "#   > R2 reads\t\t\t $R2"
		echo -e  "#   > Quality trimming\t\t $clean"
		echo -e  "#   > BAM alignments\t\t $bam"
		echo -e  "#   > Strandedness\t\t $stranded"
	fi
}
print_summary

# sanity check input
function sanity_check {
	# we need to resolve "~" in file paths supplied via text file
	bam="${bam/\~/$HOME}"
	R1="${R1/\~/$HOME}"
	R2="${R2/\~/$HOME}"
	
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
	#if [ "$ann" == "" ]; then
	#	echo "ERROR: You must specify a GFF/GTF annotation file! Exiting."
	#	exit 1
	#elif
	if [ ! -f $ann ]; then
		echo "ERROR: The annotation file $ann does not exist! Exiting."
		exit 1
	fi
	if [ "$sls" == "" ]; then
		echo "ERROR: You must specify a FASTA file with spliced leaders! Exiting."
		exit 1
	elif [ ! -f $sls ]; then
		echo "ERROR: The spliced-leader file $sls does not exist! Exiting."
		exit 1
	fi
	if [ ! "$sl2s" == "" ] && [ ! -f $sl2s ]; then
		echo "ERROR: The SL2-type SL file $sl2s does not exist! Exiting."
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
	if [ ! "$refops" == "" ] && [ ! -f $refops ]; then
		echo "ERROR: The reference operon GFF/GTF file $refops does not exist! Exiting."
		exit 1
	fi
	if [ ! "$bam" == "" ] && ([ ! "$R1" == "" ] || [ ! "$R2" == "" ]); then
		echo "WARNING: You have specified a BAM file and reads. The reads will be ignored."
	fi
}
# first sweep of checks
sanity_check

# if library config file was supplied, sanity check that, too:
if [ ! "$design" == "" ] && ([ ! "$bam" == "" ] || [ ! "$R1" == "" ] || [ ! "$R2" == "" ]); then
	echo "WARNING: You have specified a library configuration file. Specified BAM or reads files will be ignored."
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
# 1) PREPARE GENOME AND GTF
#
echo "$(timestamp) >>> STAGE 1: SL read identification"

# copy SL file to output directory
cp $sls "$outdir/SLs.fa"

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
		echo "$(timestamp) Processing annotations ..."
		gropts="--force-exons --gene2exon --keep-genes -M -K -Q -C -T"
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
fi

# non-redundant exons
if [ ! "$ann" == "" ] && [ ! -f $outdir/annotations.unique_exons.gtf ]; then
	echo "$(timestamp) Extracting unique ${featureid}s ..."
	grep "\s$featureid\s" $ann \
		| bedtools sort -i stdin \
		| bedtools merge -s -c 4,5,7,9 -o min,max,first,first -i stdin \
		| awk -F'\t' -v f="$featureid" '{print $1, ".", f, $4, $5, ".", $6, ".", $7}' OFS="\t" \
		> $outdir/annotations.unique_exons.gtf
fi

# make hisat2 index for genome
if [ ! "$genome" == "" ] && [ ! -d "$outdir/hisat2_index" ]; then
	echo "$(timestamp) Generating HISAT2 index ..."
	mkdir -p $outdir/hisat2_index
	# extract splice sites for hisat2
	hisat2_extract_splice_sites.py $ann > $outdir/hisat2_index/hisat2_splicesites.txt
	hisat2_extract_exons.py $ann > $outdir/hisat2_index/hisat2_exons.txt
	if [ -s $outdir/hisat2_index/hisat2_splicesites.txt ] && [ -s $outdir/hisat2_index/hisat2_exons.txt ]; then
		hisatgffopts+="--ss $outdir/hisat2_index/hisat2_splicesites.txt --exon $outdir/hisat2_index/hisat2_exons.txt"
	fi
	hisat2-build -p $threads $hisatgffopts $ref $outdir/hisat2_index/genome > $outdir/hisat2_index/log_hisat2-build.txt 2>&1
fi

# make bowtie2 index for transcriptome
if [ ! "$transcriptome" == "" ] && [ ! -d "$outdir/bowtie2_index" ]; then
	echo "$(timestamp) Generating BOWTIE2 index ..."
	mkdir -p "$outdir/bowtie2_index"
	bowtie2-build --threads $threads \
		"$ref" "$outdir/bowtie2_index/transcriptome" \
		> "$outdir/bowtie2_index/log_bowtie2-build.txt" 2>&1
fi

#
# 2) GENOME ALIGNMENT and SL READ IDENTIFICATION
#
#
# if a library-design file is specified, loop through libraries
# for each library, check what needs done
# new: check whether library.OK file exists. If not, submit job to check what needs done

function library_pipeline {
	echo "$(timestamp) >>> Processing library $lib"
	lib="1-library_$lib"
		
	# does library.OK file exist? If yes, library has already been processed
	lok="$outdir/$lib/library-n$slength-e$err-f$featureid-F$metafeatureid.OK"
	if [ -f "$lok" ]; then 
		echo "$(timestamp) >>> Nothing to do"
	else
		# process via HPC or local
		if [ "$hpc" == "yes" ]; then
			# HPC
			jobID=$(eval "$hpc_submit sloppr_library.sh -var=\"$threads\" -var=\"$outdir\" -var=\"$lib\" -var=\"$stranded\" -var=\"$R1\" -var=\"$R2\" -var=\"$clean\" -var=\"$bam\" -var=\"$slength\" -var=\"$err\" -var=\"$featureid\" -var=\"$metafeatureid\"")
			echo "$(timestamp) >>> Submitted HPC job $jobID"
			jobdeps+="$jobID$hpc_sep"
		else
			# local
			sloppr_library.sh -var="$threads" -var="$outdir" -var="$lib" -var="$stranded" -var="$R1" -var="$R2" -var="$clean" -var="$bam" -var="$slength" -var="$err" -var="$featureid" -var="$metafeatureid"
		fi
	fi
}

if [ ! "$design" == "" ]; then
	while IFS='#' read lib stranded R1 R2 clean bam
	do
		library_pipeline
	done < <(tr -d '#' < $design | tr '\t' '#' )
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





# generate merged BAM files for convenient viewing in IGV or similar
echo "$(timestamp) Merging SL BAM files ..."
mergedir=$outdir/1-merged_SL_BAM
mkdir -p $mergedir
# for each SL
for SL in $(grep "^>" $sls | tr -d ">")
do
	if [ ! -f $mergedir/$SL-n$slength-e$err.all-libraries.bam ]; then
		samtools merge $mergedir/$SL-n$slength-e$err.all-libraries.bam $outdir/1-library_*/$SL-n$slength-e$err.fc.bam
		samtools flagstat $mergedir/$SL-n$slength-e$err.all-libraries.bam > $mergedir/$SL-n$slength-e$err.all-libraries.bam.flagstat
		samtools index $mergedir/$SL-n$slength-e$err.all-libraries.bam
	fi
done
# for each SL type (if unknown, merge all)
if [ ! "$sl2s" == "" ] && [ ! -f $mergedir/allSL2type-n$slength-e$err.all-libraries.bam ]; then
	# SL1-type
	find $outdir/1-library_* -name "*-n$slength-e$err.fc.bam" \
		| grep -Fw -f <(grep "^>" $sls | tr -d ">" | grep -Fwv -f <(tr -d ">" < $sl2s)) \
		| xargs samtools merge $mergedir/allSL1type-n$slength-e$err.all-libraries.bam
	samtools flagstat $mergedir/allSL1type-n$slength-e$err.all-libraries.bam > $mergedir/allSL1type-n$slength-e$err.all-libraries.bam.flagstat
	samtools index $mergedir/allSL1type-n$slength-e$err.all-libraries.bam
	# SL2-type
	find $outdir/1-library_* -name "*-n$slength-e$err.fc.bam" \
		| grep -Fw -f <(grep "^>" $sls | tr -d ">" | grep -Fw -f <(tr -d ">" < $sl2s)) \
		| xargs samtools merge $mergedir/allSL2type-n$slength-e$err.all-libraries.bam
	samtools flagstat $mergedir/allSL2type-n$slength-e$err.all-libraries.bam > $mergedir/allSL2type-n$slength-e$err.all-libraries.bam.flagstat
	samtools index $mergedir/allSL2type-n$slength-e$err.all-libraries.bam
elif [ "$sl2s" == "" ] && [ ! -f $mergedir/allSLs-n$slength-e$err.all-libraries.bam ]; then
	find $outdir/1-library_* -name "*-n$slength-e$err.fc.bam" \
		| xargs samtools merge $mergedir/allSLs-n$slength-e$err.all-libraries.bam
	samtools flagstat $mergedir/allSLs-n$slength-e$err.all-libraries.bam > $mergedir/allSLs-n$slength-e$err.all-libraries.bam.flagstat
	samtools index $mergedir/allSLs-n$slength-e$err.all-libraries.bam
fi

#
# 3) READ QUANTIFICATION
#
echo "$(timestamp) >>> STAGE 2: SL quantification"
ctdir=$outdir/2-counts

# SL counts
mkdir -p $ctdir
if [ ! -f $ctdir/SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.genes.raw.txt ]; then
	echo "$(timestamp) Quantifying SL reads against genes ..."
	featureCounts -a $ann \
		-o $ctdir/SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.genes.raw.txt \
		-t $featureid -g $metafeatureid -s 1 -O -M -T $threads \
		$outdir/1-lib*/*-n$slength-e$err.fc.bam > $ctdir/log_featureCounts.SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.genes.txt 2>&1
fi
if [ ! -f $ctdir/SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.exons.raw.txt ]; then
	echo "$(timestamp) Quantifying SL reads against ${featureid}s ..."	
	featureCounts -a $outdir/annotations.unique_exons.gtf -o $ctdir/SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.exons.raw.txt \
		-t $featureid -g $metafeatureid -f -s 1 -O -M -T $threads \
		$outdir/1-lib*/*-n$slength-e$err.fc.bam > $ctdir/log_featureCounts.SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.exons.txt 2>&1
fi

# summarise background counts
if [ ! -f $ctdir/bg.counts-f$featureid-F$metafeatureid.genes.raw.txt ]; then
	# paste 7th column for all background count files
	fcbg=($outdir/*/end2end_pre-align.bam.fc-f$featureid-F$metafeatureid.txt)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(tail -n +2 $a | cut -f 7)"
		done | paste -s -d ' ')
	eval "paste <(tail -n +2 $fcbg | cut -f 1-6) $x" > $ctdir/bg.counts-f$featureid-F$metafeatureid.genes.raw.txt
	# also grab all summary files
	fcbg=($outdir/*/end2end_pre-align.bam.fc-f$featureid-F$metafeatureid.txt.summary)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(cut -f 2 $a)"
		done | paste -s -d ' ')
	eval "paste <(cut -f 1 $fcbg) $x" > $ctdir/bg.counts-f$featureid-F$metafeatureid.genes.raw.txt.summary
fi
# summarise rescued read counts
if [ ! -f $ctdir/un.featureCounts.genes.raw.txt ]; then
	# paste 7th column for all background count files
	fcbg=($outdir/*/rescued-n$slength-e$err.bam.fc-f$featureid-F$metafeatureid.txt)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(tail -n +2 $a | cut -f 7)"
		done | paste -s -d ' ')
	eval "paste <(tail -n +2 $fcbg | cut -f 1-6) $x" > $ctdir/rescued-n$slength-e$err.counts-f$featureid-F$metafeatureid.genes.raw.txt
	# also grab all summary files
	fcbg=($outdir/*/rescued-n$slength-e$err.bam.fc-f$featureid-F$metafeatureid.txt.summary)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(cut -f 2 $a)"
		done | paste -s -d ' ')
	eval "paste <(cut -f 1 $fcbg) $x" > $ctdir/rescued-n$slength-e$err.counts-f$featureid-F$metafeatureid.genes.raw.txt.summary
fi

#
# 4) CURATE GENE ANNOTATIONS
# Correct gene annotations using SL read peaks at internal exons
# Arguments:	1) directory containing featureCounts output files
#				2) number of threads
#				3) minimum number of reads required for a peak

# read length
rlen=$(samtools cat $outdir/*/*.fc.bam | samtools stats | grep ^RL | tail -n 1 | cut -f 2)

if [ ! -f $ctdir/SL.featureCounts.genes.clean.txt ] || \
   [ ! -f $ctdir/SL.featureCounts.exons.clean.txt ] || \
   [ ! -f $ctdir/bg.featureCounts.genes.clean.txt ] || \
   [ ! -f $ctdir/un.featureCounts.genes.clean.txt ]; then
	echo "$(timestamp) Curating gene annotations ..."	
	sloppr_curate_counts.R $ctdir/SL-n$slength-e$err.counts-f$featureid-F$metafeatureid $threads $rlen 4
fi

#
# 5) PREDICT OPERONS
# Arguments:	1) featureCounts output file
#				2) SL FASTA file
#				3) text file with SL2-type SL names
#				4) SL2:SL1 read-ratio threshold to classify gene as downstream (default: infinity)
#				5) Operon prefix (default: OP)
#				6) number of threads
echo "$(timestamp) >>> STAGE 3: Operon prediction"

opdir="$outdir/3-operons-$(date +"%Y-%m-%d-%H:%M:%S")-z$agg-0$zero-d$slrr-u$upstream-i$dist"
mkdir -p $opdir
#if [ -f $outdir/operons/operons.gff3 ]; then
	echo "$(timestamp) Predicting operons ..."	

	sloppr_predict_operons.R $ctdir/SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.genes.clean.txt $sls " $sl2s" $slrr $upstream $dist $opp "$agg" $zero $threads $opdir
	sloppr_predict_operons.R $ctdir/SL-n$slength-e$err.counts-f$featureid-F$metafeatureid.exons.clean.txt $sls " $sl2s" $slrr $upstream $dist $opp "$agg" $zero $threads $opdir
#fi

#
# 6) COMPARE AGAINST REFERENCE OPERONS
#
if [ ! "$refops" == "" ]; then 
	#grep "operon" $ann > $outdir/annotations.operons.gff3
	echo "$(timestamp) Comparing predictions against reference operons ..."
	for gff in $opdir/*.gff3
	do
		echo "Comparing reference operons against $gff ..."
		outf=$gff.vs.refoperons.txt
		nops=$(grep -i -c "\soperon\s" $gff)
		nref=$(wc -l < $refops)
		vops=$(bedtools intersect -u -s -a <(grep -i "\soperon\s" $gff) -b $refops | wc -l)
		echo -e "\tPredicted $nops operons, of which $vops ($((100*$vops/$nops))%) are among $nref reference operons ($((100*$vops/$nref))%)" > $outf

		nopg=$(grep -c "\sgene\s" $gff)
		vopg=$(bedtools intersect -u -s -a <(grep "\sgene\s" $gff) -b $refops | wc -l)
		echo -e "\tPredicted $nopg operonic genes, of which $vopg ($((100*$vopg/$nopg))%) match reference operons" >> $outf
		cat $outf
	done
fi

echo "$(timestamp) Finished!"
exit 0