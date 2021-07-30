#!/usr/bin/env bash
set -e

function timestamp { 
	date +"[%Y-%m-%d %H:%M:%S]"
}

function load_dependencies {
	# ensure that all dependencies are in PATH
	# MODIFY/DELETE THESE CALLS ACCORDING TO YOUR MACHINE
	echo -e "#\n$(timestamp) >>> Loading dependencies"
	module load hisat2-2.1.0
	module load samtools-1.9
	module load bedtools-2.28.0
	module load seqtk-1.3
	module load r-3.6.0
	module load gffread-0.11.4
	module load cutadapt-2.3
	module load subread-1.6.2
}

title="#\n# SLOPPR - Spliced-Leader-informed Operon Prediction from RNA-Seq data\n# Version 1.1.4\n#"

function printhelp {
	echo -e "$title"
	echo -e "# Usage: sloppr.sh [-o <>] [-c <>] -g <> -a <>"
	echo -e "# \t\t   {-1 <> [-2 <>] | -b <> | -m <>} [-r 0|1|2|x] [-q]"
	echo -e "# \t\t   -s <> [-n <>] [-e <>] [f <>]"
	echo -e "# \t\t   [-z geo|sum|median] [-0] [-S <>] [-d <>] [-u] [-i <>] [-p <>] [-x <>]"
	echo -e "#\n# General options:"
	echo -e '#    -o <dir>\tPath to output directory (default: "SLOPPR_[date+time]")'
	echo -e "#    -p <chr>\tPrefix for predicted operons (default: OP)"
	echo -e '#    -c <num>\tThreads (default: all available cores)'
	echo -e '#    -T <num>\tTEMP directory (default: $TMPDIR)'
	echo -e '#    -h\t\tPrint this help page'
	echo -e "#\n# Reference assembly:"
	echo -e '#    -g <file>\tPath to genome assembly (FASTA[.gz])'
	echo -e '#    -a <file>\tPath to genome annotations (GFF/GTF[.gz])'
	echo -e "#\n# Single RNA-Seq library:"
	echo -e '#    -1 <file>\tPath to R1 reads (FASTQ[.gz])'
	echo -e '#    -2 <file>\tPath to R2 reads (FASTQ[.gz])'
	echo -e '#    -b <file>\tPath to read alignments (BAM)'
	echo -e '#    -r <num>\tRead strandedness (0=unstranded; 1=stranded; 2=reverse-stranded; x=infer; default: x)'
	echo -e "#    -q\t\tQuality-trim 3' ends of reads and remove adapters"
	echo -e "#\n# Multiple RNA-Seq libraries:"
	echo -e '#    -m <file>\tPath to tab-delimited file describing library configurations:'
	echo -e '#    \t\tColumn 1: Library name'
	echo -e '#    \t\tColumn 2: Library strandedness (= -r option)'
	echo -e '#    \t\tColumn 3: Path to R1 reads (= -1 option)'
	echo -e '#    \t\tColumn 4: Path to R2 reads (= -2 option)'
	echo -e '#    \t\tColumn 5: "trim" to quality-trim reads (= -q option)'
	echo -e '#    \t\tColumn 6: Path to read alignments (= -b option)'
	echo -e "#\n# SL quantification:"
	echo -e '#    -s <file>\tPath to SL sequences (FASTA)'
	echo -e "#    -n <num>\tMinimum 3' SL overlap with 5' read end (default: 8)"
	echo -e "#    -e <num>\tMaximum error rate for nucleotide matching (default: 0.09)"
	echo -e "#    -f <chr>\tFeature ID for read quantification (default: exon)"
	echo -e "#    -F <chr>\tMeta-feature ID for read quantification (default: gene_id)"
	echo -e "#\n# Operon prediction:"
	echo -e "#    -z <chr>\t[geo|sum|median] Method for aggregating SL counts across libraries (default: geometric mean [geo])"
	echo -e "#    -0\t\tKeep libraries with zero counts when aggregating SL counts (default: remove zeros)"
	echo -e "#    -S <file>\tPath to list of SL names (from -s) that resolve polycistrons (SL2-type; default: all SLs)"
	echo -e "#    -d <num>\tMinimum SL2:SL1 read ratio required for downstream operonic genes (default: infinity, i.e. no SL1 reads)"
	echo -e "#    -u\t\tEnforce SL2-bias (from -d) at upstream operonic genes"
	echo -e "#    -i <num>\tMaximum intercistronic distance (default: infinity; x = infer)"
	echo -e "#    -x <file>\tReference operon annotations (GFF/GTF)"
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
threads=$(nproc)
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
		-T)		export TMPDIR=$2
				cmdline+=" -T $2"
		shift; shift;;
	    -)	            # unknown option
		shift; shift;;
	esac
done

# print input summary on screen and log file
function print_summary {
	echo -e "$title"
	echo -e "# General:"
	echo -e "#   > Output directory\t\t $outdir"
	echo -e "#   > Threads\t\t\t $threads"
	echo -e "#   > Operon prefix\t\t $opp"
	echo -e "#   > Reference operons\t\t $refops"
	echo -e "#   > TEMP directory\t\t $TMPDIR"
	echo -e "#\n# Reference:"
	echo -e "#   > Genome\t\t\t $genome"
	echo -e "#   > Annotations\t\t $ann"
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
	
	if [ "$genome" == "" ]; then
		echo "ERROR: You must specify a genome! Exiting."
		exit 1
	elif [ ! -f $genome ]; then
		echo "ERROR: The genome $genome does not exist! Exiting."
		exit 1
	fi
	if [ "$ann" == "" ]; then
		echo "ERROR: You must specify a GFF/GTF annotation file! Exiting."
		exit 1
	elif [ ! -f $ann ]; then
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

# if genome is zipped, copy and unzip
if [ "$(file $genome | grep -i gzip)" ]; then
	if [ ! -f $outdir/genome.fa ]; then
		gunzip -c $genome > $outdir/genome.fa
	fi
	genome=$outdir/genome.fa
fi

# copy annotations
# unzip if required
# convert GFF to GTF if required
gtf="$outdir/annotations.gtf"
if [ ! -f $gtf ]; then
	echo "$(timestamp) Processing annotations ..."
	gropts="--force-exons --gene2exon --keep-genes -M -K -Q -C -T"
	fform="$(file $ann)"
	if [[ "$fform" =~ "gzip" ]] && [[ "$fform" =~ ".gtf" ]]; then
		# gzipped GTF
		gunzip -c $ann > $gtf
	elif [[ "$fform" =~ "gzip" ]] && [[ ! "$fform" =~ ".gtf" ]]; then
		# gzipped GFF
		gunzip -c $ann | gffread $gropts -o $gtf
	elif [[ ! "$fform" =~ "gzip" ]] && [[ "$fform" =~ ".gtf" ]]; then
		# unzipped GTF
		cp $ann $gtf
	else
		# unzipped GFF
		gffread $ann $gropts -o $gtf
	fi
fi
ann=$gtf

# non-redundant exons
if [ ! -f $outdir/annotations.unique_exons.gtf ]; then
	echo "$(timestamp) Extracting unique ${featureid}s ..."
	grep "\s$featureid\s" $ann \
		| bedtools sort -i stdin \
		| bedtools merge -s -c 4,5,7,9 -o min,max,first,first -i stdin \
		| awk -F'\t' -v f="$featureid" '{print $1, ".", f, $4, $5, ".", $6, ".", $7}' OFS="\t" \
		> $outdir/annotations.unique_exons.gtf
fi

# make hisat2 index 
if [ ! -d $outdir/hisat2_index ]; then
	echo "$(timestamp) Generating HISAT2 index ..."
	mkdir -p $outdir/hisat2_index
	# extract splice sites for hisat2
	hisat2_extract_splice_sites.py $ann > $outdir/hisat2_index/hisat2_splicesites.txt
	hisat2_extract_exons.py $ann > $outdir/hisat2_index/hisat2_exons.txt
	if [ -s $outdir/hisat2_index/hisat2_splicesites.txt ] && [ -s $outdir/hisat2_index/hisat2_exons.txt ]; then
		hisatgffopts+="--ss $outdir/hisat2_index/hisat2_splicesites.txt --exon $outdir/hisat2_index/hisat2_exons.txt"
	fi
	hisat2-build -p $threads $hisatgffopts $genome $outdir/hisat2_index/genome > $outdir/hisat2_index/log_hisat2-build.txt 2>&1
fi

#
# 2) GENOME ALIGNMENT and SL READ IDENTIFICATION
#

function trim_reads {
	trR1=$outdir/$lib/${lib/1-library_/}.R1.trimmed.fq.gz
	trR2=$outdir/$lib/${lib/1-library_/}.R2.trimmed.fq.gz
	if [ "$clean" == "trim" ] && [ ! -f $trR1 ]; then
		echo "$(timestamp) Trimming adapters and poor-quality bases from reads ..."
		if [ "$R2" == "" ]; then
			# single-end data
			raw=$R1
			R1=$trR1
			cutadapt -a AGATCGGAAGAGC -q 20 -m 20 \
				-o $R1 -j $threads $raw \
				> $outdir/$lib/log_cutadapt.txt 2>&1
		else
			# paired-end data
			raw1=$R1
			raw2=$R2
			R1=$trR1
			R2=$trR2
			cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -m 20 \
				-o $R1 -p $R2 -j $threads $raw1 $raw2 \
				> $outdir/$lib/log_cutadapt.txt 2>&1
		fi
	fi
}

function genome_align {
	if [ ! -f $outbam ]; then
		echo "$(timestamp) Aligning reads to genome ..."
		hisat_options=""
		if [ "$R2" == "" ]; then
			hisat_options+="-U $R1"
		else
			hisat_options+="-1 $R1 -2 $R2"
		fi
		hisat2 --no-softclip --no-discordant --dta \
			-x $outdir/hisat2_index/genome $hisat_options -p $threads 2> $outdir/$lib/log_hisat2.txt \
			| samtools sort -O bam -@ $threads -o $outbam - 2>> $outdir/$lib/log_hisat2.txt
		samtools flagstat $outbam > $outbam.flagstat
		samtools index -@ $threads -b $outbam
		echo "$(timestamp) $(grep 'overall alignment rate' $outdir/$lib/log_hisat2.txt)"
	fi
}

function infer_strandedness {
	if [ -f $outbam.strandedness.txt ]; then
		stranded=$(cat $outbam.strandedness.txt)
	else
		# check BAM to infer whether we have SE or PE reads
		# if SE, must treat as unstranded data
		# if PE, infer strandedness if requested
		npe=$({ samtools view -H $bam ; samtools view $bam | head -n 10000; } | samtools view -c -f 1 -)	
		if [ "$npe" == "0" ]; then
			stranded=0
		elif [[ ! "$stranded" =~ [012] ]]; then
			# infer strandedness
			# + strand
			gtfsample=$(awk '$3~"transcript|gene|mRNA" && $7=="+"' $ann \
				| bedtools sort -i - | bedtools merge -s -i - | bedtools sample -n 1000 -i -)
			pR1=$(echo "$gtfsample" | samtools view -c -F 20 -f 64 -M -L - $bam)
			pR2=$(echo "$gtfsample" | samtools view -c -F 20 -f 128 -M -L - $bam)
				# - strand
			gtfsample=$(awk '$3~"transcript|gene|mRNA" && $7=="-"' $ann \
				| bedtools sort -i - | bedtools merge -s -i - | bedtools sample -n 1000 -i -)
			mR1=$(echo "$gtfsample" | samtools view -c -F 4 -f 80 -M -L - $bam)
			mR2=$(echo "$gtfsample" | samtools view -c -F 4 -f 144 -M -L - $bam)
				# summarise
			R1=$((pR1+mR1))
			R2=$((pR2+mR2))
			if [ $((R1+R2)) == 0 ]; then
				# double-check for single-end data; must be analysed as unstranded
				stranded=0
			else
				# need at least 30:70 bias to infer as stranded
				sR1=$((100*R1/(R1+R2)))
				sR2=$((100*R2/(R1+R2)))
				stranded=$(( 1 + 2*(sR1<=30) - (sR1<=70)))
			fi
			echo "$stranded" > "$outbam.strandedness.txt"	
		fi
		echo "$(timestamp) Inferred library strandedness: $stranded"
	fi
}

function read_extract {
	# For PE data, we want unmapped reads with mapped mates. The unmapped mate is the SL candidate
	# only use primary alignments 
	# rev-stranded data (-r 2): R2 read unmapped, R1 mate mapped	-f 132	-F 264		-f 72	-F 260	
	# fwd-stranded data (-r 1): R1 read unmapped, R2 mate mapped 	-f 68	-F 264		-f 136	-F 260
	# unstranded data (-r 0): No distinction between R1 and R2		-f 4	-F 264		-f 8	-F 260
	# For SE data, we simply want unmapped reads (-f 4)

	if [ ! -f $outdir/$lib/SL.candidates.fastq.gz ]; then
		echo "$(timestamp) Extracting candidate reads ..."
		if [ "$npe" == "0" ]; then		# SE reads
			# unmapped SE reads
			samtools view -h -f 4 $bam \
				| samtools sort -n - \
				| samtools fastq -n - 2> $outdir/$lib/log_samtoolsfastq.txt | gzip > $outdir/$lib/SL.candidates.fastq.gz
		else	# PE reads
			# unmapped reads
			samtools view -h -f $((4+ (64*$stranded))) -F 264 $bam \
				| samtools sort -n - \
				| samtools fastq -n -s $outdir/$lib/SL.candidates.fastq.gz - \
				> $outdir/$lib/log_samtoolsfastq.txt 2>&1
			# mapped mates
			samtools view -h -f $((8+ (192-64*$stranded)*($stranded>0))) -F 260 $bam \
				| samtools sort -n - \
				| samtools fastq -n -s $outdir/$lib/SL.mates.fastq.gz - \
				> $outdir/$lib/log_samtoolsfastq_mates.txt 2>&1
		fi
		nr=$(grep -o "[0-9]* reads" $outdir/$lib/log_samtoolsfastq.txt | grep -o "[0-9]*")
		echo "$(timestamp) Found $nr candidate read(pair)s ..."
	fi
}

function sl_screen {
	if [ $(find -L $outdir/$lib -name "*5p.fastq.gz" | wc -l) == 0 ]; then
		echo "$(timestamp) Screening candidate reads for SLs ..."
		if [ "$npe" == "0" ]; then
			# SE
			# We also need to check the rev complement of the reads (unstranded data)
			{ seqtk seq -r $outdir/$lib/SL.candidates.fastq.gz; 
			  seqtk seq $outdir/$lib/SL.candidates.fastq.gz; } \
				| cutadapt -g file:$sls \
				-O $slength -e $err -m 15 \
				--untrimmed-output $outdir/$lib/unsuccessful.5p.fastq.gz \
				-o $outdir/$lib/{name}.5p.fastq.gz - > $outdir/$lib/log_SL.cutadapt.txt 2>&1
		else
			# PE
			cutadapt -g file:$sls \
				-O $slength -e $err -m 15 \
				--untrimmed-output $outdir/$lib/unsuccessful.5p.fastq.gz \
				--untrimmed-paired-output $outdir/$lib/unsuccessful.mates.fastq.gz \
				-o $outdir/$lib/{name}.5p.fastq.gz \
				-p $outdir/$lib/{name}.mates.fastq.gz \
				$outdir/$lib/SL.candidates.fastq.gz \
				$outdir/$lib/SL.mates.fastq.gz > $outdir/$lib/log_SL.cutadapt.txt 2>&1
		fi
		nadapt=$(grep "with adapter" $outdir/$lib/log_SL.cutadapt.txt | head -n 1 | sed 's/.*  //g')
		echo "$(timestamp) $nadapt SL tails found among candidate reads"
	fi
}

function sl_realign {
	# align SL reads (and unsuccessful reads) back to the genome
	for SL in $(grep "^>" $sls | tr -d ">") unsuccessful
	do
		if [ ! -f $outdir/$lib/$SL.*bam ]; then
			# SE or PE data
			hisat_options=""
			if [ "$npe" == "0" ]; then
				hisat_options+="-U $outdir/$lib/$SL.5p.fastq.gz"
			else
				hisat_options+="-1 $outdir/$lib/$SL.5p.fastq.gz -2 $outdir/$lib/$SL.mates.fastq.gz"
			fi
			# unsuccessful or SL
			if [ "$SL" == "unsuccessful" ]; then
				# unsuccessful reads are allowed softclipping
				# keep mates if we have PE data
				suff=""
				echo "$(timestamp) Aligning $SL reads ..."
				hisat2 -x $outdir/hisat2_index/genome $hisat_options -p $threads 2> $outdir/$lib/log_$SL.hisat2.txt \
					| samtools sort -O bam -@ $threads -o $outdir/$lib/$SL$suff.bam - 2>> $outdir/$lib/log_$SL.hisat2.txt
			else
				# trimmed SL reads must align end-to-end
				# align mates (if available) to aid alignment,
				# but remove them from the BAM file afterwards (we want to quantify the 5' end only)
				suff=".fc"
				nr=" $(grep -A 2 "$SL" $outdir/$lib/log_SL.cutadapt.txt | tail -n 1 | grep -o "[0-9]* times" | grep -o "[0-9]*")"
				echo "$(timestamp) Aligning$nr $SL reads ..."
				hisat2 -x $outdir/hisat2_index/genome --no-softclip $hisat_options -p $threads 2> $outdir/$lib/log_$SL.hisat2.txt \
					| samtools view -@ $threads -h -f $((64*($npe>0))) - 2>> $outdir/$lib/log_$SL.hisat2.txt \
					| samtools sort -@ $threads -n - 2>> $outdir/$lib/log_$SL.hisat2.txt \
					| samtools fixmate -@ $threads - - 2>> $outdir/$lib/log_$SL.hisat2.txt \
					| samtools sort -O bam -@ $threads -o $outdir/$lib/$SL$suff.bam - 2>> $outdir/$lib/log_$SL.hisat2.txt
			fi
			samtools flagstat $outdir/$lib/$SL$suff.bam > $outdir/$lib/$SL$suff.bam.flagstat 
			samtools index -@ $threads -b $outdir/$lib/$SL$suff.bam
		fi
	done
}


function quantify_background {
	# PE data needs -p flag in featureCounts
	if [ "$npe" == "0" ]; then
		fcpaired=""
	else
		fcpaired="-p"
	fi

	# quantify end2end alignments
	# use actual strandedness
	if [ ! -f $outbam.fc.txt ]; then
		echo "$(timestamp) Quantifying background gene coverage ..."	
		if [ ! -f $outbam ]; then
			ln -fs $(pwd)/$bam $outbam
		fi
		featureCounts -a $ann -o $outbam.fc.txt \
			-t $featureid -g $metafeatureid -s $stranded $fcpaired -C -O -M -T $threads \
			$outbam > $outdir/$lib/log_featureCounts.bg.txt 2>&1
	fi
	
	# quantify unsuccessful candidate reads (not end2end and no SL)
	# due to the pseudo-stranding done during candidate read extraction, strandedness is 1 (or 0 for unstranded data)
	if [ ! -f $outdir/$lib/unsuccessful.bam.txt ]; then
		echo "$(timestamp) Quantifying unsuccessful reads ..."	
		featureCounts -a $ann -o $outdir/$lib/unsuccessful.bam.txt \
			-t $featureid -g $metafeatureid -s $(($stranded>0)) $fcpaired -C -O -M -T $threads \
			$outdir/$lib/unsuccessful.bam > $outdir/$lib/log_featureCounts.unsuccessful.txt 2>&1
	fi
}

#
# if a library-design file is specified, loop through libraries
function library_pipeline {
	echo "$(timestamp) >>> Processing library $lib"
	lib="1-library_$lib"
	outbam="$outdir/$lib/end2end_pre-align.bam"
	mkdir -p $outdir/$lib
	if [ "$bam" == "" ]; then
		trim_reads
		genome_align
		bam=$outbam
	fi
	infer_strandedness
	read_extract
	sl_screen
	sl_realign
	quantify_background	
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

# generate merged BAM files for convenient viewing in IGV or similar
echo "$(timestamp) Merging SL BAM files ..."
mergedir=$outdir/1-merged_SL_BAM
mkdir -p $mergedir
# for each SL
for SL in $(grep "^>" $sls | tr -d ">")
do
	if [ ! -f $mergedir/$SL.all-libraries.bam ]; then
		samtools merge $mergedir/$SL.all-libraries.bam $outdir/1-library_*/$SL.fc.bam
		samtools flagstat $mergedir/$SL.all-libraries.bam > $mergedir/$SL.all-libraries.bam.flagstat
		samtools index $mergedir/$SL.all-libraries.bam
	fi
done
# for each SL type (if unknown, merge all)
if [ ! "$sl2s" == "" ] && [ ! -f $mergedir/allSL2-type.all-libraries.bam ]; then
	# SL1-type
	find $outdir/1-library_* -name "*.bam" \
		| grep -Fw -f <(grep "^>" $sls | tr -d ">" | grep -Fwv -f <(tr -d ">" < $sl2s)) \
		| xargs samtools merge $mergedir/allSL1-type.all-libraries.bam
	samtools flagstat $mergedir/allSL1-type.all-libraries.bam > $mergedir/allSL1-type.all-libraries.bam.flagstat
	samtools index $mergedir/allSL1-type.all-libraries.bam
	# SL2-type
	find $outdir/1-library_* -name "*.bam" \
		| grep -Fw -f <(grep "^>" $sls | tr -d ">" | grep -Fw -f <(tr -d ">" < $sl2s)) \
		| xargs samtools merge $mergedir/allSL2-type.all-libraries.bam
	samtools flagstat $mergedir/allSL2-type.all-libraries.bam > $mergedir/allSL2-type.all-libraries.bam.flagstat
	samtools index $mergedir/allSL2-type.all-libraries.bam
elif [ "$sl2s" == "" ] && [ ! -f $mergedir/allSLs.all-libraries.bam ]; then
	find $outdir/1-library_* -name "*.fc.bam" \
		| xargs samtools merge $mergedir/allSLs.all-libraries.bam
	samtools flagstat $mergedir/allSLs.all-libraries.bam > $mergedir/allSLs.all-libraries.bam.flagstat
	samtools index $mergedir/allSLs.all-libraries.bam
fi

#
# 3) READ QUANTIFICATION
#
echo "$(timestamp) >>> STAGE 2: SL quantification"
ctdir=$outdir/2-counts

# SL counts
mkdir -p $ctdir
if [ ! -f $ctdir/SL.featureCounts.genes.raw.txt ]; then
	echo "$(timestamp) Quantifying SL reads against genes ..."
	featureCounts -a $ann -o $ctdir/SL.featureCounts.genes.raw.txt \
		-t $featureid -g $metafeatureid -s 1 -O -M -T $threads \
		$outdir/*/*.fc.bam > $ctdir/log_featureCounts.genes.txt 2>&1
fi
if [ ! -f $ctdir/SL.featureCounts.exons.raw.txt ]; then
	echo "$(timestamp) Quantifying SL reads against ${featureid}s ..."	
	featureCounts -a $outdir/annotations.unique_exons.gtf -o $ctdir/SL.featureCounts.exons.raw.txt \
		-t $featureid -g $metafeatureid -f -s 1 -O -M -T $threads \
		$outdir/*/*.fc.bam > $ctdir/log_featureCounts.exons.txt 2>&1
fi

# summarise background counts
if [ ! -f $ctdir/bg.featureCounts.genes.raw.txt ]; then
	# paste 7th column for all background count files
	fcbg=($outdir/*/end2end_pre-align.bam.fc.txt)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(tail -n +2 $a | cut -f 7)"
		done | paste -s -d ' ')
	eval "paste <(tail -n +2 $fcbg | cut -f 1-6) $x" > $ctdir/bg.featureCounts.genes.raw.txt
	# also grab all summary files
	fcbg=($outdir/*/end2end_pre-align.bam.fc.txt.summary)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(cut -f 2 $a)"
		done | paste -s -d ' ')
	eval "paste <(cut -f 1 $fcbg) $x" > $ctdir/bg.featureCounts.genes.raw.txt.summary
fi
if [ ! -f $ctdir/un.featureCounts.genes.raw.txt ]; then
	# paste 7th column for all background count files
	fcbg=($outdir/*/unsuccessful.bam.txt)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(tail -n +2 $a | cut -f 7)"
		done | paste -s -d ' ')
	eval "paste <(tail -n +2 $fcbg | cut -f 1-6) $x" > $ctdir/un.featureCounts.genes.raw.txt
	# also grab all summary files
	fcbg=($outdir/*/unsuccessful.bam.txt.summary)
	x=$(for a in ${fcbg[@]}
		do
			echo "<(cut -f 2 $a)"
		done | paste -s -d ' ')
	eval "paste <(cut -f 1 $fcbg) $x" > $ctdir/un.featureCounts.genes.raw.txt.summary
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
	sloppr_curate_counts.R $ctdir $threads $rlen 4
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

opdir="$outdir/3-operons-z$agg-0$zero-d$slrr-u$upstream-i$dist"
mkdir -p $opdir
#if [ -f $outdir/operons/operons.gff3 ]; then
	echo "$(timestamp) Predicting operons ..."	

	sloppr_predict_operons.R $ctdir/SL.featureCounts.genes.clean.txt $sls " $sl2s" $slrr $upstream $dist $opp "$agg" $zero $threads $opdir
	sloppr_predict_operons.R $ctdir/SL.featureCounts.exons.clean.txt $sls " $sl2s" $slrr $upstream $dist $opp "$agg" $zero $threads $opdir
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