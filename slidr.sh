#!/bin/bash

function timestamp { 
	date +"[%Y-%m-%d %H:%M:%S]"
}

function load_dependencies {
	# ensure that all dependencies are in PATH
	# also add SLIDR installation directory to PATH
	# export PATH=$PATH:/path/to/bin
	echo -e "\n$(timestamp) >>> Loading dependencies"
	module load hisat2-2.1.0
	module load samtools-1.9
	module load bedtools-2.28.0
	module load seqtk-1.3
	module load cdhit-4.8.1
	module load vsearch-2.4.3
	module load blast-2.9.0
	module load r-data.table
	module load r-3.6.1
	#module load star-2.7.1a
	module load gffread-0.11.4
	module load cutadapt-2.3
	module load bowtie2-2.3.5
	export PATH="$PATH:~/sharedscratch/apps/ViennaRNA-2.4.14/bin"
}

title="#\n# SLIDR - Spliced Leader IDentification from RNA-Seq data\n# Version 1.1\n#"

function printhelp {
	echo -e "$title"
	echo -e "Usage: bash slidr.sh -g <> | -t <> [-a <>] {-1 <> [-2 <>] | -b <> | m <>} [-r 0|1|2] [-q]"
	echo -e "\t\t[-l <>] [-x <>] [-e <>] [-D <>] [-S <>] [-A <>] [-R <>] [-O <>] [-L <>]"
	echo -e "\t\t[-o <>] [-p <>] [-c <>]"
	echo -e "Alignment reference:"
	echo -e '\t-g\tPath to genome assembly (FASTA)'
	echo -e '\t-a\tPath to genome annotations (GFF/GTF)'
	echo -e '\t-t\tPath to transcriptome assembly (FASTA)'
	echo -e "Single library:"
	echo -e '\t-1\tPath to R1 reads (FASTQ.GZ)'
	echo -e '\t-2\tPath to R2 reads (FASTQ.GZ)'
	echo -e '\t-b\tPath to read alignments (BAM)'
	echo -e '\t-r\tRead strandedness (0=unstranded; 1=stranded; 2=reverse-stranded; x=infer; default: x)'
	echo -e "\t-q\tQuality-trim 3' ends of reads and remove Illumina adapters"
	echo -e "Multiple libraries:"
	echo -e '\t-m\tPath to tab-delimited file describing library configurations:'
	echo -e '\t\tColumn 1: Library name'
	echo -e '\t\tColumn 2: Library strandedness (= -r option)'
	echo -e '\t\tColumn 3: Path to R1 reads (= -1 option)'
	echo -e '\t\tColumn 4: Path to R2 reads (= -2 option)'
	echo -e '\t\tColumn 5: "trim" to quality-trim reads (= -q option)'
	echo -e '\t\tColumn 6: Path to read alignments (= -b option)'
	echo -e "Tail filtering:"
	echo -e "\t-l\tMinimum tail length (default: 8)"
	echo -e "\t-x\tScale factor for upper tail length limit (default: 1.0)"
	echo -e "\t-e\tBLASTN E-value (default: 1)"
	echo -e "SL RNA filtering:"
	echo -e "\t-D\tSplice donor site regex (default: GT)"
	echo -e "\t-S\tSm binding site regex (default: .{40,55}AT{4,6}G)"
	echo -e "\t-A\tSplice acceptor site regex (default: AG)"
	echo -e "\t-R\tMaximum SL RNA length excluding SL (default: 80)"
	echo -e "\t-O\tMaximum SL overlap with trans-splice acceptor site (default: 10)"
	echo -e "\t-L\tMaximum base-pair span within stem-loop (default: 35)"
	echo -e "General:"
	echo -e '\t-o\tPath to output directory (default: "SLIDR_[date+time]")'
	echo -e '\t-p\tPrefix for predicted SL sequences (default: SL)'
	echo -e '\t-c\tThreads (default: all available cores)'
	echo -e '\t-h\tPrint this help page'
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
threads=$(nproc)
stranded="x"
aligner="hisat2"
clean="no"
tlength=8
tscale=1.0
evalue=1
donor="GT"
sm=".{40,55}AT{4,6}G"
acceptor="AG"
rlength=80
overlap=10
sloop=25
slprefix="SL"

# parse command-line options
cmdline="bash slidr.sh"
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
				cmdline+=" -o $2"
		shift; shift;;
		-p)		slprefix=$2
				cmdline+=" -p $2"
		shift; shift;;
	    -c)		threads=$2
				cmdline+=" -c $2"
		shift; shift;;
	    -)	            # unknown option
		shift; shift;;
	esac
done

# print input summary on screen and log file
function print_summary {
	echo -e "$title"
	echo -e "---------------- General"
	printf "%24s %s\n" "Output directory:" $outdir
	printf "%24s %s\n" "SL prefix:" $slprefix
	printf "%24s %s\n" "Threads:" $threads	
	echo -e "\n-------------- Reference"
	printf "%24s %s\n" "Genome:" $genome
	printf "%24s %s\n" "Annotations:" $ann
	printf "%24s %s\n" "Transcriptome:" $transcriptome
	echo -e "\n--------- Tail filtering"
	printf "%24s %s\n" "Min. length (bp):" $tlength
	printf "%24s %s\n" "Max. length (scale):" $tscale
	printf "%24s %s\n" "BLASTN E-value:" $evalue
	echo -e "\n------- SL RNA filtering"
	printf "%24s %s\n" "Splice donor site:" $donor
	printf "%24s %s\n" "Sm binding site:" $sm
	printf "%24s %s\n" "Splice acceptor site:" $acceptor
	printf "%24s %s\n" "Max. SL RNA length:" $rlength
	printf "%24s %s\n" "Max. acceptor overlap:" $overlap
	printf "%24s %s\n" "Max. stem-loop span:" $sloop
	echo -e "\n----------- RNA-Seq data"
	if [ ! "$design" == "" ] ; then
		printf "%24s %s\n" "Libraries:" $design
	else
		printf "%24s %s\n" "R1 reads:" $R1
		printf "%24s %s\n" "R2 reads:" $R2
		printf "%24s %s\n" "Quality trimming:" $clean
		printf "%24s %s\n" "BAM alignments:" $bam
		printf "%24s %s\n" "Strandedness:" $stranded
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
			echo "WARNING: Splice acceptor sites cannot be identified in a transcriptome reference."
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
		echo "* Library $lib (-r $stranded) -1 $R1 -2 $R2 -q $clean -b $bam"
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
echo "$(timestamp) >>> STAGE 1: Read tail extraction"

# convert GFF to GTF if required
if [ $(echo "$ann" | grep -o "...$") == "gtf" ]; then
	gtf=$ann
else
	if [ ! -f $outdir/annotations.gtf ]; then
		echo "$(timestamp) Converting GFF to GTF ..."
		gffread $ann -M -K -Q -C -T -o "$outdir/annotations.gtf"
	fi
	gtf="$outdir/annotations.gtf"	
fi

# make hisat2 index for genome
if [ ! "$genome" == "" ] && [ ! -d $outdir/hisat2_index ]; then
	mkdir -p $outdir/hisat2_index
	# if we have genome annotations, use them
	if [ ! "$ann" == "" ]; then
		# extract splice sites for hisat2
		hisat2_extract_splice_sites.py $gtf > $outdir/hisat2_index/hisat2_splicesites.txt
		hisat2_extract_exons.py $gtf > $outdir/hisat2_index/hisat2_exons.txt
		hisatgffopts="--ss $outdir/hisat2_index/hisat2_splicesites.txt
				--exon $outdir/hisat2_index/hisat2_exons.txt"
	fi
	echo "$(timestamp) Generating HISAT2 index ..."
	hisat2-build -p $threads $hisatgffopts $genome $outdir/hisat2_index/genome \
		> $outdir/hisat2_index/log_hisat2-build.txt 2>&1
fi
# make bowtie2 index for transcriptome
if [ ! "$transcriptome" == "" ] && [ ! -d $outdir/bowtie2_index ]; then
	mkdir -p $outdir/bowtie2_index
	# remove redundant transcripts
	echo "$(timestamp) Removing redundant transcripts ..."
	cd-hit-est -M 0 -c 1.0 -T $threads -i $transcriptome -o $outdir/transcriptome.cdhit \
		> $outdir/bowtie2_index/log_cd-hit-est.txt 2>&1
	transcriptome=$outdir/transcriptome.cdhit
	# make index
	echo "$(timestamp) Generating BOWTIE2 index ..."
	bowtie2-build --threads $threads \
		$transcriptome $outdir/bowtie2_index/transcriptome \
		> $outdir/bowtie2_index/log_bowtie2-build.txt 2>&1
fi

#
# 2) GENOME/TRANSCRIPTOME ALIGNMENT
#
# only when no BAM file was specified!

function trim_reads {
	trR1=$outdir/$lib/${lib/library_/}.R1.trimmed.fq.gz
	trR2=$outdir/$lib/${lib/library_/}.R2.trimmed.fq.gz
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

function star_align { #DEPRECATED
	# make STAR index
	if [ ! -d $outdir/STAR_index ]; then
		echo "$(timestamp) Generating STAR index ..."
		mkdir -p $outdir/STAR_index
		NB=12 #$(echo "l($(grep -v ">" $genome | wc -c))/l(2)/2-1" | bc -l | sed 's/\..*//g')
		if [ ! "$ann" == "" ]; then
			rlength=$(seqtk seq -l 0 $R1 | head -n 4000 | wc -L)
			stargtf="--sjdbGTFfile $ann --sjdbOverhang $(($rlength-1))" 
		fi
		STAR --runMode genomeGenerate \
			--genomeDir $outdir/STAR_index \
			--genomeFastaFiles $genome \
			--runThreadN $threads \
			--genomeSAindexNbases $([ $NB -gt 14 ] && echo "14" || echo "$NB") \
			$stargtf \
			> $outdir/log_star-genomeGenerate.txt 2>&1
	fi

	echo "$(timestamp) Aligning reads to genome ..."
	STAR --genomeDir $outdir/STAR_index \
		--readFilesIn $R1 $R2 \
		--readFilesCommand zless \
		--outSAMtype BAM Unsorted \
		--outFileNamePrefix $outdir/STAR \
		--runThreadN $threads \
		--outFilterMultimapNmax 20 \
		--alignSJoverhangMin 8 \
		--alignSJDBoverhangMin 1 \
		--outFilterMismatchNmax 999 \
		--outFilterMismatchNoverReadLmax 0.04 \
		--alignIntronMin 10 \
		--alignIntronMax 1000000 \
		--alignMatesGapMax 1000000 \
		--alignEndsType Local \
		> $outdir/log_star_align.txt
}

function genome_align {
	bam="$outdir/$lib/alignments-x$tscale.bam"
	if [ ! -f $bam ]; then
		echo "$(timestamp) Aligning reads to genome ..."
		hisat_options=""
		if [ "$R2" == "" ]; then
			hisat_options+="-U $R1"
		else
			hisat_options+="-1 $R1 -2 $R2"
		fi
		
		hisat2 --sp 1,0 --mp 3,1 --score-min L,5,$(echo "$tscale" | awk '{print -0.4*$0}') \
			-x $outdir/hisat2_index/genome $hisat_options -p $threads 2> $bam.log_hisat2.txt \
			| samtools sort -O bam -@ $threads -o $bam - 2>> $bam.log_hisat2.txt
		samtools flagstat $bam > $bam.flagstat
		samtools index -@ $threads -b $bam
		echo "$(timestamp) $(grep 'overall alignment rate' $bam.log_hisat2.txt)"
	fi
}

function transcriptome_align {
	bam="$outdir/$lib/alignments-x$tscale.bam"
	if [ ! -f $bam ]; then
		echo "$(timestamp) Aligning reads to transcriptome ..."
		bowtie_options=""
		if [ "$R2" == "" ]; then
			bowtie_options+="-U $R1"
		else
			bowtie_options+="-1 $R1 -2 $R2"
		fi
		bowtie2 --local --score-min L,5,$(echo "$tscale" | awk '{print 1-($0*0.4)}') --ma 1 --mp 3,1 -k 5 \
			-x $outdir/bowtie2_index/transcriptome $bowtie_options -p $threads 2> $bam.log_bowtie2.txt \
			| samtools sort -O bam -@ $threads -o $bam - 2>> $bam.log_bowtie2.txt
		samtools flagstat $bam > $bam.flagstat
		samtools index -@ $threads -b $bam
		echo "$(timestamp) $(grep 'overall alignment rate' $bam.log_bowtie2.txt)"
	fi
}

#
# 3) TAIL EXTRACTION
#

# function for rev compl
function revcomp {
	tr [:upper:] [:lower:] | rev | tr 'a' 'T' | tr 'c' 'G' | tr 'g' 'C' | tr 't' 'A' | tr 'n' 'N'
}
# function to convert to upper-case nucleotides
function uppercasenuc {
	tr [:lower:] [:upper:]
}

function infer_strandedness {
	if [ "$stranded" == "x" ]; then
		# infer strandedness
		if [ "$gtf" == "" ]; then
			stranded=0
		else
			# + strand
			gtfsample=$(awk '$3~"transcript|gene|mRNA" && $7=="+"' $gtf | bedtools sort -i - | bedtools merge -s -i - | bedtools sample -n 100 -i -)
			pR1=$(echo "$gtfsample" | samtools view -c -F 20 -f 64 -M -L - $bam)
			pR2=$(echo "$gtfsample" | samtools view -c -F 20 -f 128 -M -L - $bam)

			# - strand
			gtfsample=$(awk '$3~"transcript|gene|mRNA" && $7=="-"' $gtf | bedtools sort -i - | bedtools merge -s -i - | bedtools sample -n 100 -i -)
			mR1=$(echo "$gtfsample" | samtools view -c -F 4 -f 80 -M -L - $bam)
			mR2=$(echo "$gtfsample" | samtools view -c -F 4 -f 144 -M -L - $bam)

			# summarise
			R1=$((pR1+mR1))
			R2=$((pR2+mR2))
			if [ $((R1+R2)) == 0 ]; then
				# single-end data must be analysed as unstranded
				stranded=0
			else
				sR1=$((100*R1/(R1+R2)))
				sR2=$((100*R2/(R1+R2)))
				stranded=$(( 1 + 2*(sR1<40) - (sR1<60)))
			fi
			#echo -e "R1: $R1 \t $sR1 %"
			#echo -e "R2: $R2 \t $sR2 %"
			echo "$(timestamp) Inferred library strandedness: $stranded"
		fi
	fi
}

function tail_extract {
	### rev-stranded data (-s 2):	R2 are from 5' end of transcript
	###								R1 are only useful if the fragment is short;
	###								the 3' end of R1 can then read into the rev compl of the SL
	#
	# mapped to + strand: -F 20			R2: -f 128		take 5' tail from read (10S90M)
	# mapped to - strand: -F 4 -f 16	R2: -f 128		take 5' tail from read (rev compl) (90M10S)
	# mapped to + strand: -F 20			R1: -f 64		take 3' tail from read (rev compl) (90M10S)
	# mapped to - strand: -F 4 -f 16	R1: -f 64		take 3' tail from read (10S90M)

	### fwd-stranded data (-s 1):	As above, but swap R1 and R2

	### unstranded data (-s 0): 	No distinction between R1 and R2
	#
	# mapped to + strand: -F 20			take 5' tail from read (10S90M)
	# mapped to - strand: -F 4 -f 16	take 5' tail from read (rev compl) (90MS10)
	# mapped to + strand: -F 20			take 3' tail from read (rev compl) (90MS10)
	# mapped to - strand: -F 4 -f 16	take 3' tail from read (10S90M)
	
	# check if we have single-end reads:
	# if so, treat like unstranded data
	npe=$({ samtools view -H $bam ; samtools view $bam | head -n 10000; } | samtools view -c -f 1 -)
	if [ "$npe" == "0" ]; then
		stranded=0
	fi

	echo "$(timestamp) Extracting soft-clipped tails from alignments ..."

	if [ ! -f "$bam.tails+5p.fa.gz" ]; then
		samtools view -H $bam > "$bam.tails+5p.sam"
		samtools view -F 20 -f $((0+ (64*$stranded))) $bam \
			| awk '$6 ~ "^[0-9]*S"' | tee -a $bam.tails+5p.sam \
			| awk '{gsub("S.*", "", $6); t=substr($10,1,$6); print ">"$1"\n"t}' \
			| seqtk seq -U > "$bam.tails+5p.fa"
		gzip $bam.tails+5p.sam
		gzip $bam.tails+5p.fa
	fi
	if [ ! -f "$bam.tails-5p.fa.gz" ]; then
		samtools view -H $bam > "$bam.tails-5p.sam"	
		samtools view -F 4 -f $((16+ (64*$stranded))) $bam \
			| awk '$6 ~ "[0-9]*S$"' | tee -a $bam.tails-5p.sam \
			| awk '{gsub("S", "", $6); gsub(".*[A-Z]", "", $6); 
				t=substr($10,length($10)+1-$6,$6); print ">"$1"\n"t}' \
			| seqtk seq -Ur > "$bam.tails-5p.fa"
		gzip $bam.tails-5p.sam
		gzip $bam.tails-5p.fa
	fi		
	if [ ! -f "$bam.tails+3p.fa.gz" ]; then	
		samtools view -H $bam > "$bam.tails+3p.sam"
		samtools view -F 20 -f $((0+ (192-64*$stranded)*($stranded>0))) $bam \
			| awk '$6 ~ "[0-9]*S$"' | tee -a $bam.tails+3p.sam \
			| awk '{gsub("S", "", $6); gsub(".*[A-Z]", "", $6); 
				t=substr($10,length($10)+1-$6,$6); print ">"$1"\n"t}' \
			| seqtk seq -Ur > "$bam.tails+3p.fa"
		gzip $bam.tails+3p.sam
		gzip $bam.tails+3p.fa
	fi		
	if [ ! -f "$bam.tails-3p.fa.gz" ]; then
		samtools view -H $bam > "$bam.tails-3p.sam"
		samtools view -F 4 -f $((16+ (192-64*$stranded)*($stranded>0))) $bam \
			| awk '$6 ~ "^[0-9]*S"' | tee -a $bam.tails-3p.sam \
			| awk '{gsub("S.*", "", $6); t=substr($10,1,$6); print ">"$1"\n"t}' \
			| seqtk seq -U > "$bam.tails-3p.fa"
		gzip $bam.tails-3p.sam
		gzip $bam.tails-3p.fa
	fi		
	echo "$(timestamp) Extracted $(zgrep -h '^>' $bam.tails*p.fa.gz | grep -c '^>') tails ..."
}

# we need to extract tails for each library
# if a library-design file is specified, loop through libraries
# for each library, check what needs done
function library_pipeline {
	echo "$(timestamp) >>> Processing library $lib"
		lib="library_$lib"
		mkdir -p $outdir/$lib
		if [ "$bam" == "" ]; then
			trim_reads
			if [ ! "$genome" == "" ]; then
				genome_align
			elif [ ! "$transcriptome" == "" ]; then
				transcriptome_align
			fi
		fi
		infer_strandedness
		tail_extract
}
if [ ! "$design" == "" ]; then
	while IFS='#' read lib stranded R1 R2 clean bam
	do
		library_pipeline
	done < <(tr '\t' '#' < $design)
	else
		# just do it once
		lib="data"
		library_pipeline
fi

### from now on, everything will be done on all libraries together
mkdir -p $outdir/1-tails

# collect all tails
if [ ! -f $outdir/1-tails/tails-x$tscale.fa.gz ]; then
	echo "$(timestamp) Collecting tails ..."
	cat $outdir/library_*/alignments-x$tscale.bam.tails*p.fa.gz > $outdir/1-tails/tails-x$tscale.fa.gz
fi

# dereplicate 
tailprefix=$outdir/1-tails/tails-x$tscale-l$tlength
if [ ! -f $tailprefix.derep.fa.gz ]; then
	echo "$(timestamp) Dereplicating tails ..."
	vsearch --derep_fulllength $outdir/1-tails/tails-x$tscale.fa.gz \
		--output $tailprefix.derep.fa \
		--uc $tailprefix.derep.out.txt \
		--minseqlength $tlength --sizeout --threads $threads \
		> $tailprefix.derep.log-vsearch.txt 2>&1
	#awk '$1=="C"{print $9, $3}' OFS="\t" $outdir/tails/tails.derep.out.txt \
	#	> $outdir/tails/tails.derep.clusterinfo.txt
	awk '$1=="H"{print $10, $9}$1=="S"{print $9, $9}' OFS="\t" $tailprefix.derep.out.txt \
		| gzip > $tailprefix.derep.clusterinfo.txt.gz
	gzip $tailprefix.derep.fa
	gzip $tailprefix.derep.out.txt
	
	ns=$(grep -o '[0-9]* seqs' $tailprefix.derep.log-vsearch.txt | grep -o '[0-9]*')
	nt=$(grep -o '[0-9]* unique' $tailprefix.derep.log-vsearch.txt | grep -o '[0-9]*')
	echo "$(timestamp) Dereplicated $ns tails (>= $tlength bp) to $nt unique tails"	
fi

# cluster
clusterprefix=$tailprefix.cluster 
if [ ! -f $clusterprefix.centroids.fa.gz ]; then
	echo "$(timestamp) Clustering tails ..."
	vsearch --cluster_fast $tailprefix.derep.fa.gz --id 1.0 --rightjust --minwordmatches 100 \
		--uc $clusterprefix.out --centroids $clusterprefix.centroids.fa \
		--consout $clusterprefix.consensus.fa \
		--sizein --clusterout_sort -xsize \
		--minseqlength $tlength --threads $threads --iddef 0 --qmask none \
		> $clusterprefix.log-vsearch.txt 2>&1
	awk '$1=="H"{print $10, $9}$1=="S"{print $9, $9}' OFS="\t" $clusterprefix.out \
		| sed 's/;size=[0-9]\+;//g' | gzip > $clusterprefix.clusterinfo.txt.gz
	gzip $clusterprefix.out
	gzip $clusterprefix.centroids.fa
	gzip $clusterprefix.consensus.fa
	
	#cd-hit-est -i $outdir/tails/tails.unique.fa -o "$outdir/tails/02_LQ_clusters.fa" \
	#	-c 1.0 -G 0 -aS 0.1 -aL 0.1 \
	#	-d 0 -T $threads -M 0 \
	#	-l $(($tlength-1)) -g 1 -p 1 -r 0 -n 8 \
	#	> "$outdir/tails/log_cdhitest_cluster.txt" 2>&1
	#awk '/^>/{c=$2}!/^>/{gsub("^>|...$", "", $3); print c, $3}' OFS="\t" "$outdir/tails/02_LQ_clusters.fa.clstr" \
	#	> "$outdir/tails/02_LQ_clusters.fa.clusterinfo"
	echo "$(timestamp) Identified $(zgrep -c '^>' $clusterprefix.centroids.fa.gz) tail clusters"
fi

#
# 4) TAIL FILTERING
#
mkdir -p $outdir/2-RNA_filters
echo "$(timestamp) >>> STAGE 2: SL RNA identification"
# map to the genome or transcriptome
# ensure that the 3' end of the tail is aligned (5' end can contain noise, so not important)

if [ ! "$genome" == "" ]; then
	ref=$genome
elif [ ! "$transcriptome" == "" ]; then
	ref=$transcriptome
fi
echo "$(timestamp) Aligning tail cluster centroids to reference ..."
if [ ! -d $outdir/blast_db ]; then
	mkdir -p $outdir/blast_db
	makeblastdb -dbtype nucl -in $ref -out $outdir/blast_db/reference > $outdir/blast_db/log_makeblastdb.txt 2>&1
fi

blastprefix=$outdir/2-RNA_filters/centroids-x$tscale-l$tlength.blastn-e$evalue
if [ ! -f $blastprefix.out.gz ]; then
	zless $clusterprefix.centroids.fa.gz \
		| blastn -db $outdir/blast_db/reference -dust yes \
			-query - -word_size 8 -num_threads $threads -perc_identity 100 -evalue $evalue \
			-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq qlen' \
			2> $blastprefix.log_blastn.txt \
			| awk '$8==$14' > $blastprefix.out
	gzip $blastprefix.out
	#rm $outdir/1-tails/alignments.blastn
	nbl=$(zcat $blastprefix.out.gz | wc -l)
	echo "$(timestamp) Generated $nbl alignment records"
fi

#
# 5) SPLICE DONOR AND SM BINDING SITES
#
# find splice donor and Sm binding site
# if we only have a transcriptome, this could still work if we have full SL RNAs
# this is robust to empty donor or Sm sites; it will simply return empty overlaps and the full $rlength bp.
echo "$(timestamp) Identifying splice donor and Sm binding sites (${donor}${sm}) ..."

smout="$blastprefix.slRNA-R$rlength-D$donor-S$sm.txt.gz"
if [ ! -f "$smout" ]; then
	cp $ref $outdir/reference.fa
	samtools faidx $outdir/reference.fa
	zcat $blastprefix.out.gz \
	 | awk '$9<$10{print $2, $10, $10+1, $13"::"$1, ".", "+"} 
			$9>$10{print $2, ($10-2>0?$10-2:0), ($10-1>1?$10-1:1), $13"::"$1, ".", "-"}' OFS="\t" \
	 | bedtools slop -g $outdir/reference.fa.fai -i - -s -l 0 -r $(($overlap+$rlength)) \
	 | bedtools getfasta -name+ -tab -s -fi $outdir/reference.fa -bed - \
	 | awk -F'\t|::' -v p="${donor}${sm}" -v l="$rlength" -v x="$overlap" \
		'{  s=match(toupper($4), toupper(p)); o=toupper(substr($4, 1, s-1)); if(o==""){o=""}; gsub(":", "\t", $3); 
			if(s>0 && s<=(x+1)){print $2, $3, $1, o, $1""o""substr(toupper($4),RSTART,RLENGTH+l)}}' OFS="\t" \
		| gzip > "$smout"
	nsm=$(zcat "$smout" | wc -l)
	nbl=$(zcat $blastprefix.out.gz | wc -l)
	echo "$(timestamp) Identified $nsm donor splice sites and Sm binding sites ($((100*$nsm/$nbl))% of alignments)"
fi

#
# 6) SPLICE ACCEPTOR SITES
#
# we can find acceptor sites when we have a genome.
# this won't work with a transcriptome.
# we still need to collect read alignment locations.
# filter by acceptor site motif.
tailbams=$outdir/2-RNA_filters/alignments-x$tscale
spliceout=$tailbams.spliceacceptor-$acceptor-O$overlap.txt.gz
if [ ! -f "$spliceout" ]; then
	echo "$(timestamp) Extracting splice acceptor sites ($acceptor) ..."
	
	# collect read alignments
	if [ ! -f "$tailbams.tails+5p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails+5p.bam $outdir/library_*/*x$tscale.bam.tails+5p.sam.gz
	fi
	if [ ! -f "$tailbams.tails-3p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails-3p.bam $outdir/library_*/*x$tscale.bam.tails-3p.sam.gz
	fi
	if [ ! -f "$tailbams.tails-5p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails-5p.bam $outdir/library_*/*x$tscale.bam.tails-5p.sam.gz
	fi
	if [ ! -f "$tailbams.tails+3p.bam" ]; then
		samtools merge -f -@ $threads -O bam $tailbams.tails+3p.bam $outdir/library_*/*x$tscale.bam.tails+3p.sam.gz
	fi
	
	# length of acceptor site (AG: 2bp)
	# if regex is used, try and work out how many characters using bracket pairs
	# e.g., G[TA][AC] would be three characters
	acclen=$(($(echo "$acceptor" | sed 's/\[[^][]*\]/./g' | wc -m)-1))
	
	# collect acceptor regions and filter by acceptor site
	# only use reads that were clustered (shorter reads than $tlength were dropped)
	#readids=$(zless $outdir/1-tails/tails.derep.clusterinfo.txt.gz | cut -f 2)
	clinfo=$tailprefix.derep.clusterinfo.txt.gz

	{ bedtools bamtobed -i $tailbams.tails+5p.bam \
		| grep -Fw -f <(zless $clinfo | cut -f 2) - \
		| awk '{$3=$2+1; print}' OFS="\t";
	  bedtools bamtobed -i $tailbams.tails-3p.bam \
		| grep -Fw -f <(zless $clinfo | cut -f 2) - \
		| awk '{$3=$2+1; $6="+"; print}' OFS="\t";
	  bedtools bamtobed -i $tailbams.tails-5p.bam \
		| grep -Fw -f <(zless $clinfo | cut -f 2) - \
		| awk '{$2=$3-1; print}' OFS="\t";
	  bedtools bamtobed -i $tailbams.tails+3p.bam \
		| grep -Fw -f <(zless $clinfo | cut -f 2) - \
		| awk '{$2=$3-1; $6="-"; print}' OFS="\t"; } \
	| bedtools slop -g $outdir/reference.fa.fai -i - -s -l $acclen -r $overlap \
	| bedtools getfasta -fi $outdir/reference.fa -bed - -s -tab -name+ \
	| awk -F'\t|::' -v a="$acceptor" 'toupper($3)~toupper(a){gsub(":", "\t", $2); 
				gsub("/[12]$", "", $1); print $1, $2, toupper($3)}' OFS="\t" \
	| gzip > $spliceout
	
	nas=$(zcat $spliceout | wc -l)
	#nasc=$(zcat $spliceout | awk -v a="$acceptor" 'toupper($4)~toupper(a)' | wc -l)
	
	#gzip $spliceout
	#gzip $outdir/2-RNA_filters/splice_acceptor_sites_candidates.txt
	
	#echo "$(timestamp) Extracted $nasc potential splice acceptor sites ($((100*$nasc/$nas))% of $nas locations)"
	echo "$(timestamp) Extracted $nasc splice acceptor sites"
fi

#
# 7) RNAfold on unique donor sequences
#
rnafout="$smout.RNAfold-L$sloop.txt.gz"
if [ ! -f "$rnafout" ]; then
	echo "$(timestamp) Predicting RNA secondary structures ..."
	# we want the Sm site to be unpaired, so need to add a constraint string (...xxxx) below the sequence
	donlen=$(($(echo "$donor" | sed 's/\[[^][]*\]/./g' | wc -m)-1))
	smlen=$(($(echo "$sm" | sed -e 's/\[[^][]*\]/./g' -e 's/{[^}{]*}/./g' | wc -m)-2))
	zcat "$smout" | cut -f 6 | sort | uniq \
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
	#		while(c<(RSTART+d)) { printf "x"; c++ }
fi
ud=$(zcat "$rnafout" | wc -l)
echo "$(timestamp) Processed $ud unique RNAs"

#
# 8) reconstruct 3' end of SL tails and cluster
#
echo "$(timestamp) >>> STAGE 3: Consensus SL construction"

# collect all filter results
cand="$outdir/2-RNA_filters/SL.filters-x$tscale-l$tlength-e$evalue-R$rlength-D$donor-S$sm-L$sloop-A$acceptor-O$overlap.txt.gz"
if [ ! -f $cand ]; then
	echo "$(timestamp) Collecting SL RNA filter results ..."
	# merge all filter results into a single table
	LC_ALL=C
	export LC_ALL
	#echo -e "Donor_Region\tRead\tCentroid\tRepresentative\tChrom_Tail\tPos_Tail\tBLASTN_Match\tOutron_Overlap\tChrom_Gene\tPos_Gene\tAcceptor_Region\tLoops\tMFE_Frequency\tEnsemble_Diversity" \
	#	> $outdir/2-RNA_filters/SL_merged_filters.txt
	join -e '-' -t "$(printf '\t')" -1 1 -2 2 \
		<(zless $tailprefix.derep.clusterinfo.txt.gz | sort -S50% -k1,1 | uniq) \
		<(zless $clusterprefix.clusterinfo.txt.gz | sort -S50% -k 2,2 | uniq) \
		| sort -S50% -k3,3 | uniq \
	| join -e '-' -t "$(printf '\t')" -1 3 -2 1 - <(zcat "$smout" | sort -S50% -k1,1 | uniq) \
		| sort -S50% -k3,3 | uniq \
	| join -e '-' -t "$(printf '\t')" -1 3 -2 1 - <(zcat "$spliceout" | sort -S50% -k1,1 | uniq) \
		| sort -S50% -k 8,8 | uniq \
	| join -e '-' -t "$(printf '\t')" -1 8 -2 1 - <(zcat "$rnafout" | sort -S50% -k1,1 | uniq) \
	| gzip > $cand
fi
# construct SL tails and cluster them
if [ ! -f $cand.centroids.fa.gz ]; then
	echo "$(timestamp) Re-clustering SL candidates ..."
	zcat $cand | awk -F'\t' '{gsub("-", "", $8); print ">"$2"::"$7"::"$8"\n"$7$8}' \
		| gzip > $cand.SLtails.fa.gz
	#zcat $outdir/2-RNA_filters/splice_donor_sites.txt | awk -F'\t' '{print ">"$1"::"$4"::"$5"\n"$4$5}' \
	#	| gzip > $outdir/2-RNA_filters/SL_tails.fa.gz
	
	vsearch --cluster_fast $cand.SLtails.fa.gz \
		--uc $cand.clusters.out \
		--centroids $cand.centroids.fa \
		--id 1.0 --iddef 0 --minseqlength $tlength --qmask none \
		--sizeorder --maxaccepts 100 --wordlength 15 \
		--threads $threads --rightjust --minwordmatches 1000 \
		> "$cand.log_vsearch-cluster.txt" 2>&1
		
	awk -F'::|\t' '$1~"S|H"{print $2, $9, $10, $11}' OFS="\t" $cand.clusters.out \
		| gzip > $cand.clusterinfo.txt.gz
	gzip $cand.clusters.out
	gzip $cand.centroids.fa
	#gzip $outdir/2-RNA_filters/SL_merged_filters.txt
		
	#grep "^>" $outdir/2-RNA_filters/SL_clusters.msa | \
	#	awk -F'::' 'BEGIN{c=0}/consensus/{c++}!/consensus/{gsub(">\\*?", "", $1); print c, $1,$2,$3}' OFS='\t' \
	#	> $outdir/2-RNA_filters/SL_clusters.txt
fi
nc=$(grep -o "Clusters: [0-9]*" $cand.log_vsearch-cluster.txt | grep -o "[0-9]*")
echo "$(timestamp) Generated $nc clusters"

#
# 9) final filtering and consensus construction in R
#
echo "$(timestamp) Constructing final consensus SLs ..."
resultsdir=$outdir/3-results-x$tscale-l$tlength-e$evalue-R$rlength-D$donor-S$sm-L$sloop-A$acceptor-O$overlap
mkdir -p $resultsdir
Rscript slidr_analyse.R $cand "$resultsdir" " $acceptor" $slprefix $threads

echo "$(timestamp) Finished!"
exit 0

