#!/usr/bin/env bash
set -e

function timestamp { 
	date +"[%Y-%m-%d %H:%M:%S]"
}

#
# process a single library
#
# 1) trim reads if requested
# 2) align to reference
# 3) infer read strandedness
# 4) extract read tails
# 5) write checkpoint file

# get variables from arguments: 
threads=${1/-var=/}
outdir=${2/-var=/}
lib=${3/-var=/}
stranded=${4/-var=/}
tscale=${5/-var=/}
R1=${6/-var=/}
R2=${7/-var=/}
clean=${8/-var=/}
bam=${9/-var=/}

# all file paths have already been verified by parent script
# debug:
# echo "Have received: threads=$threads outdir=$outdir lib=$lib stranded=$stranded R1=$R1 R2=$R2 clean=$clean bam=$bam"

##########################
## function definitions ##
##########################

function trim_reads {
	trR1=$outdir/$lib/${lib/1-library_/}.R1.trimmed.fq.gz
	trR2=$outdir/$lib/${lib/1-library_/}.R2.trimmed.fq.gz
	if [ "$clean" == "trim" ] && [ ! -f $trR1 ]; then
		echo "$(timestamp) Trimming adapters and poor-quality bases from reads ..."
		if [ "$R2" == "" ]; then
			# single-end data
			raw=$R1
			R1=$trR1
			cutadapt -a AGATCGGAAGAGC -a "A{100}" -g "T{100}" \
				-n 3 -q 20 -m 20 -l 150 -O 1 \
				-o $R1 -j $threads $raw \
				> $outdir/$lib/log_cutadapt.txt 2>&1
		else
			# paired-end data
			raw1=$R1
			raw2=$R2
			R1=$trR1
			R2=$trR2
			cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
				-a "A{100}" -A "A{100}" -g "T{100}" -G "T{100}" \
				-n 3 -q 20 -m 20 -l 150 -O 1 \
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
		hisat2 --sp 1,0 --mp 3,1 -k 5 --score-min L,5,$(echo "$tscale" | awk '{print -0.4*$0}') \
			-x $outdir/hisat2_index/genome $hisat_options -p $threads 2> $outbam.log_hisat2.txt \
			| samtools sort -m 1G -@ $threads -O bam -o $outbam - 2>> $outbam.log_hisat2.txt
		samtools flagstat -@ $threads $outbam > $outbam.flagstat
		samtools index -@ $threads -b $outbam
		echo "$(timestamp) $(grep 'overall alignment rate' $outbam.log_hisat2.txt)"
	fi
}

function transcriptome_align {
	if [ ! -f $outbam ]; then
		echo "$(timestamp) Aligning reads to transcriptome ..."
		bowtie_options=""
		if [ "$R2" == "" ]; then
			bowtie_options+="-U $R1"
		else
			bowtie_options+="-1 $R1 -2 $R2"
		fi
		bowtie2 --local --score-min L,5,$(echo "$tscale" | awk '{print 1-($0*0.4)}') --ma 1 --mp 3,1 -k 5 \
			-x $outdir/bowtie2_index/transcriptome $bowtie_options -p $threads 2> $outbam.log_bowtie2.txt \
			| samtools sort -m 1G -@ $threads -O bam -o $outbam - 2>> $outbam.log_bowtie2.txt
		samtools flagstat -@ $threads $outbam > $outbam.flagstat
		samtools index -@ $threads -b $outbam
		echo "$(timestamp) $(grep 'overall alignment rate' $outbam.log_bowtie2.txt)"
	fi
}


function infer_strandedness {
	# strandedness may already exist 
	if [ -f $outbam.strandedness.txt ]; then
		stranded=$(cat $outbam.strandedness.txt)
	else
		# strandedness must be inferred only if no strandedness is given
		if [[ ! "$stranded" =~ [012] ]]; then
			# infer strandedness if annotations are available
			# otherwise, strandedness must be 0
			ann="$outdir/annotations.gtf.genes"
			if [ ! -f "$ann" ]; then
				stranded=0
			else
				# check if we have single-end reads:
				# if so, must treat like unstranded data
				npe=$({ samtools view -H $bam ; samtools view $bam | head -n 10000; } | samtools view -c -f 1 -)
				if [ "$npe" == "0" ]; then
					stranded=0
				else
					# + strand
					gtfsample=$(awk '$7=="+"' $ann \
						| bedtools sort -i - | bedtools merge -s -i - | bedtools sample -n 1000 -i -)
					pR1=$(echo "$gtfsample" | samtools view -c -F 20 -f 64 -M -L - $bam)
					pR2=$(echo "$gtfsample" | samtools view -c -F 20 -f 128 -M -L - $bam)

					# - strand
					gtfsample=$(awk '$7=="-"' $ann \
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
					
				fi
			fi
			echo "$(timestamp) Inferred library strandedness: $stranded"
		fi
		echo "$stranded" > $outbam.strandedness.txt			
	fi
}

function tail_extract {
	### rev-stranded data (-r 2):	R2 are from 5' end of transcript
	###								R1 are only useful if the fragment is short;
	###								the 3' end of R1 can then read into the rev compl of the SL
	#
	# mapped to + strand: -F 20			R2: -f 128		take 5' tail from read (10S90M)
	# mapped to - strand: -F 4 -f 16	R2: -f 128		take 5' tail from read (rev compl) (90M10S)
	# mapped to + strand: -F 20			R1: -f 64		take 3' tail from read (rev compl) (90M10S)
	# mapped to - strand: -F 4 -f 16	R1: -f 64		take 3' tail from read (10S90M)

	### fwd-stranded data (-r 1):	As above, but swap R1 and R2

	### unstranded data (-r 0) and single-end reads: 	No distinction between R1 and R2
	#
	# mapped to + strand: -F 20			take 5' tail from read (10S90M)
	# mapped to - strand: -F 4 -f 16	take 5' tail from read (rev compl) (90M10S)
	# mapped to + strand: -F 20			take 3' tail from read (rev compl) (90M10S)
	# mapped to - strand: -F 4 -f 16	take 3' tail from read (10S90M)
	
	# we want to keep secondary alignments
	# could be changed to ignore secondary alignments:
	# -F 276 instead of -F 20
	# -F 260 instead of -F 4
	
	echo "$(timestamp) Extracting soft-clipped tails from alignments ..."

	{ if [ ! -f "$outbam.tails+5p.fa.gz" ]; then
		samtools view -H $bam > "$outbam.tails+5p.sam"
		samtools view -F 20 -f $((0+ (64*$stranded))) $bam \
			| awk '$6 ~ "^[0-9]*S"{$1="aln+5p"NR"_"$1; print}' OFS="\t" | tee -a $outbam.tails+5p.sam \
			| awk '{gsub("S.*", "", $6); t=substr($10,1,$6); print ">"$1"\n"t}' \
			| seqtk seq -U > "$outbam.tails+5p.fa"
		gzip $outbam.tails+5p.sam
		gzip $outbam.tails+5p.fa
	fi; } &
	{ if [ ! -f "$outbam.tails-5p.fa.gz" ]; then
		samtools view -H $bam > "$outbam.tails-5p.sam"	
		samtools view -F 4 -f $((16+ (64*$stranded))) $bam \
			| awk '$6 ~ "[0-9]*S$"{$1="aln-5p"NR"_"$1; print}' OFS="\t" | tee -a $outbam.tails-5p.sam \
			| awk '{gsub("S", "", $6); gsub(".*[A-Z]", "", $6); 
				t=substr($10,length($10)+1-$6,$6); print ">"$1"\n"t}' \
			| seqtk seq -Ur > "$outbam.tails-5p.fa"
		gzip $outbam.tails-5p.sam
		gzip $outbam.tails-5p.fa
	fi; } &		
	{ if [ ! -f "$outbam.tails+3p.fa.gz" ]; then	
		samtools view -H $bam > "$outbam.tails+3p.sam"
		samtools view -F 20 -f $((0+ (192-64*$stranded)*($stranded>0))) $bam \
			| awk '$6 ~ "[0-9]*S$"{$1="aln+3p"NR"_"$1; print}' OFS="\t" | tee -a $outbam.tails+3p.sam \
			| awk '{gsub("S", "", $6); gsub(".*[A-Z]", "", $6); 
				t=substr($10,length($10)+1-$6,$6); print ">"$1"\n"t}' \
			| seqtk seq -Ur > "$outbam.tails+3p.fa"
		gzip $outbam.tails+3p.sam
		gzip $outbam.tails+3p.fa
	fi; } &
	{ if [ ! -f "$outbam.tails-3p.fa.gz" ]; then
		samtools view -H $bam > "$outbam.tails-3p.sam"
		samtools view -F 4 -f $((16+ (192-64*$stranded)*($stranded>0))) $bam \
			| awk '$6 ~ "^[0-9]*S"{$1="aln-3p"NR"_"$1; print}' OFS="\t" | tee -a $outbam.tails-3p.sam \
			| awk '{gsub("S.*", "", $6); t=substr($10,1,$6); print ">"$1"\n"t}' \
			| seqtk seq -U > "$outbam.tails-3p.fa"
		gzip $outbam.tails-3p.sam
		gzip $outbam.tails-3p.fa
	fi; } & wait	
	echo "$(timestamp) Extracted $(zgrep -h '^>' $outbam.tails*p.fa.gz | grep -c '^>') tails"
}

##########################
## main pipeline        ##
##########################

mkdir -p $outdir/$lib

# expected alignment output file when processing reads
outbam="$outdir/$lib/alignments-x$tscale.bam"

# when BAM file is provided via $bam, no alignments to be done
if [ "$bam" == "" ]; then
	#	
	# 1) trim reads
	#
	trim_reads
	
	#
	# 2) align reads to reference
	#
	if [ -d "$outdir/hisat2_index" ]; then
		genome_align
	elif [ -d "$outdir/bowtie2_index" ]; then
		transcriptome_align
	fi
	bam=$outbam
fi
#
# 3) infer strandedness
#
infer_strandedness

#
# 4) extract tails
#
tail_extract

#
# 5) write checkpoint file
#
touch "$outdir/$lib/library-x$tscale.OK"

exit 0