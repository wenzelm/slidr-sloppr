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
# 4) extract SL candidate read
# 5) identify SL tails
# 6) re-align SL reads
# 7) quantify background gene expression
# 8) write checkpoint file

# get variables from arguments: 
threads=${1/-var=/}
outdir=${2/-var=/}
lib=${3/-var=/}
stranded=${4/-var=/}
R1=${5/-var=/}
R2=${6/-var=/}
clean=${7/-var=/}
bam=${8/-var=/}
slength=${9/-var=/}
err=${10/-var=/}
featureid=${11/-var=/}
metafeatureid=${12/-var=/}

# all file paths are already verified by parent script
#echo "Have received: threads=$threads outdir=$outdir lib=$lib stranded=$stranded R1=$R1 R2=$R2 clean=$clean bam=$bam"

# annotations already exist
ann="$outdir/annotations.gtf"

# SL FASTA
sls="$outdir/SLs.fa"

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
			cutadapt -a AGATCGGAAGAGC -q 20 -m 20 -O 1 \
				-o $R1 -j $threads $raw \
				> $outdir/$lib/log_cutadapt.txt 2>&1
		else
			# paired-end data
			raw1=$R1
			raw2=$R2
			R1=$trR1
			R2=$trR2
			cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -m 20 -O 1 \
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
		samtools flagstat -@ $threads $outbam > $outbam.flagstat
		samtools index -@ $threads -b $outbam
		echo "$(timestamp) $(grep 'overall alignment rate' $outdir/$lib/log_hisat2.txt)"
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
		bowtie2 --end-to-end --no-discordant \
			-x $outdir/bowtie2_index/transcriptome $bowtie_options -p $threads 2> $outbam.log_bowtie2.txt \
			| samtools sort -m 3G -@ $threads -O bam -o $outbam - 2>> $outbam.log_bowtie2.txt
		samtools flagstat -@ $threads $outbam > $outbam.flagstat
		samtools index -@ $threads -b $outbam
		echo "$(timestamp) $(grep 'overall alignment rate' $outbam.log_bowtie2.txt)"
	fi
}

function infer_strandedness {
	# check whether we have any PE reads; we'll need this info a few times
	npe=$({ samtools view -H $bam ; samtools view $bam | head -n 10000; } | samtools view -c -f 1 -)
	
	if [[ ! "$stranded" =~ [012] ]]; then
		if [ -f $outbam.strandedness.txt ]; then
			stranded=$(cat $outbam.strandedness.txt)
		else
			# infer strandedness if annotations are available
			# otherwise, strandedness must be 0
			ann="$outdir/annotations.gtf"
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
				fi
			fi
			echo "$stranded" > $outbam.strandedness.txt			
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
			samtools view -@ $threads -h -f 4 $bam \
				| samtools sort -@ $threads -n - \
				| samtools fastq -@ $threads -n - 2> $outdir/$lib/log_samtoolsfastq.txt | gzip > $outdir/$lib/SL.candidates.fastq.gz
		else	# PE reads
			# unmapped reads
			samtools view -@ $threads -h -f $((4+ (64*$stranded))) -F 264 $bam \
				| samtools sort -@ $threads -n - \
				| samtools fastq -@ $threads -n -s $outdir/$lib/SL.candidates.fastq.gz - \
				> $outdir/$lib/log_samtoolsfastq.txt 2>&1
			# mapped mates
			samtools view -@ $threads -h -f $((8+ (192-64*$stranded)*($stranded>0))) -F 260 $bam \
				| samtools sort -@ $threads -n - \
				| samtools fastq -@ $threads -n -s $outdir/$lib/SL.mates.fastq.gz - \
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
				--untrimmed-output $outdir/$lib/rescued-n$slength-e$err.5p.fastq.gz \
				-o $outdir/$lib/{name}-n$slength-e$err.5p.fastq.gz - > $outdir/$lib/log_SL.cutadapt-n$slength-e$err.txt 2>&1
		else
			# PE
			cutadapt -g file:$sls \
				-O $slength -e $err -m 15 \
				--untrimmed-output $outdir/$lib/rescued-n$slength-e$err.5p.fastq.gz \
				--untrimmed-paired-output $outdir/$lib/rescued-n$slength-e$err.mates.fastq.gz \
				-o $outdir/$lib/{name}-n$slength-e$err.5p.fastq.gz \
				-p $outdir/$lib/{name}-n$slength-e$err.mates.fastq.gz \
				$outdir/$lib/SL.candidates.fastq.gz \
				$outdir/$lib/SL.mates.fastq.gz > $outdir/$lib/log_SL.cutadapt-n$slength-e$err.txt 2>&1
		fi
		nadapt=$(grep "with adapter" $outdir/$lib/log_SL.cutadapt-n$slength-e$err.txt | head -n 1 | sed 's/.*  //g')
		echo "$(timestamp) $nadapt SL tails found among candidate reads"
	fi
}

function sl_realign {
	# align SL reads (and rescued reads) back to the genome
	for SL in $(grep "^>" $sls | tr -d ">") rescued
	do
		if [ ! -f $outdir/$lib/$SL-n$slength-e$err.*bam ]; then
			if [ -d "$outdir/hisat2_index" ]; then
				aligner_e2e="hisat2 --no-softclip -x $outdir/hisat2_index/genome"
				aligner_soft="hisat2 -x $outdir/hisat2_index/genome"
			elif [ -d "$outdir/bowtie2_index" ]; then
				aligner_e2e="bowtie2 --end-to-end -x $outdir/bowtie2_index/transcriptome"
				aligner_soft="bowtie2 --local -x $outdir/bowtie2_index/transcriptome"
			fi
			
			# SE or PE data
			hisat_options=""
			if [ "$npe" == "0" ]; then
				hisat_options+="-U $outdir/$lib/$SL-n$slength-e$err.5p.fastq.gz"
			else
				hisat_options+="-1 $outdir/$lib/$SL-n$slength-e$err.5p.fastq.gz -2 $outdir/$lib/$SL-n$slength-e$err.mates.fastq.gz"
			fi
			# rescued or SL
			if [ "$SL" == "rescued" ]; then
				# rescued reads are allowed softclipping
				# keep mates if we have PE data
				suff=""
				echo "$(timestamp) Aligning $SL reads ..."
				eval "$aligner_soft $hisat_options -p $threads 2> $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt \
					| samtools sort -O bam -@ $threads -o $outdir/$lib/$SL-n$slength-e$err$suff.bam - 2>> $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt"
			else
				# trimmed SL reads must align end-to-end
				# align mates (if available) to aid alignment,
				# but remove them from the BAM file afterwards (we want to quantify the 5' end only)
				suff=".fc"
				nr=" $(grep -Fw -A 2 "$SL" $outdir/$lib/log_SL.cutadapt-n$slength-e$err.txt | tail -n 1 | grep -o "[0-9]* times" | grep -o "[0-9]*")"
				echo "$(timestamp) Aligning$nr $SL reads ..."
				eval "$aligner_e2e $hisat_options -p $threads 2> $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt \
					| samtools view --remove-flags 1 -@ $threads -h -f $((64*($npe>0))) - 2>> $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt \
					| samtools sort -O bam -@ $threads -o $outdir/$lib/$SL-n$slength-e$err$suff.bam - 2>> $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt"
					#| samtools sort -@ $threads -n - 2>> $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt \
					#| samtools fixmate -@ $threads - - 2>> $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt \
					#
			fi
			samtools flagstat $outdir/$lib/$SL-n$slength-e$err$suff.bam > $outdir/$lib/$SL-n$slength-e$err$suff.bam.flagstat 
			samtools index -b $outdir/$lib/$SL-n$slength-e$err$suff.bam
			echo "$(timestamp) $(grep 'overall alignment rate' $outdir/$lib/log_$SL-n$slength-e$err.hisat2.txt)"
		fi
	done
}

function quantify_background {
	# PE data needs -p flag in featureCounts
	if [ "$npe" == "0" ]; then
		fcpaired=""
	else
		fcpaired="-p --countReadPairs -C"
	fi

	# quantify end2end alignments
	# use actual strandedness
	if [ ! -f "$outbam.fc-f$featureid-F$metafeatureid.txt" ]; then
		echo "$(timestamp) Quantifying background gene coverage ..."	
		if [ ! -f $outbam ]; then
			ln -fs $(pwd)/$bam $outbam
		fi
		featureCounts -a $ann -o $outbam.fc-f$featureid-F$metafeatureid.txt \
			-t $featureid -g $metafeatureid -s $stranded $fcpaired -O -M -T $threads \
			$outbam > $outdir/$lib/log_featureCounts.bg-f$featureid-F$metafeatureid.txt 2>&1
	fi
	
	# quantify rescued candidate reads (not end2end and no SL)
	# due to the pseudo-stranding done during candidate read extraction, strandedness is 1 (or 0 for unstranded data)
	if [ ! -f $outdir/$lib/rescued-n$slength-e$err.bam.fc-f$featureid-F$metafeatureid.txt ]; then
		echo "$(timestamp) Quantifying rescued reads ..."	
		featureCounts -a $ann -o $outdir/$lib/rescued-n$slength-e$err.bam.fc-f$featureid-F$metafeatureid.txt \
			-t $featureid -g $metafeatureid -s $(($stranded>0)) $fcpaired -O -M -T $threads \
			$outdir/$lib/rescued-n$slength-e$err.bam > $outdir/$lib/log_featureCounts.rescued-n$slength-e$err-f$featureid-F$metafeatureid.txt 2>&1
	fi
}

##########################
## main pipeline        ##
##########################

mkdir -p $outdir/$lib

# expected alignment output file when processing reads
outbam="$outdir/$lib/end2end_pre-align.bam"

# when BAM file is provided via $bam, no alignments to be done
if [ "$bam" == "" ]; then
	#	
	# 1) trim reads
	#
	trim_reads
	
	#
	# 2) align reads to reference (note: implement transcriptome mode)
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
# 4) extract candidate reads
#
read_extract

#
# 5) screen for SL tails
#
sl_screen

#
# 6) re-align SL reads
#
sl_realign

#
# 7) quantify background gene expression
#
quantify_background	
	
#
# 8) write checkpoint file
#
touch "$outdir/$lib/library-n$slength-e$err-f$featureid-F$metafeatureid.OK"

exit 0