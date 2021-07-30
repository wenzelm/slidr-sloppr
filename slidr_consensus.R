#!/usr/bin/env Rscript

slidr_version <- "1.1.5"

library("parallel")
library("data.table")

args = commandArgs(trailingOnly=TRUE)
wdir <- args[1]
acc <- gsub(" ", "", args[2])
sl.prefix <- args[3]
ncores <- as.numeric(args[4])

# data table threads
setDTthreads(ncores)

#
# 1) read all filter files into data tables
#

cat("... ")

# read filters.fofn
fofn <- readLines(file.path(wdir, "filters.fofn"))

# read info
ri <- fread(cmd=paste0("zcat '", fofn[1], "'"),
			sep="\t", header=F, showProgress=F, 
			col.names=c("Read","Centroid"))

# acceptor info
ai <- fread(cmd=paste0("zcat '", fofn[2], "'"),
			sep="\t", header=F, showProgress=F, 
			col.names=c("Read","Chrom_Gene","Pos_Gene", "Acceptor"))
ai <- merge(ri, ai, by="Read", allow.cartesian=TRUE, sort=F)
ai[,Read:=gsub("aln[^_]*_", "", Read)]
rm(ri)
setkeyv(ai, "Centroid")

# donor info
di <- fread(cmd=paste0("zcat '", fofn[3], "'"),
			sep="\t", header=F, showProgress=F, 
			col.names=c("Centroid","Size","Chrom_Tail","Pos_Tail",
					"BLASTN_Match","Outron_Overlap","Donor_Region"))
set(di,,"Size", NULL)					
di[Outron_Overlap=="-",Outron_Overlap:=""]
di[,Tail:=paste0(BLASTN_Match, Outron_Overlap)]

# add SL tail clusterinfo
si <- fread(cmd=paste0("zcat '", fofn[5], "'"),
			sep="\t", header=F, showProgress=F, 
			col.names=c("Cluster", "Centroid","BLASTN_Match", "Outron_Overlap"))
di <- merge(di, si, by=c("Centroid","BLASTN_Match", "Outron_Overlap"), allow.cartesian=TRUE, sort=F)
rm(si)
set(di,,"BLASTN_Match", NULL)

ai <- merge(di, ai, by="Centroid", allow.cartesian=TRUE, sort=F)
cat(paste0(uniqueN(ai[,Read]), " reads (", ai[,.N], " records) loaded\n... "))

#
# 2) resolve splice acceptor sites
#
acc.len <- nchar(gsub("\\[[^][]*\\]", ".", acc))
# upstream acceptor sequence must be identical to Outron_Overlap from the tail
ai[, acc.upstr:=substr(Acceptor, acc.len+1, acc.len+nchar(Outron_Overlap))]
ai[, Acceptor:=substr(Acceptor, nchar(Outron_Overlap)+1, nchar(Outron_Overlap)+acc.len)]
ai[acc.upstr!=Outron_Overlap, Acceptor:="#NA#"]
ai <- ai[grep(paste0("^", acc, "$"), Acceptor)]
set(ai,,"acc.upstr", NULL)
set(ai,,"Acceptor", NULL)

cat(paste0(uniqueN(ai[,Read]), " reads (", ai[,.N]," records) with splice acceptor site '", acc, "'\n... "))

#
# 3) correct genomic positions by Outron_Overlap
# resolve 3' position of SL 
# Pos_Tail is the 3' position of the tail (5' position can be truncated depending on tail)
ai[, Strand_Tail:=substr(Pos_Tail, nchar(Pos_Tail)-1, nchar(Pos_Tail)-1)]
ai[, offs:=as.integer(nchar(Outron_Overlap))]
ai[, Start_Tail:=as.integer(substr(Pos_Tail, 1, regexpr("-", Pos_Tail)-1))]
ai[, Stop_Tail:=as.integer(substr(Pos_Tail, regexpr("-", Pos_Tail)+1, regexpr("\\(", Pos_Tail)-1))]
# +: Start=Tailstart, Stop=Stop+offs
ai[Strand_Tail=="+", Stop_Tail:=Start_Tail+offs]
ai[Strand_Tail=="+", Start_Tail:=Stop_Tail-nchar(Tail)+1]
ai[Strand_Tail=="+", Donor_Pos:=Stop_Tail+1]
# -: Start=Start-offs, Stop=Tailend
ai[Strand_Tail=="-", Start_Tail:=Stop_Tail-offs+1]
ai[Strand_Tail=="-", Stop_Tail:=Start_Tail+nchar(Tail)-1]
ai[Strand_Tail=="-", Donor_Pos:=Start_Tail-1]

set(ai,,"Pos_Tail", NULL)

# need to shift 5' gene position by overlap
# new position is location of splice acceptor site
ai[, Strand_Gene:=substr(Pos_Gene, nchar(Pos_Gene)-1, nchar(Pos_Gene)-1)]
ai[, Start_Gene:=as.integer(substr(Pos_Gene, 1, regexpr("-", Pos_Gene)-1))]
ai[, Stop_Gene:=as.integer(substr(Pos_Gene, regexpr("-", Pos_Gene)+1, regexpr("\\(", Pos_Gene)-1))]

# +: Start=Start+offs, Stop=Start+acc.len
ai[Strand_Gene=="+", Start_Gene:=Start_Gene+offs+1]
ai[Strand_Gene=="+", Stop_Gene:=Start_Gene+acc.len-1]
# -: Start=Stop-acc.len, Stop=Stop-offs
ai[Strand_Gene=="-", Stop_Gene:=Stop_Gene-offs]
ai[Strand_Gene=="-", Start_Gene:=Stop_Gene-acc.len+1]

set(ai,,"offs", NULL)
set(ai,,"Pos_Gene", NULL)
set(ai,,"Outron_Overlap", NULL)

#
# 4) keep reads where the SL gene is at least 1000bp away from trans-spliced gene
#
# distance between tail and gene
ai[Strand_Tail=="+", Distance:=abs(Start_Gene-Stop_Tail)]
ai[Strand_Tail=="-", Distance:=abs(Start_Tail-Stop_Gene)]

# different chromosomes or different strands => Infite distance
ai[((Chrom_Tail != Chrom_Gene) | (Strand_Tail != Strand_Gene)), Distance:=-1]

ai <- ai[Distance<0 | Distance>1000]

cat(paste0(uniqueN(ai[,Read]), " reads (", ai[,.N]," records) with plausible SLTS pattern (>1 kbp between donor and acceptor site)\n... "))

#
# 5) resolve cluster consensus
#
clustcons <- function(tl){
	tail.len <- nchar(tl)
	cutoff <- median(max(tail.len)-tail.len)
	tail.cons <- gsub(paste0("(^.{", cutoff, "})"), "\\L\\1", tl[which.max(tail.len)[1]], perl=T)
	maj.cons <- (tail.len>=max(tail.len)-cutoff)
	
	#tail.align <- paste(sapply(max(nchar(tl))-nchar(tl), 
		#function(x) paste(rep("-", times=x), collapse="")),	tl, sep="")
	#tail.cons <- paste(apply(do.call(rbind, strsplit(tail.align, "")), 2, 
		#function(x) {
			#cons <- sort(table(x[x!="-"]), decreasing=T)[1]
			#ifelse(cons/length(x)>=0.5, names(cons), tolower(names(cons)))
		#}), collapse="")
	#maj.len <- length(gregexpr("[A-Z]", tail.cons[1])[[1]])
	#maj.cons <- (nchar(tl)>=maj.len)
	list(tail.cons, maj.cons)
}
ai[, c('Consensus','Maj_Consensus'):=clustcons(Tail), by=Cluster]

cat(paste(uniqueN(ai[,Consensus]), "consensus SLs\n... ")) 

#
# 6) summarise
#
# for each full tail, we want to know:
# - how many reads support the tail
# - how many genes the tail is spliced to

# RNAfold info
rna <- fread(cmd=paste0("zcat '", fofn[4], "'"),
			sep="\t", header=F, showProgress=F, 
			col.names=c("Donor_Region", "Loops","MFE_Frequency","Ensemble_Diversity"))
rna[,Loops:=sapply(gregexpr("\\(\\.+\\)", Loops), function(x) length(attr(x, "match.length")))]

# first, generate summary table to prioritise SLs by coverage
summ <- merge(ai, rna, by="Donor_Region", allow.cartesian=T, sort=F)[, 
	.(Length=unique(nchar(Consensus)),	
	  Reads=uniqueN(Read),	
	  SL_RNA_Genes=uniqueN(.SD[Maj_Consensus==TRUE, ], by=c("Chrom_Tail", "Donor_Pos", "Strand_Tail")),
	  SLTS_Sites=uniqueN(.SD, by=c("Chrom_Gene", "Start_Gene", "Stop_Gene", "Strand_Gene")),				
	  Stem_Loops=paste(.SD[Maj_Consensus==TRUE, sort(unique(Loops))], collapse=";"),
	  MFE_Frequency=paste(.SD[Maj_Consensus==TRUE, sort(unique(MFE_Frequency))], collapse=";"),
	  Ensemble_Diversity=paste(.SD[Maj_Consensus==TRUE, sort(unique(Ensemble_Diversity))], collapse=";")
	), by=Consensus]

setorderv(summ, c("SLTS_Sites", "Reads"), order=-1)
write.table(summ, file.path(wdir, "raw.tsv"), row.names=F, col.names=T, quote=F, sep="\t")

summ <- summ[SLTS_Sites>1 & Reads>1]
summ[,Name:=paste0(sl.prefix, 1:.N)]
setcolorder(summ, c("Name", names(summ)[-ncol(summ)]))
write.table(summ, file.path(wdir, "SLs.tsv"), row.names=F, col.names=T, quote=F, sep="\t")

cat(paste(nrow(summ), "SLs with at least 2 reads and spliced to at least 2 genes\n... ")) 
cat("Top 10 SLs:\n")
print(head(summ[,c("Consensus", "Reads", "SL_RNA_Genes", "SLTS_Sites")], n=10), row.names=F)

#
# 7) write FASTA and GFF files
#
cat("Writing output files ...\n")
dir.create(file.path(wdir, "SL_RNA_genes"))

# SL FASTA
writeLines(summ[, c(paste0(">", Name), Consensus), by=seq_len(nrow(summ))][[2]], file.path(wdir, "SLs.fa"))

# SLTS genes GFF
# don't filter by majority consensus - every read matters
writeSLTS <- function(d){
	sl.name <- summ[Consensus==d[1,Consensus], Name]
	if(length(sl.name)>0) {
	gff <- unique(d, by=c("Chrom_Gene", "Start_Gene", "Stop_Gene", "Strand_Gene"))[order(Chrom_Gene, Start_Gene),
		paste(Chrom_Gene, "SLIDR", "trans_splice_acceptor_site", Start_Gene, Stop_Gene, ".", Strand_Gene, ".", 
		paste0("ID=", sl.name, ".trans-splice-site-", seq(1:.N)), sep="\t")]
	writeLines(c("##gff-version 3",	paste("# predicted using SLIDR", slidr_version), gff), 
		file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".trans-splice-sites.gff3")))	
	}
	return(NA)
}
invisible(ai[, writeSLTS(.SD), by=Consensus, .SDcols=c("Consensus", "Chrom_Gene", "Start_Gene", "Stop_Gene", "Strand_Gene")])

# SL RNA genes
# write only those satisfying majority consensus (reduced risk of pseudogenes)
ai <- ai[Maj_Consensus==TRUE, ]
ai[, Len:=nchar(Donor_Region)]
setorderv(ai, c("Chrom_Tail", "Donor_Pos", "Len"), order=c(1,1,-1))			
ai <- unique(ai, by=c("Consensus", "Chrom_Tail", "Donor_Pos"))

writeSLRNA <- function(d){
	sl.name <- summ[Consensus==d[1,Consensus], Name]
	if(length(sl.name)>0) {
		slRNA <- d[, 
			.(ID=paste0(">", sl.name, ".", 1:.N), Donor_Region,
			  GFF=paste(Chrom_Tail, "SLIDR", "spliced_leader_RNA", 
				Start_Tail,	Stop_Tail, ".", Strand_Tail, ".", 
				paste0("ID=", sl.name, ".", seq(1:.N)), sep="\t"))]
		writeLines(as.vector(apply(slRNA[,-3], 1, paste)), file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".RNA_genes.fa")))
		writeLines(c("##gff-version 3",	paste("# predicted using SLIDR", slidr_version), slRNA[,GFF]), 
			file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".RNA_genes.gff3")))	
		return(slRNA)
	} else {
		return(data.table(ID="", Donor_Region="", GFF=""))
	}
}

SLRNA <- ai[Maj_Consensus==TRUE, writeSLRNA(.SD), by=Consensus, 
	.SDcols=c("Consensus", "Len", "Donor_Region", "Chrom_Tail", "Start_Tail", "Stop_Tail", "Strand_Tail", "Donor_Pos")]
SLRNA[,Order:=as.integer(gsub(paste0(">", sl.prefix), "", gsub("\\.[0-9]*", "", ID)))]
setorderv(SLRNA, "Order", order=1)

# FASTA
writeLines(as.vector(apply(SLRNA[!is.na(Order),c(2,3)], 1, paste)), 
			file.path(wdir, "SL_RNA_genes.fa"))	
# GFF
writeLines(c("##gff-version 3",	paste("# predicted using SLIDR", slidr_version), 
	SLRNA[!is.na(Order), GFF]), file.path(wdir, "SL_RNA_genes.gff3"))

q()
