#!/usr/bin/env Rscript

slidr_version <- "1.1.3"

library("parallel")
library("data.table")

args = commandArgs(trailingOnly=TRUE)
#filterin <- args[1]
wdir <- args[1]
acc <- gsub(" ", "", args[2])
sl.prefix <- args[3]
ncores <- as.numeric(args[4])

wdir <- "slidr_Cel_SLf_final/3-results-x1.0-l8-DGC-e1-R80-DGT-S.{40,55}AC?T{4,6}G-L35-AAG-O10/"
acc="AG"
ncores=1
sl.prefix="SL"

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
ai <- merge(ri, ai, by="Read", allow.cartesian=TRUE)
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
di <- merge(di, si, by=c("Centroid","BLASTN_Match", "Outron_Overlap"), allow.cartesian=TRUE)
rm(si)
set(di,,"BLASTN_Match", NULL)

ai <- merge(di, ai, by="Centroid", allow.cartesian=TRUE)
cat(paste(uniqueN(ai[,.(Read)]), "reads loaded\n... "))

#
# 2) resolve splice acceptor sites
#
acc.len <- nchar(gsub("\\[[^][]*\\]", ".", acc))
# upstream acceptor sequence must be identical to Outron_Overlap from the tail
ai[, acc.upstr:=substring(Acceptor, acc.len+1, acc.len+nchar(Outron_Overlap))]
ai[, Acceptor:=substring(Acceptor, nchar(Outron_Overlap)+1, nchar(Outron_Overlap)+acc.len)]
ai[acc.upstr!=Outron_Overlap, Acceptor:="#NA#"]
ai <- ai[grep(paste0("^", acc, "$"), Acceptor)]
set(ai,,"acc.upstr", NULL)
set(ai,,"Acceptor", NULL)

cat(paste0(uniqueN(ai[,.(Read)]), " reads with splice acceptor site (", acc, ")\n... "))

#
# 3) correct genomic positions by Outron_Overlap
# resolve 3' position of SL 
# Pos_Tail is the 3' position of the tail (5' position can be truncated depending on tail)
startpos <- function(p) {
	s <- strsplit(p, "[-(]")
	s <- mclapply(s, mc.cores=ncores, function(x) {
		if(length(x)==3) {
			return(x[1])
		} else if(length(x)==4) {
			return(x[2])
		}
	})
	unlist(s)
}
ai[, strand:=regmatches(Pos_Tail, regexpr("\\([\\+\\-]\\)", Pos_Tail))]
ai[, offs:=nchar(Outron_Overlap)]
ai[, pos:=startpos(Pos_Tail)]
ai[, Pos_Tail:=paste0(as.numeric(pos) + (1-grepl("\\-", strand)*2)*offs, strand)]

# need to shift 5' gene position by overlap
# new position is location of splice acceptor site
ai[, strand:=regmatches(Pos_Gene, regexpr("\\([\\+\\-]\\)", Pos_Gene))] 
ai[, pos:=startpos(Pos_Gene)] #, by = seq_len(nrow(ai))
ai[, neg:=grepl("\\-", strand)]
ai[, Pos_Gene:=paste0(as.numeric(pos) + (1-neg*2)*offs + 1 -neg*(acc.len+1),
					"-",
					as.numeric(pos) + (1-neg*2)*offs + (1-neg)*acc.len,
					strand)]

set(ai,,"strand", NULL)
set(ai,,"offs", NULL)
set(ai,,"pos", NULL)
set(ai,,"neg", NULL)
set(ai,,"Outron_Overlap", NULL)

#
# 4) keep reads where the SL gene is at least 1000bp away from trans-spliced gene
#
# distance between tail and gene
ai[,Distance:=abs(as.numeric(gsub("\\(.*", "", Pos_Tail))-as.numeric(gsub("\\-.*", "", Pos_Gene)))]

# different chromosomes or different strands => Infite distance
ai[((Chrom_Tail!=Chrom_Gene) | 
	  (regmatches(Pos_Tail, regexpr("\\([\\+\\-]\\)", Pos_Tail)) != 
	   regmatches(Pos_Gene, regexpr("\\([\\+\\-]\\)", Pos_Gene)))),
   Distance:=Inf]

ai <- ai[Distance>1000]

cat(paste(uniqueN(ai[,.(Read)]), "reads with plausible SLTS pattern (>1 kbp between donor and acceptor site)\n... "))

#
# 5) resolve cluster consensus
#
clustcons <- function(tl){
	tail.align <- paste(sapply(max(nchar(tl))-nchar(tl), 
		function(x) paste(rep("-", times=x), collapse="")),	tl, sep="")
	tail.cons <- paste(apply(do.call(rbind, strsplit(tail.align, "")), 2, 
		function(x) {
			cons <- sort(table(x[x!="-"]), decreasing=T)[1]
			ifelse(cons/length(x)>=0.5, names(cons), tolower(names(cons)))
		}), collapse="")
	maj.len <- length(gregexpr("[A-Z]", tail.cons[1])[[1]])
	maj.cons <- (nchar(tl)>=maj.len)
	list(tail.cons, maj.cons)
}
ai[, c('Consensus','Maj_Consensus'):=clustcons(Tail), by=Cluster]

# Check if all Consensuses are unique to a single cluster (if not, there is a clustering problem)
#cc <- table(apply(table(ai$Consensus, ai$Cluster), 1, function(x) sum(x>0)))
#if(length(cc)>1) cat("WARNING: Some consensus tails occur in more than one tail cluster!\n")
#cc <- table(apply(table(ai$Consensus, ai$Cluster), 2, function(x) sum(x>0)))
#if(length(cc)>1) cat("WARNING: Some consensus tails occur in more than one tail cluster!\n")

cat(paste(uniqueN(ai[,.(Consensus)]), "consensus SLs\n... ")) 

#
# 6) summarise
#
# for each full tail, we want to know:
# - how many reads support the tail
# - how many genes the tail is spliced to

################################ leave RNAfold for later
#setkeyv(m2, c("Centroid", "BLASTN_Match", "Outron_Overlap"))
# RNAfold info
rna <- fread(cmd=paste0("zcat '", fofn[4], "'"),
			sep="\t", header=F, showProgress=F, 
			col.names=c("Donor_Region", "Loops","MFE_Frequency","Ensemble_Diversity"))
rna[,Loops:=sapply(gregexpr("\\(\\.+\\)", Loops), function(x) length(attr(x, "match.length")))]

# add RNAfold to donor
#setkeyv(di, "Donor_Region")
#setkeyv(rna, "Donor_Region")
#di <- di[rna]
#rm(rna)
#############################################


# first, generate summary table to prioritise SLs by coverage
rnastats <- function(d, v){
	paste(unique(rna[d[Maj_Consensus==TRUE, .(Donor_Region)], on="Donor_Region"][[v]]), collapse=";")
}
summ <- ai[, .(Length=unique(nchar(Consensus)),	
				Reads=length(unique(Read)),	
				SL_RNA_Genes=uniqueN(.SD[Maj_Consensus==TRUE, .(Chrom_Tail, Pos_Tail)]),
				SLTS_Genes=uniqueN(.SD[, .(Chrom_Gene, Pos_Gene)]), 
				Stem_Loops=rnastats(.SD, "Loops"), 
				MFE_Frequency=rnastats(.SD, "MFE_Frequency"),
				Ensemble_Diversity=rnastats(.SD, "Ensemble_Diversity")), by=Consensus]

summ <- summ[order(SLTS_Genes, Reads, decreasing=TRUE)]
write.table(summ, file.path(wdir, "raw.tsv"), row.names=F, col.names=T, quote=F, sep="\t")

summ <- summ[SLTS_Genes>1 & Reads>1]
summ[,Name:=paste0(sl.prefix, 1:.N)]
setcolorder(summ, c("Name", names(summ)[-ncol(summ)]))
write.table(summ, file.path(wdir, "SLs.tsv"), row.names=F, col.names=T, quote=F, sep="\t")

cat(paste(nrow(summ), "SLs with at least 2 reads and spliced to at least 2 genes\n... ")) 
cat("Top 10 SLs:\n")
print(head(summ[,c("Consensus", "Reads", "SL_RNA_Genes", "SLTS_Genes")], n=10), row.names=F)

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
	sl.name <- d[,Name]
	gff <- unique(ai[Consensus==d[,Consensus], .(Chrom_Gene, Pos_Gene,
		Strand=substr(Pos_Gene, nchar(Pos_Gene)-1, nchar(Pos_Gene)-1),
		Pos1=as.integer(sapply(strsplit(Pos_Gene, "[-(]"), function(x) x[1])),
		Pos2=as.integer(sapply(strsplit(Pos_Gene, "[-(]"), function(x) x[2])))])[order(Chrom_Gene, Pos1),
		.(GFF=(paste(Chrom_Gene, "SLIDR", "trans_splice_acceptor_site", Pos1, Pos2, ".", Strand, ".", 
		paste0("ID=", sl.name, ".trans-splice-site-", seq(1:.N)), sep="\t")))][,GFF]
	writeLines(c("##gff-version 3",	paste("# predicted using SLIDR", slidr_version), gff), file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".trans-splice-sites.gff3")))
	return(NA)
}
invisible(summ[, writeSLTS(.SD), by=seq_len(nrow(summ))])

# SL RNA genes
# write only those satisfying majority consensus (reduced risk of pseudogenes)
# FASTA
SLRNA_FASTA <- function(d){
	sl.name <- d[,Name]
	fasta <- unique(ai[Consensus==d[,Consensus] & Maj_Consensus==TRUE, 
		.(Chrom_Tail, Pos_Tail, Donor_Region, Len=nchar(Donor_Region))][order(Chrom_Tail, Pos_Tail, -Len)], 
		by=c("Chrom_Tail", "Pos_Tail"))[, .(ID=paste0(">", sl.name, ".", 1:.N), Donor_Region)]
	fasta <- as.vector(apply(fasta, 1, paste))
	writeLines(fasta, file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".RNA_genes.fa")))
	return(fasta)
}
writeLines(summ[, .(fasta=SLRNA_FASTA(.SD)), by=seq_len(nrow(summ))][,fasta], 
			file.path(wdir, "SL_RNA_genes.fa"))
# GFF
SLRNA_GFF <- function(d){
	sl.name <- d[,Name]
	gff <- unique(ai[Consensus==d[,Consensus] & Maj_Consensus==TRUE, 
		.(Tail, Chrom_Tail, Pos_Tail, Donor_Region, Len=nchar(Donor_Region))][order(Chrom_Tail, Pos_Tail, -Len)], 
		by=c("Chrom_Tail", "Pos_Tail"))[, c('Strand', 'Stop'):=list(
			substr(Pos_Tail, nchar(Pos_Tail)-1, nchar(Pos_Tail)-1),
		    as.integer(gsub("\\(.\\)", "", Pos_Tail)))][,
		c('ID', 'Start'):=list(seq(1:.N),
		Stop-(as.integer(paste0(Strand, "1"))*nchar(Tail)))][][,
		.(GFF=(paste(Chrom_Tail, "SLIDR", "spliced_leader_RNA", 
		max(1, min(c(Start, Stop))+1), 
		max(c(Start, Stop)), 
		".", Strand, ".", 
		paste0("ID=", sl.name, ".", ID), sep="\t"))), by=ID][,GFF]
	writeLines(c("##gff-version 3",	paste("# predicted using SLIDR", slidr_version), gff), file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".RNA_genes.gff3")))	
	return(gff)
}
writeLines(c("##gff-version 3",	paste("# predicted using SLIDR", slidr_version), 
	summ[, .(gff=SLRNA_GFF(.SD)), by=seq_len(nrow(summ))][,gff]), file.path(wdir, "SL_RNA_genes.gff3"))

q()
