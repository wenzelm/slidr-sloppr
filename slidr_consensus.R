#!/usr/bin/env Rscript

library("parallel")
library("data.table")

args = commandArgs(trailingOnly=TRUE)
filterin <- args[1]
wdir <- args[2]
acc <- gsub(" ", "", args[3])
sl.prefix <- args[4]
ncores <- as.numeric(args[5])

#wdir <- "slidr_spliced_Cel_SRR1585277_STAR"
#wdir <- "slidr_SLfinder_transtest"
#wdir <- "slidr_Cint_all"
#wdir <- "slidr_Hvul_SLf"
#wdir <- "slidr_Ppunc"
#wdir <- "slidr_Cbrigs_Uyar"
#acc <- "AG"
#ncores=16
#sl.prefix <- "SL"

#clust <- makeCluster(ncores)
setDTthreads(ncores)

# 1) read merged filter file
#m <- fread(cmd=paste0("gunzip -c ", file.path(wdir, "2-RNA_filters", "SL_merged_filters.txt.gz")), 
m <- fread(cmd=paste0("zcat '", filterin, "'"), 
		sep="\t", header=F, showProgress=F, 
		col.names=c("Donor_Region","Read","Centroid","Representative","Size", "Chrom_Tail","Pos_Tail",
					"BLASTN_Match","Outron_Overlap","Chrom_Gene","Pos_Gene","Acceptor_Region",
					"Loops","MFE_Frequency","Ensemble_Diversity"))
#m$Read <- gsub("aln[^_]*_", "", m$Read)
m$Outron_Overlap[m$Outron_Overlap=="-"] <- ""
m$Tail <- paste0(m$BLASTN_Match, m$Outron_Overlap)
m$Loops <- sapply(gregexpr("\\(\\.+\\)", m$Loops), function(x) length(attr(x, "match.length")))
setkeyv(m, c("Centroid", "BLASTN_Match", "Outron_Overlap"))


# 2) add SL cluster ID
hq <- fread(cmd=paste0("zcat '", filterin, ".clusterinfo.txt.gz'"), sep="\t", header=F, showProgress=F,
					colClasses="character",
					col.names=c("Cluster", "Centroid", "BLASTN_Match", "Outron_Overlap"),
					key=c("Centroid", "BLASTN_Match", "Outron_Overlap"))
m <- merge(unique(hq), m, by=c("Centroid", "BLASTN_Match", "Outron_Overlap"), sort=F, allow.cartesian=TRUE)
rm(hq)
cat(paste(length(unique(m$Read)), "reads loaded\n"))

# 3) resolve splice acceptor site
acc.len <- nchar(gsub("\\[[^][]*\\]", ".", acc))
m$Acceptor.Splicesite <- substring(m$Acceptor_Region, nchar(m$Outron_Overlap)+1, nchar(m$Outron_Overlap)+acc.len)

# additional sanity check: is the upstream outron sequence identical to Outron_Overlap from the tail?
ovl <- substring(m$Acceptor_Region, acc.len+1, acc.len+nchar(m$Outron_Overlap))
m$Acceptor.Splicesite[ovl!=m$Outron_Overlap] <- "#NA#"
m <- m[grep(paste0("^", acc, "$"), m$Acceptor.Splicesite),]
cat(paste0(length(unique(m$Read)), " reads with splice acceptor site (", acc, ")\n"))

# 4) correct genomic positions using Outron_Overlap
# resolve 3' position of SL 
# Pos_Tail is the 3' position of the tail (5' position can be truncated depending on tail)
strand <- regmatches(m$Pos_Tail, regexpr("\\([\\+\\-]\\)", m$Pos_Tail))
offs <- nchar(m$Outron_Overlap)
pos <- sapply(strsplit(m$Pos_Tail, "[-(]"), function(x) {
	if(length(x)==3) return(x[1])
	if(length(x)==4) return(x[2])
})
m$Pos_Tail <- paste0(as.numeric(pos) + (1-grepl("\\-", strand)*2)*offs, strand)

# need to shift 5' gene position by overlap
# new position is location of splice acceptor site 
strand <- regmatches(m$Pos_Gene, regexpr("\\([\\+\\-]\\)", m$Pos_Gene))
pos <- sapply(strsplit(m$Pos_Gene, "[-(]"), function(x) {
	if(length(x)==3) return(x[1])
	if(length(x)==4) return(x[2])
})
neg <- grepl("\\-", strand)
m$Pos_Gene <- paste0(as.numeric(pos) + (1-neg*2)*offs + 1 -neg*(acc.len+1),
					"-",
					as.numeric(pos) + (1-neg*2)*offs + (1-neg)*acc.len,
					strand)

old <- function() {
m$Pos_Tail <- apply(m, 1, function(rd){
#m$Pos_Tail <- do.call(c, mclapply(mc.cores=ncores, 1:nrow(m), function(rd){
#	rd <- unlist(m[rd,])
	strand <- regmatches(rd["Pos_Tail"], regexpr("\\([\\+\\-]\\)", rd["Pos_Tail"]))
	pos <- unlist(strsplit(rd["Pos_Tail"], "[-(]"))[1:2]
	offs <- nchar(rd["Outron_Overlap"])
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, strand)
	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs, strand)
	#if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, "-", as.numeric(pos[1])+offs+2, strand)
	#if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs-2, "-", as.numeric(pos[2])-offs, strand)
	newpos
})
m$Pos_Gene <- apply(m, 1, function(rd){
#m$Pos_Gene <- do.call(c, mclapply(mc.cores=ncores, 1:nrow(m), function(rd){
#	rd <- unlist(m[rd,])
	strand <- regmatches(rd["Pos_Gene"], regexpr("\\([\\+\\-]\\)", rd["Pos_Gene"]))
	pos <- unlist(strsplit(rd["Pos_Gene"], "[-(]"))[1:2]
	offs <- nchar(rd["Outron_Overlap"])
	#if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, strand)
	#if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs, strand)
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs+1, "-", as.numeric(pos[1])+offs+acc.len, strand)
	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs-acc.len, "-", as.numeric(pos[2])-offs, strand)
	newpos
})
# 5) keep reads where the SL gene is at least 1000bp away from trans-spliced gene
m$Distance <- apply(m, 1, function(rd){
#m$Distance <- do.call(c, mclapply(mc.cores=ncores, 1:nrow(m), function(rd){
#	rd <- unlist(m[rd,])
	if(rd["Chrom_Tail"] != rd["Chrom_Gene"]) {	# different chromosomes; easy
		di <- Inf
	} else {	# same chromosomes
		# check strands
		s1 <- regmatches(rd["Pos_Tail"], regexpr("\\([\\+\\-]\\)", rd["Pos_Tail"]))
		s2 <- regmatches(rd["Pos_Gene"], regexpr("\\([\\+\\-]\\)", rd["Pos_Gene"]))
		if(s1 != s2) {	# different strands
			di <- Inf
		} else {	# same strands
			p1 <- unlist(strsplit(rd["Pos_Tail"], "[-(]"))
			p2 <- unlist(strsplit(rd["Pos_Gene"], "[-(]"))
			if(s1=="(+)") di <- abs(as.numeric(p1[1])-as.numeric(p2[1]))
			if(s1=="(-)") di <- abs(as.numeric(p1[2])-as.numeric(p2[2]))
		}
	}
	di
})
}

# 5) keep reads where the SL gene is at least 1000bp away from trans-spliced gene
# distance between tail and gene
m$Distance <- abs(as.numeric(gsub("\\(.*", "", m$Pos_Tail))-as.numeric(gsub("\\-.*", "", m$Pos_Gene)))
# different chromosomes or different strands => Infite distance
dif <- (m$Chrom_Tail!=m$Chrom_Gene) | 
	  (regmatches(m$Pos_Tail, regexpr("\\([\\+\\-]\\)", m$Pos_Tail)) != 
	   regmatches(m$Pos_Gene, regexpr("\\([\\+\\-]\\)", m$Pos_Gene)))
m$Distance[dif] <- Inf

m <- subset(m, Distance>1000)
cat(paste(length(unique(m$Read)), "reads with plausible SLTS pattern (>1 kbp between donor and acceptor site)\n"))

# 6) resolve cluster consensus
m <- mclapply(mc.cores=ncores, split(m, m$Cluster), function(cls){
	#cls$Consensus <- rep(cls$Tail[which.max(nchar(cls$Tail))], times=nrow(cls))
	tail.align <- paste(sapply(max(nchar(cls$Tail))-nchar(cls$Tail), 
		function(x) paste(rep("-", times=x), collapse="")), 
		cls$Tail, sep="")
	cls$Consensus <- paste(apply(do.call(rbind, strsplit(tail.align, "")), 2, 
		function(x) {
			cons <- sort(table(x[x!="-"]), decreasing=T)[1]
			ifelse(cons/length(x)>=0.5, names(cons), tolower(names(cons)))
		}), collapse="")
	maj.len <- length(gregexpr("[A-Z]", cls$Consensus[1])[[1]])	
	#cls$Maj_Consensus <- sapply(cls$Tail, function(x) grepl(x, cls$Consensus[1], ignore.case=F))
	cls$Maj_Consensus <- nchar(cls$Tail)>=maj.len
	return(cls)
})
m <- do.call(rbind, m)
# Check if all Consensuses are unique to a single cluster (if not, there is a clustering problem)
cc <- table(apply(table(m$Consensus, m$Cluster), 1, function(x) sum(x>0)))
if(length(cc)>1) cat("WARNING: Some consensus tails occur in more than one tail cluster!\n")
cc <- table(apply(table(m$Consensus, m$Cluster), 2, function(x) sum(x>0)))
if(length(cc)>1) cat("WARNING: Some consensus tails occur in more than one tail cluster!\n")

cat(paste(length(unique(m$Consensus)), "consensus SLs\n")) 

# 7) summarise
# for each full tail, we want to know:
# - how many reads support the tail
# - how many genes the tail is spliced to

# first, generate summary table to prioritise SLs by coverage
cons <- split(m, m$Consensus)
#cons <- split(m[m$Loops==2,], m[m$Loops==2,]$Consensus)

summ <- lapply(cons, function(tl){
	data.frame(Sequence=unique(tl$Consensus),
	Length=unique(nchar(tl$Consensus)),
	Reads=length(unique(tl$Read)),
	#Unique_Reads=length(unique(tl$Representative)),
	#Unique_Tails=length(unique(tl$Tail)),
	SL_RNA_Genes=nrow(unique(tl[tl$Maj_Consensus==TRUE,c("Chrom_Tail", "Pos_Tail")])),
	SLTS_Genes=nrow(unique(tl[,c("Chrom_Gene", "Pos_Gene")])),
	Stem_Loops=paste(unique(tl$Loops), collapse=";"),
	MFE_Frequency=paste(unique(tl$MFE_Frequency), collapse=";"),
	Ensemble_Diversity=paste(unique(tl$Ensemble_Diversity), collapse=";"),
	stringsAsFactors=F)
})
summ <- do.call(rbind, summ)
summ <- summ[rev(order(summ$SLTS_Genes, summ$Reads)),]
write.table(summ, file.path(wdir, "raw.tsv"), row.names=F, col.names=T, quote=F, sep="\t")

summ <- subset(summ, SLTS_Genes>1 & Reads>1)
summ <- data.frame(Name=paste0(sl.prefix, 1:nrow(summ)), summ)
write.table(summ, file.path(wdir, "SLs.tsv"), row.names=F, col.names=T, quote=F, sep="\t")

cat(paste(nrow(summ), "SLs with at least 2 reads and spliced to at least 2 genes\n")) 
cat("Top 10 SLs:\n")
print(head(summ[,c("Sequence", "Reads", "SL_RNA_Genes", "SLTS_Genes")], n=10), row.names=F)

# then, write donatron FASTA and GFF files
cat("Writing output files ...\n")
dir.create(file.path(wdir, "SL_RNA_genes"))
sls <- apply(summ, 1, function(su) {
	sl.name <- su["Name"]
	donatrons <- cons[[su["Sequence"]]][,c("Consensus", "Tail", "Read", "Donor_Region", "Chrom_Tail", "Pos_Tail", "Maj_Consensus")]
	donatrons <- subset(as.data.frame(donatrons), Maj_Consensus==TRUE)
	donatrons <- donatrons[rev(order(nchar(donatrons$Donor_Region))),]
	donatrons <- donatrons[order(donatrons$Chrom_Tail, donatrons$Pos_Tail),]
	donatrons <- donatrons[!duplicated(donatrons[,c("Chrom_Tail", "Pos_Tail")]),]
	donatrons$Strand <- regmatches(donatrons[,"Pos_Tail"], regexpr("[+-]", donatrons[,"Pos_Tail"]))
	donatrons$Stop <- as.numeric(gsub("\\(.\\)", "", donatrons$Pos_Tail))
	donatrons$Start=donatrons$Stop-(as.numeric(paste0(donatrons$Strand, "1"))*nchar(donatrons$Tail))
	
	# donatron FASTA
	fa <- lapply(1:nrow(donatrons), function(don){
		c(paste0(">", sl.name, ".", don), #" [", donatrons[don,"Chrom_Tail"], ":", donatrons[don,"Pos_Tail"], "]"),
		paste0(donatrons[don,"Donor_Region"]))
	})
	fa <- do.call(c, fa)
	writeLines(fa, file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".RNA_genes.fa")))
	# donatron GFF
	gff <- lapply(1:nrow(donatrons), function(don){
		# chrom, source, type. start, end, score, strand, phase, attributes 
		paste(donatrons[don,"Chrom_Tail"],
		"SLIDR", "spliced_leader_RNA", 
		max(1, min(donatrons[don,c("Start", "Stop")])+1), 
		max(donatrons[don,c("Start", "Stop")]), ".", donatrons[don,"Strand"], ".", 
		paste0("ID=", sl.name, ".", don),
		sep="\t")
		
		#c(paste0(">", sl.name, ".", don," [", donatrons[don,"Chrom_Tail"], ":", donatrons[don,"Pos_Tail"], "]"),
		#paste0(donatrons[don,"Consensus"], donatrons[don,"Donor_Region"]))
	})
	gff <- unname(do.call(c, gff))
	#gff <- gff[order(gff[,1], gff[,4]),]
	writeLines(c("##gff-version 3",	"# predicted using SLIDR 1.0", gff), file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".RNA_genes.gff3")))
	
	# write trans-spliced genes to separate GFF file
	tsg <- unique(cons[[su["Sequence"]]][,c("Chrom_Gene", "Pos_Gene")])
	tsg <- as.data.frame(tsg)
	tsg$Strand <- substr(tsg[,"Pos_Gene"], nchar(tsg[,"Pos_Gene"])-1, nchar(tsg[,"Pos_Gene"])-1)
	tsg.gff <- lapply(1:nrow(tsg), function(don){
		# chrom, source, type. start, end, score, strand, phase, attributes 
		pos <- unlist(strsplit(tsg[don, "Pos_Gene"], "[-(]"))[1:2]
		paste(tsg[don,"Chrom_Gene"],
		"SLIDR", "trans_splice_acceptor_site", 
		pos[1], pos[2], ".", tsg[don,"Strand"], ".", 
		paste0("ID=", sl.name, ".trans-splice-site-", don),
		sep="\t")
	})
	tsg.gff <- unname(do.call(c, tsg.gff))
	writeLines(c("##gff-version 3",	"# predicted using SLIDR 1.0", tsg.gff), file.path(wdir, "SL_RNA_genes", paste0(sl.name, ".trans-splice-sites.gff3")))
	
	# return FASTA consensus for SL
	# and the gff output
	# and the donatron FASTA
	return(list(c(paste0(">", sl.name), unique(cons[[su["Sequence"]]]$Consensus, "]")), gff, fa))
})
sl.fa <- do.call(c, lapply(sls, function(x) x[[1]]))
sl.gff <- do.call(c, lapply(sls, function(x) x[[2]]))
sl.don.fa <- do.call(c, lapply(sls, function(x) x[[3]]))
writeLines(sl.fa, file.path(wdir, "SLs.fa"))
writeLines(c("##gff-version 3",	"# predicted using SLIDR 1.0", sl.gff), file.path(wdir, "SL_RNA_genes.gff3"))
writeLines(sl.don.fa, file.path(wdir, "SL_RNA_genes.fa"))
