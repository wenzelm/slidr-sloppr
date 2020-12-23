library(parallel)
library(data.table)

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
		col.names=c("Donor_Region","Read","Centroid","Representative","Chrom_Tail","Pos_Tail",
					"BLASTN_Match","Outron_Overlap","Chrom_Gene","Pos_Gene","Acceptor_Region",
					"Loops","MFE_Frequency","Ensemble_Diversity"))
m$Outron_Overlap[m$Outron_Overlap=="-"] <- ""
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
#if(acc!="") {
	acc.len <- nchar(gsub("\\[[^][]*\\]", ".", acc))
	m$Acceptor.Splicesite <- substring(m$Acceptor_Region, nchar(m$Outron_Overlap)+1, nchar(m$Outron_Overlap)+acc.len)
	# additional sanity check: is the upstream outron sequence identical to Outron_Overlap from the tail?
	#clusterExport(cl=clust, varlist="acc.len")
	ovl <- substring(m$Acceptor_Region, acc.len+1, acc.len+nchar(m$Outron_Overlap))
	m$Acceptor.Splicesite[ovl!=m$Outron_Overlap] <- "#NA#"
	m <- m[grep(paste0("^", acc, "$"), m$Acceptor.Splicesite),]
	cat(paste0(length(unique(m$Read)), " reads with splice acceptor site (", acc, ")\n"))
#}

# 4) correct genomic positions using Outron_Overlap
# resolve 3' position of SL 
# Pos_Tail is the 3' position of the tail (5' position can be truncated depending on tail)
m$Tail <- paste0(m$BLASTN_Match, m$Outron_Overlap)
m$Pos_Tail <- apply(m, 1, function(rd){
	strand <- regmatches(rd["Pos_Tail"], regexpr("\\([\\+\\-]\\)", rd["Pos_Tail"]))
	pos <- unlist(strsplit(rd["Pos_Tail"], "[-(]"))[1:2]
	offs <- nchar(rd["Outron_Overlap"])
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, strand)
	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs, strand)
	#if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, "-", as.numeric(pos[1])+offs+2, strand)
	#if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs-2, "-", as.numeric(pos[2])-offs, strand)
	newpos
})

# need to shift 5' gene position by overlap
# new position is location of splice acceptor site 
m$Pos_Gene <- apply(m, 1, function(rd){
	strand <- regmatches(rd["Pos_Gene"], regexpr("\\([\\+\\-]\\)", rd["Pos_Gene"]))
	pos <- unlist(strsplit(rd["Pos_Gene"], "[-(]"))[1:2]
	offs <- nchar(rd["Outron_Overlap"])
	#if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, strand)
	#if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs, strand)
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs+1, "-", as.numeric(pos[1])+offs+length(acc)+1, strand)
	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs-length(acc), "-", as.numeric(pos[2])-offs, strand)
	newpos
})

# 5) keep reads where the SL gene is at least 1000bp away from trans-spliced gene
m$Distance <- apply(m, 1, function(rd){
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
m <- subset(m, Distance>1000)
cat(paste(length(unique(m$Read)), "reads with plausible trans-splicing pattern (>1 kbp between donor and acceptor site)\n"))

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
	SL_Genes=nrow(unique(tl[,c("Chrom_Tail", "Pos_Tail")])),
	Spliced_Genes=nrow(unique(tl[,c("Chrom_Gene", "Pos_Gene")])),
	Stem_Loops=paste(unique(tl$Loops), collapse=";"),
	MFE_Frequency=paste(unique(tl$MFE_Frequency), collapse=";"),
	Ensemble_Diversity=paste(unique(tl$Ensemble_Diversity), collapse=";"),
	stringsAsFactors=F)
})
summ <- do.call(rbind, summ)
summ <- summ[rev(order(summ$Spliced_Genes, summ$Reads)),]
write.table(summ, file.path(wdir, "3-candidate_SLs.txt"), row.names=F, col.names=T, quote=F, sep="\t")

summ <- subset(summ, Spliced_Genes>1 & Reads>1)
summ <- data.frame(Name=paste0(sl.prefix, 1:nrow(summ)), summ)
write.table(summ, file.path(wdir, "3-final_SLs.txt"), row.names=F, col.names=T, quote=F, sep="\t")

cat(paste(nrow(summ), "SLs with at least 2 reads and spliced to at least 2 genes\n")) 
cat("Top 10 SLs:\n")
print(head(summ[,c("Sequence", "Reads", "SL_Genes", "Spliced_Genes")], n=10), row.names=F)

# then, write donatron FASTA and GFF files
cat("Writing output files ...\n")
dir.create(file.path(wdir, "3-final_SL_donatrons"))
sls <- apply(summ, 1, function(su) {
	sl.name <- su["Name"]
	donatrons <- cons[[su["Sequence"]]][,c("Consensus", "Tail", "Read", "Donor_Region", "Chrom_Tail", "Pos_Tail")]
	donatrons <- as.data.frame(donatrons)
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
	writeLines(fa, file.path(wdir, "3-final_SL_donatrons", paste0(sl.name, ".donatrons.fa")))
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
	writeLines(c("##gff-version 3",	"# predicted using SLIDR 1.0", gff), file.path(wdir, "3-final_SL_donatrons", paste0(sl.name, ".donatrons.gff3")))
	
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
	writeLines(c("##gff-version 3",	"# predicted using SLIDR 1.0", tsg.gff), file.path(wdir, "3-final_SL_donatrons", paste0(sl.name, ".trans-splice-sites.gff3")))
	
	# return FASTA consensus for SL
	# and the gff output
	# and the donatron FASTA
	return(list(c(paste0(">", sl.name), unique(cons[[su["Sequence"]]]$Consensus, "]")), gff, fa))
})
sl.fa <- do.call(c, lapply(sls, function(x) x[[1]]))
sl.gff <- do.call(c, lapply(sls, function(x) x[[2]]))
sl.don.fa <- do.call(c, lapply(sls, function(x) x[[3]]))
writeLines(sl.fa, file.path(wdir, "3-final_SLs.fa"))
writeLines(c("##gff-version 3",	"# predicted using SLIDR 1.0", sl.gff), file.path(wdir, "3-final_SLs_donatrons.gff3"))
writeLines(sl.don.fa, file.path(wdir, "3-final_SLs_donatrons.fa"))

q()


###################################################################################
## old code


old_fread {
# 1) Read depth for each representative read
dr <- fread(cmd=paste0("gunzip -c ", file.path(wdir, "1-tails", "tails.derep.clusterinfo.txt")), sep="\t", header=F,
					colClasses="character",	
					col.names=c("Representative", "Read"), showProgress=F)
cn <- fread(cmd=paste0("gunzip -c ", file.path(wdir, "1-tails", "clusters.clusterinfo.txt")), sep="\t", header=F,
					colClasses="character",	
					col.names=c("Centroid", "Representative"), showProgress=F)
dr <- merge(unique(dr), unique(cn), by="Representative")
rm(cn)					
cat(paste(length(unique(dr$Read)), "raw reads\n"))
					
# 2) add donor/Sm data
don <- fread(cmd=paste0("gunzip -c ", file.path(wdir, "2-RNA_filters", "splice_donor_sites.txt.gz")), sep="\t", header=F, showProgress=F, 
					colClasses="character",	
					col.names=c("Centroid", "Chrom", "Pos", "BLASTN_Match", "Outron_Overlap", "Donor_Region"))
don$Outron_Overlap[don$Outron_Overlap=="."] <- ""
m <- merge(dr, unique(don), by="Centroid", sort=F, allow.cartesian=TRUE)
rm(don, dr)
cat(paste(length(unique(m$Read)), "reads with evidence of splice donor and/or Sm sites\n"))

# 4) add RNAfold data
rn <- fread(cmd=paste0("gunzip -c ", file.path(wdir, "2-RNA_filters", "RNAfold.txt.gz")), sep="\t", header=F, showProgress=F, 
					colClasses="character",	
					col.names=c("Donor_Region", "Loops", "MFE_Frequency", "Ensemble_Diversity"))
rn$Loops <- sapply(gregexpr("\\(\\.+\\)", rn$Loops), function(x) length(attr(x, "match.length")))

m <- merge(m, unique(rn), by="Donor_Region", sort=F, allow.cartesian=TRUE)
rm(rn)
cat(paste(length(unique(m$Read)), "reads with RNAfold results\n"))

# 5) add cluster data
hq <- fread(cmd=paste0("gunzip -c ", file.path(wdir, "2-RNA_filters", "SL_clusters.clusterinfo.txt.gz")), sep="\t", header=F, showProgress=F,
					colClasses="character",	
					col.names=c("Cluster", "Centroid", "BLASTN_Match", "Outron_Overlap"))
m <- merge(unique(hq), m, by=c("Centroid", "BLASTN_Match", "Outron_Overlap"), sort=F, allow.cartesian=TRUE)
tmp <- merge(unique(hq), tmp, by=c("Centroid", "BLASTN_Match", "Outron_Overlap"), sort=F)
rm(hq)
#cat(paste(length(unique(m$Cluster)), "clusters\n"))

# 6) add acceptor data
#if(file.exists(file.path(wdir, "SL_RNA_filters", "splice_acceptor_sites.txt"))) {
ag <- fread(cmd=paste0("gunzip -c ", file.path(wdir, "2-RNA_filters", "splice_acceptor_sites.txt.gz")), sep="\t", header=F,
				colClasses="character",	
				col.names=c("Read", "Chrom", "Pos", "Acceptor_Region"), showProgress=T)
ag$Read <- gsub(".{2}$", "", ag$Read)
m <- merge(m, unique(ag), by="Read", sort=F, suffixes=c("_Tail", "_Gene"), allow.cartesian=TRUE)
rm(ag)
cat(paste(length(unique(m$Read)), "reads with splice acceptor data\n"))
}


old_dump {
#m$Outron <- apply(m, 1, 
	#function(x) substring(x["Acceptor_Region"], acc.len+1, acc.len+nchar(x["Outron_Overlap"])))
#m$Acceptor.Splicesite <- apply(m, 1, 
	#function(x) substring(x["Acceptor_Region"], nchar(x["Outron_Overlap"])+1, nchar(x["Outron_Overlap"])+acc.len))

#m$Acceptor.Splicesite <- apply(m, 1, 
#	function(x) {
#		ss <- substring(x["Acceptor_Region"], nchar(x["Outron_Overlap"])+1, nchar(x["Outron_Overlap"])+acc.len)
#		if(substring(x["Acceptor_Region"], 1, nchar(x["Outron_Overlap"]))==x["Outron_Overlap"]) {
#			ss <- 
#		} else {
#			ss <- ""
#		}
#		return(ss)
#	})

#m <- subset(m, Acceptor=="AG")
#}

#cat("Constructing consensus SLs ...\n")

# resolve 3' position of SL 
# Pos_Tail is the 3' position of the tail
m$Tail <- paste0(m$BLASTN_Match, m$Outron_Overlap)
m$Pos_Tail <- apply(m, 1, function(rd){
	strand <- regmatches(rd["Pos_Tail"], regexpr("\\([\\+\\-]\\)", rd["Pos_Tail"]))
	pos <- unlist(strsplit(rd["Pos_Tail"], "[-(]"))[1:2]
	offs <- nchar(rd["Outron_Overlap"])
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, strand)
	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs, strand)
	#if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, "-", as.numeric(pos[1])+offs+2, strand)
	#if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs-2, "-", as.numeric(pos[2])-offs, strand)
	newpos
})

# need to shift 5' gene position by overlap
m$Pos_Gene <- apply(m, 1, function(rd){
	strand <- regmatches(rd["Pos_Gene"], regexpr("\\([\\+\\-]\\)", rd["Pos_Gene"]))
	pos <- unlist(strsplit(rd["Pos_Gene"], "[-(]"))[1:2]
	offs <- nchar(rd["Outron_Overlap"])
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, strand)
	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs, strand)
	#if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, "-", pos[2], strand)
	#if(strand=="(-)") newpos <- paste0(pos[1], "-", as.numeric(pos[2])-offs, strand)
	newpos
})

# need to shift 3' tail position by overlap
#m$Pos_Tail <- apply(m, 1, function(rd){
#	strand <- regmatches(rd["Pos_Tail"], regexpr("\\([\\+\\-]\\)", rd["Pos_Tail"]))
#	pos <- unlist(strsplit(rd["Pos_Tail"], "[-(]"))[1:2]
#	offs <- nchar(rd["Outron_Overlap"])
#	if(strand=="(+)") newpos <- paste0(pos[1], "-", as.numeric(pos[1])+offs, strand)
#	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs, "-", as.numeric(pos[2]), strand)
#	#if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, "-", as.numeric(pos[1])+offs+2, strand)
#	#if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs-2, "-", as.numeric(pos[2])-offs, strand)
#	newpos
#})
# keep reads where the SL gene is at least 500bp away from trans-spliced gene
m$Distance <- apply(m, 1, function(rd){
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
m <- subset(m, Distance>500)

cat(paste(length(unique(m$Read)), "reads with plausible trans-splicing pattern\n"))
}






for(id in summ$Name){
	#sl.name <- paste0("SL", id)
	#header <- paste0(">"sl.name, " [", unique(cons[[id]]$Consensus), "]")
	#print(head(cons[[id]]))
	donatrons <- unique(cons[[id]][,c("Consensus", "Donor_Region", "Chrom_Tail", "Pos_Tail")])
	fa <- lapply(1:nrow(donatrons), function(don){
		c(paste0(">", sl.name, ".", don," [", donatrons[don,"Chrom_Tail"], ":", donatrons[don,"Pos_Tail"], "]"),
		paste0(donatrons[don,"Consensus"], donatrons[don,"Donor_Region"]))
	})
	fa <- do.call(c, fa)
	writeLines(fa, file.path(wdir, "SL_donatrons", paste0(sl.name, ".donatrons.fa")))
	#print(nrow(cons[[id]][,c("Consensus", "Donor_Region", "Chrom_Tail", "Pos_Tail")]))
	#print(nrow(unique(cons[[id]][,c("Consensus", "Donor_Region", "Chrom_Tail", "Pos_Tail")])))
	#paste0(cons[[id]]$Consensus, cons[[id]]$Donor)
}

summ_tails <- function(x){

}
summ_tails(m)
#head(summ_tails(subset(m, Loops==2)))

write.table(summ_tails(m), file.path(wdir, "final_SLs.txt"), quote=F, row.names=T, col.names=NA, sep="\t")
write.table(summ_tails(subset(m, Loops==2)), file.path(wdir, "final_SLs_2loops.txt"), quote=F, row.names=T, col.names=NA, sep="\t")

# FASTA files for the SL and the full donatron (SL+Donor)



 
 q()


#####################

# read the following files:
# File that links representative reads to duplicate reads
# File that links cluster IDs to clustered representative reads
# File containing GT/Sm sites for each representative read
# File containing AG sites for each read alignment with soft-masked tail
# File containing RNAfold output for each GT/Sm sequence

dr <- read.table(file.path(wdir, "tails", "tails.derep.clusterinfo.txt"), stringsAsFactors=F, sep="\t",
					col.names=c("Representative", "Reads"))
#cl <- read.table(file.path(wdir, "02_LQ_clusters.fa.clusterinfo"), stringsAsFactors=F, sep="\t",
#					col.names=c("Cluster", "Representative"))
hq <- read.table(file.path(wdir, "SL_RNA_filters", "SL_clusters.txt"), stringsAsFactors=F, sep="\t",
					col.names=c("Cluster", "Representative", "Match", "Overlap"))

sl <- read.table(file.path(wdir, "SL_RNA_filters", "GT_splice_Sm_sites.txt"), stringsAsFactors=F, sep="\t",
					col.names=c("Representative", "Chrom", "Pos", "Match", "Overlap", "Donor"))

#sl <- read.table(text=gsub("::", "\t", readLines(file.path(wdir, "03_GT_splice_Sm_sites.txt"))), stringsAsFactors=F, sep="\t",
#					col.names=c("Representative", "Chrom", "Pos", "Match", "Overlap", "Donor"))
#ag <- read.table(text=gsub("::", "\t", readLines(file.path(wdir, "03_AG_splice_sites.txt"))), stringsAsFactors=F, sep="\t",
#					col.names=c("Read", "Chrom", "Pos", "Acceptor"))
# RNAfold output

# Merge and filter datasets
m <- merge(unique(dr), unique(sl), by="Representative", sort=F)
m <- merge(m, unique(ag), by="Read", sort=F, suffixes=c(".tail", ".gene"))
# resolve overlap between 3' tail end and AG splice acceptor site
m$Acceptor <- apply(m, 1, function(x) substring(x["Acceptor"], nchar(x["Overlap"])+1, nchar(x["Overlap"])+2))
m$Tail <- paste0(m$Match, m$Overlap)
m <- subset(m, Acceptor=="AG")
nrow(m)
# need to shift gene position by overlap
m$Pos.gene <- apply(m, 1, function(rd){
	strand <- regmatches(rd["Pos.gene"], regexpr("\\([\\+\\-]\\)", rd["Pos.gene"]))
	pos <- unlist(strsplit(rd["Pos.gene"], "[-(]"))[1:2]
	offs <- nchar(rd["Overlap"])
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, "-", pos[2], strand)
	if(strand=="(-)") newpos <- paste0(pos[1], "-", as.numeric(pos[2])-offs, strand)
	newpos
})
m$Pos.tail <- apply(m, 1, function(rd){
	strand <- regmatches(rd["Pos.tail"], regexpr("\\([\\+\\-]\\)", rd["Pos.tail"]))
	pos <- unlist(strsplit(rd["Pos.tail"], "[-(]"))[1:2]
	offs <- nchar(rd["Overlap"])
	if(strand=="(+)") newpos <- paste0(as.numeric(pos[1])+offs, "-", as.numeric(pos[1])+offs+2, strand)
	if(strand=="(-)") newpos <- paste0(as.numeric(pos[2])-offs-2, "-", as.numeric(pos[2])-offs, strand)
	newpos
})
# add RNAfold and HQ clusters
m <- merge(m, unique(rr), by="Donor", sort=F)
m <- merge(unique(hq), m, by=c("Representative", "Match", "Overlap"), sort=F)

# filter cis-splicing patterns
# keep reads where the SL gene is at least 500bp away from trans-spliced gene
m$Distance <- apply(m, 1, function(rd){
	if(rd["Chrom.tail"] != rd["Chrom.gene"]) {	# different chromosomes; easy
		di <- Inf
	} else {	# same chromosomes
		# check strands
		s1 <- regmatches(rd["Pos.tail"], regexpr("\\([\\+\\-]\\)", rd["Pos.tail"]))
		s2 <- regmatches(rd["Pos.gene"], regexpr("\\([\\+\\-]\\)", rd["Pos.gene"]))
		if(s1 != s2) {	# different strands
			di <- Inf
		} else {	# same strands
			p1 <- unlist(strsplit(rd["Pos.tail"], "[-(]"))
			p2 <- unlist(strsplit(rd["Pos.gene"], "[-(]"))
			if(s1=="(+)") di <- abs(as.numeric(p1[1])-as.numeric(p2[1]))
			if(s1=="(-)") di <- abs(as.numeric(p1[2])-as.numeric(p2[2]))
		}
	}
	di
})
m <- subset(m, Distance>500)
nrow(m)

# resolve cluster consensus
m <- mclapply(split(m, m$Cluster), mc.cores=function(cls){
	cls$Consensus <- rep(cls$Tail[which.max(nchar(cls$Tail))], times=nrow(cls))
	tail.align <- paste(sapply(max(nchar(cls$Tail))-nchar(cls$Tail), 
		function(x) paste(rep("-", times=x), collapse="")), 
		cls$Tail, sep="")
	cls$Consensus2 <- paste(apply(do.call(rbind, strsplit(tail.align, "")), 2, 
		function(x) names(sort(table(x[x!="-"]), decreasing=T))[1]), collapse="")
	return(cls)
})
m <- do.call(rbind, m)



#o <- unique(subset(m, Consensus2=="TCGGTTTAATTACCCAAGTTTGAG")[,c("Chrom.tail", "Pos.tail", "Match", "Overlap")])
#unique(o[order(o$Pos.tail),])


#table(m$Consensus, m$Cluster)

#m$Consensus <- unlist(sapply(split(m$Tail, m$Cluster), function(x) rep(x[which.max(nchar(x))], times=length(x))))

# play around with positive control (Cel SL)


n <- m[grep("ggtttaattacccaagtttgag|ggttttaacccagttactcaag", m$Consensus, ignore.case=T),]
table(n$Tail, n$Cluster)
table(n$Consensus, n$Cluster)

# Check if all Consensuses are unique to a single cluster (if not, there is a clustering problem)
table(apply(table(m$Consensus, m$Cluster), 1, function(x) sum(x>0)))
table(apply(table(m$Consensus, m$Cluster), 2, function(x) sum(x>0)))

# for each full tail, we want to know:
# - how many reads support the tail
# - how many genes the tail is spliced to

summ_tails <- function(x){
	x <- mclapply(split(x, x$Consensus), mc.cores=16, function(tl){
		data.frame(Length=nchar(tl$Consensus)[1],
		Reads=length(unique(tl$Read)),
		Unique_Reads=length(unique(tl$Representative)),
		Unique_Tails=length(unique(tl$Tail)),
		Tail_Genes=nrow(unique(tl[,c("Chrom.tail", "Pos.tail")])),
		Spliced_Genes=nrow(unique(tl[,c("Chrom.gene", "Pos.gene")])),
		Loops=paste(unique(tl$Loops), collapse=";"),
		MFE_Frequency=paste(unique(tl$MFE_Frequency), collapse=";"),
		Ensemble_Diversity=paste(unique(tl$Ensemble_Diversity), collapse=";")
		)
	})
	ss <- subset(do.call(rbind, x), Spliced_Genes>1 & Reads>1)
	ss[rev(order(ss$Reads, ss$Spliced_Genes)),]
}
head(summ_tails(m))
head(summ_tails(subset(m, Loops==2)))

write.table(summ_tails(m), file.path(wdir, "04_final_tails.txt"), quote=F, row.names=T, col.names=NA, sep="\t")
write.table(summ_tails(subset(m, Loops==2)), file.path(wdir, "04_final_tails_loops.txt"), quote=F, row.names=T, col.names=NA, sep="\t")
 
q()



# Merge all datafiles into master dataframe
m <- merge(sl, rr, by=c("Representative", "Donor"), sort=F)
m <- merge(unique(dr), m, by="Representative", sort=F)
m <- merge(unique(cl), m, by="Representative", sort=F)
m <- merge(unique(hq), m, by="Representative", sort=F)
m <- merge(m, ag, by="Read", sort=F); nrow(m)



# resolve consensus tail within each cluster
# any 3' truncation should already be resolved
# but 5' truncation may still be present:
# AACCGGTTAACCGGTT
# ---CGGTTAACCGGTT
# since the tails should all line up at 3' end,
# we can simply retain the longest read as a consensus.




n$Consensus <- do.call(c, sapply(split(n$Fulltail, n$Cluster), function(x)   rep(x[which.max(nchar(x))], times=length(x))))






q()

##########################








summ <- mclapply(split(bl, bl[,1]), function(x){ 
	# best hit(s)
	#x <- split(bl, bl[,1])[[1]]
	tl <- x[x[,12]==max(x[,12]),]
	# cluster id
	clid <- gsub(";", "", gsub(".*clusterid=", "", tl[1,1]))
	# find all reads in cluster
	clreads <- subset(cl, V1==clid)
	readlocs <- unique(clreads[,3:4])
	# compare each tail location against each read location
	comps <- apply(tl, 1, function(tloc){ 
		comp <- apply(readlocs, 1, function(rloc){
			if(tloc[2]==rloc[1]) {	# same contig
				distance <- as.numeric(tloc[9])-as.numeric(rloc[2])
			} else {	# different contig
				distance <- NA
			}
			data.frame(tloc[2], rloc[1], distance)
		})
		do.call(rbind, comp)
	})
	comps <- do.call(rbind, comps)
	data.frame(Tail=tl[1,13], Tail.Locations=nrow(tl), 
				Same.contig=paste(round(100*sum(!is.na(comps[,3]))/nrow(comps), digits=2), "%"), 
				Min.Distance=min(abs(comps[,3]), na.rm=T), 
				Max.Distance=max(abs(comps[,3]), na.rm=T))
}, mc.cores = 16)

summ <- do.call(rbind, summ)
summ[,4][is.infinite(summ[,4])] <- NA
summ[,5][is.infinite(summ[,5])] <- NA
write.table(summ, file.path(wdir, "candidates.txt"), quote=F, sep="\t", row.names=T, col.names=NA)

#summ[grep("=4351;|=4627;|=5033;|=5085;|=5212;|=5499;|=4607;|=4726;|=4987;|=5498;", rownames(summ)),]

