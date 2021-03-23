#!/usr/bin/env Rscript

options(scipen=999)
library(parallel)

# expect three arguments from command line:
# 1: working directory (containing featureCounts output)
# 2: number of threads
# 3: read length
# 4: minimum read threshold

args = commandArgs(trailingOnly=TRUE)
wdir <- args[1]
threads <- as.numeric(args[2])
rlen <- as.numeric(args[3])
thresh <- as.numeric(args[4])

# load and tidy up data from featureCounts
loadFC <- function(fc) {
  dd <- read.table(fc, header=T, stringsAsFactors = F, sep="\t")
  # tidy up sample names
  #colnames(dd) <- gsub(".R1.bam", "", gsub("^.*.LIB", "LIB", colnames(dd)))
  
  # tidy up chromosomes
  if(length(grep(";", dd$Chr))>0) dd$Chr <- unlist(sapply(strsplit(dd$Chr, ";"), function(x) paste(unique(x), collapse=";")))
  
  # tidy up strand info
  if(length(grep(";", dd$Strand))>0) dd$Strand <- unlist(sapply(strsplit(dd$Strand, ";"), function(x) paste(unique(x), collapse=";")))
  
  # tidy up start/end/length
  if(length(grep(";", dd$Start))>0) dd$Start <- sapply(strsplit(dd$Start, ";"), function(x) min(as.numeric(x)))
  if(length(grep(";", dd$End))>0) dd$End <- sapply(strsplit(dd$End, ";"), function(x) max(as.numeric(x)))
  dd$Length <- abs(dd$End-dd$Start)+1
  
  # remove any unresolved genes (annotation errors)
  dd <- dd[unique(grep(";", dd$Chr, invert=T), grep(";", dd$Strand, invert=T)),]
  
  # sort by start position (should already be sorted, but to be sure)
  dd <- dd[order(dd$Chr,dd$Start),]
  dd <- subset(dd, Geneid != "")
}

# split genes by exon peaks
peak_split <- function(dd){
	dd$Counts <- rowSums(dd[-c(1:6)])
	dd$Counts[dd$Counts<thresh] <- 0
	genes <- mclapply(split(dd, dd$Geneid), mc.cores = threads, FUN=function(gene){
		minus <- any(gene$Strand=="-")
		if(minus) gene <- gene[nrow(gene):1,]
		cts <- gene$Counts

		# new: need to allow single-exon splits if distance is sufficiently large
		# find exons with counts where the following exon also has counts
		o <- which(cts!=0 & c(cts[-1], 0)!=0)
		# if exon length is smaller than read length, delete counts of second exon
		del <- o[which(gene[o,"Length"]<rlen)]+1
		cts[del] <- 0
		gene[del,7:(ncol(gene)-1)] <- 0

		# get distance between 5' ends of consecutive exons
		# d5 <- abs(gene[o, 3+minus]-gene[o+1, 3+minus])
		# if distance is smaller than read length (150 bp), delete counts of second exon
		# cts[o[which(d5<150)]] <- 0
		
		# all count locations are peaks
		pos <- which(cts>0)
		
		# old: primary peaks
		#d1 <- diff(c(0,cts))
		#pr <- cts>0 & d1-cts==0
		# old: secondary peaks 
		#d2 <- c(0,0,pr[0:(length(pr)-2)])
		#d3 <- c(diff(cts),0)
		#se <- cts>0 & d2==1 & d3<=0
		#pos <- which(pr|se)
		
		# split at peak positions
		if(length(pos)>0 && any(pos>2)){	
			sp <- c(pos, nrow(gene)+1)
			if(pos[1]>1) sp <- c(1, sp)
			lbl <- rep(paste0("split", 1:(length(sp)-1)), times=c(diff(sp)))
			gene$Geneid <- paste(gene$Geneid, lbl, sep="_")
		}
		# aggregate by gene name
		gene <- do.call(rbind, lapply(split(gene, gene$Geneid), function(g){
				data.frame(Geneid=unique(g$Geneid),
						Chr=unique(g$Chr),
						Start=min(g$Start),
						End=max(g$End),
						Strand=unique(g$Strand),
						Length=sum(g$Length),
						rbind(colSums(g[-c(1:6)])))
		}))
		# revert if necessary
		if(minus) gene <- gene[nrow(gene):1,]
		gene[,-ncol(gene)]
	})
	genes <- do.call(rbind, genes)
	return(genes)
}

# background counts and gene-based SL counts need tidying only
for(fc in c("bg.featureCounts.genes.raw.txt", "un.featureCounts.genes.raw.txt", "SL.featureCounts.genes.raw.txt")){
	gbc <- loadFC(file.path(wdir, fc))
	cat(paste0(fc, ": ", nrow(gbc), " gene records comprising ", length(unique(gbc$Geneid)), " gene IDs ...\n"))
	write.table(gbc, file.path(wdir, gsub("raw", "clean", fc)), row.names=F, col.names=T, quote=F, sep="\t")
}

# exon-based counts need tidying and peak correction
ebc <- loadFC(file.path(wdir, "SL.featureCounts.exons.raw.txt"))
ng <- length(unique(ebc$Geneid))
cat(paste("Processed", nrow(ebc), "exon records comprising", ng, "gene IDs\n"))

ebc <- peak_split(ebc)
ngn <- length(unique(ebc$Geneid))
cat(paste0("Gained ", ngn-ng, " gene IDs (", ngn , " total) by splitting at internal exons with SL reads\n"))

write.table(ebc, file.path(wdir, "SL.featureCounts.exons.clean.txt"), row.names=F, col.names=T, quote=F, sep="\t")

