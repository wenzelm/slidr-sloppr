#!/usr/bin/env Rscript

sloppr_version <- "1.1.5"

# this script is designed to predict operons from a featureCount matrix

options(scipen=999)
library("parallel")
library("ggdendro")
library("MASS")
library("ggplot2")
library("grid")
library("gridExtra")
library("scales")
library("reshape2")
library("glmpca")

# Arguments:	1) featureCounts output file
#				2) SL FASTA file
#				3) text file with SL2-type SL names
#				4) SL2:SL1 read-ratio threshold to classify gene as downstream (default: infinity)
#				5) SL2-bias at upstream operonic genes (yes/no)
#				6) intercistronic distance cutoff (default: infinity)
#				7) Operon prefix (default: OP)
#				8) aggregation function (sum|median|geo)
#				9) method for dealing with zeros (default:ignore)
#				10) number of threads
#				11) output directory

args = commandArgs(trailingOnly=TRUE)
#args <- c("sloppr_toy_data//2-counts/SL.featureCounts.genes.clean.txt", 
#			"SL.fasta", 
#			"SL2.txt", 
#		"infinity", 
#			"no",
#			"infinity",
#			"COP", 
#			"geometric mean",
#			"remove",
#			"1", ".")
#print(cbind(args))
f <- args[1]
s <- args[2]
s2 <- gsub("^ ", "", args[3])
rr.thresh <- as.numeric(args[4])
upstream.bias <- args[5]
dist.cutoff <- args[6]
op.prefix <- args[7]
agg <- args[8]
zero <- args[9]
threads <- as.numeric(args[10])
wdir <- args[11]

# colours for plotting
colours.cl <- c(`1`="gray30", `2`="#FF9100", 
				Cluster1="gray30", Cluster2="#FFA500", `Cluster1+Cluster2`="red2",
				SL1="gray30", SL2="#FFA500", `SL1+SL2`="red2", `no SL`="lightgray", `not expressed`="steelblue2",
				operonic="#FFA500", `not operonic`="gray30")

# outfile prefix
outpref <- gsub(".clean.txt", "", basename(f))

# redirect warnings and messages to a log file
sink(file(file.path(wdir, paste0(outpref, ".log.txt")), open = "wt"), type = "message")

cat(paste0("Processing counts data from ", f, " ...\n"))

#
# 1) load data
#
cleanSLnames <- function(x){
	x <- gsub("^>", "", x)
	# need to remove special characters
	x.tmp <- make.names(x)  
	# in some cases an "X" is added to the beginning; remove if necessary
	x <- substr(x.tmp, nchar(x.tmp)-nchar(x)+1, nchar(x.tmp))
	return(x)
}

# SL
sl <- readLines(s)
sl <- cleanSLnames(sl[grep("^>", sl)])

# a priori SL2-type SLs
if(file.exists(s2)) {
	sl2 <- cleanSLnames(readLines(s2)) 
} else {
	sl2 <- sl
}
cat(paste("SL1-type SLs:\n\t\t"))
cat(paste(setdiff(sl, sl2), collapse="\n\t\t"))
cat(paste("\nSL2-type SLs:\n\t\t"))
cat(paste(sl2, collapse="\n\t\t"))
cat("\n")

# featureCounts results
readFC <- function(f){
  d <- read.table(f, sep="\t", header=T, stringsAsFactors=F, comment.char="")
  # simplify
  #d <- d[,grep("alternative", colnames(d), invert=T)]
  
  # extract meta information
  meta <- data.frame(ID=colnames(d)[-c(1:6)], stringsAsFactors=F)
  # resolve library names
  meta$Library <- gsub(".*library_", "", gsub(".fc.bam$", "", meta$ID))
  for(s in sl) meta$Library <- gsub(paste0(".", s, "$"), "", meta$Library)
  # resolve SL names
  meta$SL <- rep("", times=nrow(meta))
  for(s in sl) meta$SL[grep(paste0(s, "\\."), meta$ID)] <- s
  # add SL-types (if known)
  meta$SL.type <- "SL1"
  meta$SL.type[meta$SL %in% sl2] <- "SL2"
  # add raw read counts
  meta$Reads <- colSums(d[-c(1:6)])
  
  #cat("\nIdentified the following samples and SLs:\n")
  #print(meta[-1])
  
  # sort genes and calculate distances
  d <- d[order(d$Chr,d$Start),]
  d <- mclapply(split(d, interaction(d$Chr, d$Strand), drop=T), function(chrs){
    if(nrow(chrs)<2) {
      chrs$Distance <- NA
    } else {
      if(any(chrs$Strand=="+")){
        chrs$Distance <- c(sapply(1:(length(chrs$Start)-1), function(x) { chrs$Start[x+1]-chrs$End[x] }), NA)
      } else {
        #chrs <- chrs[nrow(chrs):1,]
        chrs$Distance <- c(NA, sapply(2:length(chrs$Start), function(x) { chrs$Start[x]-chrs$End[x-1] }))
      }
    }
    return(chrs)
  }, mc.cores=threads)
  d <- do.call(rbind.data.frame, d)
  d <- d[order(d$Chr,d$Start),]
  d <- cbind(d[,c(1:6)], Distance=d$Distance, d[-c(1:6, ncol(d))])
  
  # compute CPM for SL counts from raw background counts
  # compute TPM for background counts
  
  # obtain background gene expression levels
  bg <- read.table(gsub("SL.featureCounts.genes|SL.featureCounts.exons", "bg.featureCounts.genes", f), sep="\t", header=T, stringsAsFactors=F)
  names(bg) <- gsub(".*library_", "", gsub(".end2end_pre.align.bam$", "", names(bg)))
  un <- read.table(gsub("SL.featureCounts.genes|SL.featureCounts.exons", "un.featureCounts.genes", f), sep="\t", header=T, stringsAsFactors=F)
  names(un) <- gsub(".*library_", "", gsub(".untrimmed.bam$", "", names(un)))
  # add un to bg counts
  un <- un[match(bg$Geneid, un$Geneid),]
  bg[-c(1:6)] <- bg[-c(1:6)] + un[-c(1:6)]
  
  # synchronise both count tables (necessary when using exon-corrected data)
  bg <- bg[match(gsub("_split.*", "", d$Geneid), bg$Geneid),]
  bg$Geneid <- d$Geneid
  
  d.cpm <- d
  bg.counts <- bg[-c(2:6)]
  for(lib in unique(meta$Library)){
	# add SL and bg counts
	sl.counts <- as.matrix(d[subset(meta, Library==lib)$ID])
	bg.counts[,lib] <- bg.counts[,lib] + rowSums(sl.counts)
	# CPM
	total <- sum(bg.counts[,lib]) #+ sum(sl.counts)
	d.cpm[subset(meta, Library==lib)$ID] <- sl.counts*1e6/total
	# TPM
	bg.counts[,lib] <- bg.counts[,lib]/bg$Length
	bg.counts[,lib] <- bg.counts[,lib]/(sum(bg.counts[,lib])/1e6)
  }
    
  bg.counts <- bg.counts[match(d$Geneid, bg.counts$Geneid),]
  return(list(Counts=d, CPM=d.cpm, Meta=meta, bg=bg.counts))
}
fc <- readFC(f)

#
# 1) infer SL clusters
# even if known SL2-type SLs were supplied.
#
sl_clusters <- function(dd){
  # extract counts
  meta <- dd$Meta
  cts <- dd$Counts[-c(1:7)]
  
  # filter
  #cts[cts<thresh] <- 0
  keep <- colSums(cts)>0
  cts <- cts[keep]  
  meta <- meta[keep,]

  # MDS/K-means/LDA
  # catch error here
  gpca <- glmpca(cts[rowSums(cts)>0,,drop=F],2,fam="nb")
  names(gpca$factors) <- c("X1", "X2")
  
  if(nrow(gpca$factors)>=2) {
	hc <- hclust(dist(gpca$factors), method="ward.D2")
  } else {
	hc <- hclust(dist(rbind(gpca$factors, gpca$factors)), method="ward.D2") #dummy
  }
  
  # write dendrogram groups to file
  write.table(cutree(hc, k=2:nrow(hc$merge)), file.path(wdir, paste0(outpref, ".SL_dendrogram.groups.txt")), row.names=T, col.names=NA, quote=F, sep="\t")
  
  #di <- dist(t(cts), method="binary")
  #mds <- cmdscale(di)
  #hc <- hclust(di, method="ward.D2")
  
  hc <- dendro_data(hc, type="rectangle")
  km <- kmeans(gpca$factors, centers=2-(nrow(gpca$factors)<=2), nstart=25, iter.max=25)

  if(nrow(km$centers)==2) {
	cat("successful!\n")
	ld <- predict(lda(gpca$factors, km$cluster, tol=1e-25))$x[,1]
  } else {
	cat("too few samples to define clusters!\n")
	ld <- 0
  }
  
  # write results file 
  meta <- data.frame(meta, SL.cluster=factor(km$cluster), gpca$factors, LD1=ld)
  write.table(meta, file.path(wdir, paste0(outpref, ".SL_clusters.txt")), row.names=F, col.names=T, sep="\t", quote=F)
  
  # generate plots  
  meta <- merge(hc$labels, meta, by.x="label", by.y="ID")
  meta <- meta[order(meta$SL.cluster),]
  meta$label <- paste0("    ", meta$Library, "-", gsub("splice.*", "", gsub("_", "", meta$SL)), " (", meta$Reads, " reads)")
  meta$Index <- 1:nrow(meta)
  # revert MDS and LD axes if necessary, such that Cluster1 is on the left
  if(order(tapply(meta$X1, meta$SL.cluster, mean))[1] == 2) meta$X1 <- meta$X1*(-1)
  if(order(tapply(meta$LD1, meta$SL.cluster, mean))[1] == 2) meta$LD1 <- meta$LD1*(-1)
  
  # plot1: dendrogram
  seg <- segment(hc)
  lbl.y <- max(seg$y)/20
  gg.tree <- ggplot() +
    geom_segment(data=seg, aes(x=x, y=y, xend=xend, yend=yend)) +
	# first circle/square without label
	# then label
    geom_point(data=meta, aes(x=x, color=SL.cluster, shape=SL.type), y=lbl.y, size=6) +
    geom_text(data=meta, aes(x=x, label=Index), size=3, y=lbl.y, 
		color="white", fontface=2, vjust=0.5) +
    geom_text(data=meta, aes(x=x, label=label, hjust=0, color=SL.cluster), y=lbl.y, size=3, vjust=0.5) +
	#geom_label(data=meta, aes(x=x, label=Index, fill=SL.cluster), hjust=0.5, size=3, y=lbl.y/5, color="white", fontface=2, vjust=0.5,
	#	label.padding = unit(0.25, "lines"),  label.r = unit(0, "lines") ) +
    scale_colour_manual(values=colours.cl[c("1", "2")]) +
	scale_fill_manual(values=colours.cl[c("1", "2")]) +
    scale_shape_manual(values=c(19, 15)) +
    coord_flip() + 
    scale_y_reverse(expand = expansion(mult = c(0.02, 1.5))) +
    theme_void() +
    theme(legend.position="right")
  
  # plot2: MDS with convex hull
  hull.data <- do.call(rbind, lapply(split(meta, meta$SL.cluster), function(x) x[chull(x[,c("X1", "X2")]),]))
  gg.mds <- ggplot(meta, aes(x=X1, y=X2, fill=SL.cluster, colour=SL.cluster)) +
    geom_point(aes(shape=SL.type), size=8, colour="white") +
    geom_polygon(data=hull.data, alpha=0.2) +
    geom_text(aes(label=Index), color="white", size=3, fontface=2, vjust=0.5) +
	#geom_label(aes(label=Index, fill=SL.cluster), color="white", size=3, fontface=2, vjust=0.5,
	#		label.padding = unit(0.25, "lines"),  label.r = unit(0, "lines") ) +
    scale_x_continuous(labels=function(x) sprintf("%.2f", x)) +
    scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
    scale_colour_manual(values=colours.cl) +
    scale_shape_manual(values=c(21, 22)) +
    scale_fill_manual(values=colours.cl) +
    xlab("Axis1") + ylab("Axis2") +
    theme_bw() +
    theme(legend.position="none")
	
  
  # plot3: LDA density
  gg.lda <- ggplot(meta, aes(x=LD1, fill=SL.cluster, colour=SL.cluster)) +
    geom_density(colour=NA) +
    geom_rug() +
    scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
    scale_colour_manual(values=colours.cl) +
    scale_fill_manual(values=colours.cl) +
    xlab("LD1") + ylab("Density") +
    theme_bw() +
    theme(legend.position="none")
  
  sl.gg <- arrangeGrob(gg.tree, gg.mds, gg.lda, layout_matrix=matrix(c(1,1,2,3), ncol=2), widths=c(0.5,0.5))
  dd$Meta[keep, "SL.cluster"] <- paste0("Cluster", km$cluster)
  dd <- c(dd, MDSplots=list(sl.gg))
  return(dd)
}
cat("Identifying SL clusters ... ")
fc <- sl_clusters(fc)
pdf(file.path(wdir, paste0(outpref, ".SL_clusters.pdf")), width=12, height=12)
	grid.draw(fc$MDSplots)
invisible(dev.off())

# write huge dendrogram
pdf(file.path(wdir, paste0(outpref, ".SL_dendrogram.pdf")), width=14, height=72)
	grid.draw(fc$MDSplots$grobs[[1]])
invisible(dev.off())

#
# summarise SL1/SL2 proportions for each library (for diagnosis).
sl_proportions <- function(ty, sls){
	ddf <- lapply(c("*", unique(fc$Meta$Library)), function(lib){
		cts <- lapply(sls, function(s) {
			rowSums(as.matrix(fc$CPM[,fc$Meta$ID[which(grepl(lib, fc$Meta$Library) & fc$Meta[,ty]==s)]]))
		})
		cts <- do.call(cbind.data.frame, cts)
		colnames(cts) <- sls
		cts$Status <- "no SL"
		cts$Status[rowSums(cbind(fc$bg[,grep(lib,colnames(fc$bg)[-1])+1]))==0] <- "not expressed"
		cts$Status[cts[,sls[1]]>0] <- sls[1]
		cts$Status[cts[,sls[2]]>0] <- sls[2]
		cts$Status[cts[,sls[1]]>0 & cts[,sls[2]]>0] <- paste(sls, collapse="+")
		cts$Status <- factor(cts$Status, c("not expressed", "no SL", sls[1], sls[2], paste(sls, collapse="+")))
		prop <- table(cts$Status)
		prop <- data.frame(prop, Percent=paste(round(100*prop/sum(prop), 2), "%"))
		data.frame(Type=ty, Library=gsub("^\\*$", "Total", lib), prop)
	})
	do.call(rbind, ddf)	
}
ddf <- rbind(sl_proportions("SL.type", c("SL1", "SL2")),
			sl_proportions("SL.cluster", c("Cluster1", "Cluster2")))
cat("\nNumbers of genes receiving SLs:")
print(unname(subset(ddf, Library=="Total")[-c(6,7),3:5]), row.names=F)
write.table(ddf, file.path(wdir, paste0(outpref, ".SL_genes.txt")), row.names=F, col.names=T, sep="\t", quote=F)

# plot
ddf$Var1 <- factor(ddf$Var1, c("not expressed", "no SL", "SL1", "Cluster1", "SL1+SL2", "Cluster1+Cluster2", "SL2", "Cluster2"))
pdf(file.path(wdir, paste0(outpref, ".SL_genes.pdf")), width=10, height=12)
	gg.prop <- ggplot(ddf, aes(x=Library, y=Freq, fill=Var1)) + 
		geom_bar(position="fill", stat="identity") + 
		facet_wrap(~Type) +
		scale_fill_manual(name="Gene status", values=colours.cl[levels(ddf$Var1)]) +
		scale_y_continuous(labels=percent) +
		xlab("Library") + ylab("Percent") +
		coord_flip() + 
		theme_classic() +
		theme(legend.position="top")
	gg.prop
invisible(dev.off())

#
# 2) resolve biological replicates for each SL type
#	 compute SL2:SL1 read ratios
#
resolve_reps <- function(ty, sls){
	# sum CPM within each library across SLs
	cpm.sum <- lapply(sls, function(sl){
		cpm.sum <- lapply(unique(fc$Meta$Library), function(lib){
			rowSums(as.matrix(fc$CPM[fc$Meta$ID[which(fc$Meta$Library == lib & fc$Meta[,ty]==sl)]]))
		})
		do.call(cbind, cpm.sum)
	})
	names(cpm.sum) <- sls
	# resolve replicates
	cts.res <- sapply(cpm.sum, function(x) apply(x, 1, function(x){
		# remove zeros?
		if(zero!="keep") x <- x[x>0]
		if(agg=="sum") {			# sum counts
			x <- sum(x)
		} else if(agg=="median"){	# median counts
			x <- median(x)
		} else {					# default to geometric mean
			x <- exp(mean(log(x)))
		}
		x[is.nan(x)] <- 0
		x
	}))
	colnames(cts.res) <- sls
	cts.res
}
cts.ty <- resolve_reps("SL.type", c("SL1", "SL2"))
cts.cl <- resolve_reps("SL.cluster", c("Cluster1", "Cluster2"))

# add read ratios (all possible ways)
cts.res <- data.frame(fc$CPM[,c(1:7)], meanTPM=rowMeans(fc$bg[-1]), cts.ty, SL2v1=cts.ty[,"SL2"]/cts.ty[,"SL1"],
			cts.cl, Cluster1v2=cts.cl[,"Cluster1"]/cts.cl[,"Cluster2"], 
			Cluster2v1=cts.cl[,"Cluster2"]/cts.cl[,"Cluster1"])
		
# read ratio plot
ddf <- rbind(data.frame(Comp="Cluster1-vs-2", Ratio=cts.res$Cluster1v2), 
			data.frame(Comp="Cluster2-vs-1", Ratio=cts.res$Cluster2v1),
			data.frame(Comp="SL2-vs-SL1", Ratio=cts.res$SL2v1))
ddf <- subset(ddf, is.finite(Ratio) & Ratio>0)
ddf$Comp <- factor(ddf$Comp, levels(ddf$Comp)[c(3,1,2)])
rr.sum <- do.call(rbind, as.list(tapply(ddf$Ratio, ddf$Comp, summary)))
cat("\nSummary of SL2:SL1 read ratio:\n")
print(rr.sum)
write.table(rr.sum, file.path(wdir, paste0(outpref, ".SL_readratio.txt")), row.names=T, col.names=NA, sep="\t", quote=F)
#rr.thresh <- 2
ddf$Comp <- factor(ddf$Comp, rev(levels(ddf$Comp)))
pdf(file.path(wdir, paste0(outpref, ".SL_readratio.pdf")), width=8, height=8)
	gg.rr <- ggplot(ddf, aes(x=Comp, y=Ratio, fill=Comp)) + 
	  annotate("rect", xmin=-Inf, xmax=Inf, ymin=rr.thresh, ymax=Inf, fill=colours.cl["SL2"], alpha=0.8) +
	  geom_violin(trim=F, bw=0.2, alpha=0.5, color="white") + 
	  #geom_jitter(shape=21, width=0.1, alpha=0.5, color="white") + 
	  geom_boxplot(width=0.1, alpha=1, color="white", outlier.shape=NA) + 
	  geom_hline(yintercept=1, color="black", linetype="dashed") +
	  scale_fill_manual(values=c("red2", "red2", "red2")) +
	  scale_y_log10(labels=label_number(accuracy = 1, big.mark = ",")) +
	  coord_flip() +
	  ylab("Read ratio") + xlab("SL contrast") +
	  theme_bw() +
	  theme(legend.position="none")
	gg.rr
invisible(dev.off())

# need to infer classes and operons for all three read-ratios
# and work out which one is best
# if SL2-types are known, ignore the clusters
# else, do it both ways and check which one gives plausible intercistronic distances.

#
# 3) infer gene classes and operons
#
# classes
inferClasses <- function(r){
	rr <- grep(paste0(r, "v"), colnames(cts.res))
	r.compl <- gsub("[12]", 2/as.integer(substr(r, nchar(r), nchar(r))), r)
	cl <- rep("SL1+SL2", times=nrow(cts.res))
	cl[cts.res[,rr]==0] <- "SL1"
	cl[is.infinite(cts.res[,rr])] <- "SL2"
	cl[is.na(cts.res[,rr])] <- "no SL"
	data.frame(cts.res[,c(1:8)], cts.res[,c(r.compl,r)], Ratio=cts.res[,rr], Class=cl, stringsAsFactors = F)
}
cts.op <- lapply(c("SL2", "Cluster1", "Cluster2"), inferClasses)

# SL CPM vs background mean TPM
da <- data.frame(TPM=rowMeans(fc$bg[-1])[match(cts.op[[1]]$Geneid, fc$bg$Geneid)], 
				SL=cts.op[[1]]$Class,
				CPM=rowSums(cts.res[,c("SL1", "SL2")]))

gg.tpm1 <- ggplot(da, aes(x=SL, y=TPM, fill=SL)) + 
	geom_violin(trim=F, bw=0.1, alpha=1, color="white") +
	geom_boxplot(color="white", width=0.2, outlier.shape=NA) + 
	scale_y_log10(labels=comma) + 
	scale_fill_manual(values=colours.cl) +
	ylab("Background (TPM)") +
	theme_classic() +
	theme(legend.position="none")

gg.tpm2 <- ggplot(subset(da, SL!="no SL"), aes(y=TPM, x=CPM, fill=SL)) +
	geom_point(shape=21, colour="white") +
	stat_smooth() +
	facet_wrap(~SL, scales="free") +
	scale_fill_manual(values=colours.cl) +
	scale_x_continuous(labels=comma) +
	scale_y_continuous(labels=comma) +
	xlab("SL (CPM)") + ylab("Background (TPM)") +
	theme_classic() +
	theme(legend.position="none")
gg.tpm <- arrangeGrob(gg.tpm1, gg.tpm2, widths=c(0.3, 0.7), ncol=2)

pdf(file.path(wdir, paste0(outpref, ".CPM_vs_TPM.pdf")), width=12, height=5)
	grid.draw(gg.tpm)
invisible(dev.off())

# operons
inferOperons <- function(cts, cutoff=Inf){
	# set cutoff to maximum intergenic distance to enable comparisons with Infinity
	if(cutoff==Inf) cutoff <- max(cts$Distance, na.rm=T)
	cts$Status <- rep("not trans-spliced", times=nrow(cts))
	cts$Quality <- "pass"
	cts$Operon <- rep(NA, times=nrow(cts))
	cts <- mclapply(split(cts, interaction(cts$Chr, cts$Strand), drop=T), function(chrs){
		# find runs of downstream genes on each chromosome and strand
		# revert if negative strand
		if(any(chrs$Strand=="-")) { chrs <- chrs[nrow(chrs):1,] }

		# all SL2-classed genes are downstream genes by default
		chrs$Status[chrs$Class=="SL2"] <- "downstream"
		  
		# SL1 and SL1+SL2 genes are monocistronic by default
		chrs$Status[chrs$Class %in% c("SL1", "SL1+SL2")] <- "monocistronic"
		  
		# designate SL1+SL2 genes with SL2:SL1 ratio above user-defined threshold as downstream
		# if set to zero, any trans-spliced gene will be designated as downstream
		chrs$Status[chrs$Ratio >= rr.thresh] <- "downstream"

		# if intergenic distance constraints have been provided, correct false positive downstream genes
		# in this case, operons have to have at least two genes (2xdown or 1up+1down)
		# classify downstream genes as monocistronic if distance is > cutoff or NA (final gene on contig) ...
		# ...AND distance of previous gene is also > cutoff (to avoid last downstream genes in operons being removed)
		d <- which(chrs$Status=="downstream" & (chrs$Distance>cutoff | is.na(chrs$Distance)))
		m <- d[c(Inf, chrs$Distance)[d]>cutoff]
		chrs$Status[m] <- "monocistronic"

		# if previous gene is <= cutoff AND not another downstream 
		# AND upstream genes must have SL2-bias, also classify as monocistronic
		d <- which(chrs$Status=="downstream" & (chrs$Distance>cutoff | is.na(chrs$Distance)))
		m <- d[c(Inf, chrs$Distance)[d]<=cutoff & c("dummy", chrs$Status)[d]!="downstream" & upstream.bias=="yes"]
		chrs$Status[m] <- "monocistronic"

		# need to break apart runs of downstream genes that are multiple operons back to back
		# this only happens when upstream genes must have SL2 bias. Otherwise, the runs are broken up by upstream genes anyway.
		d <- which(chrs$Status=="downstream" & chrs$Distance>cutoff)
		b <- d[c(chrs$Status, "dummy")[d+1]=="downstream"]
		b <- sort(c(1:nrow(chrs), b))
		chrs <- chrs[b,]
		chrs$Status[duplicated(b)] <- "dummy"
		
		# downstream genes <= cutoff that are not followed by another downstream gene suggest incomplete operons (SL missing)
		# if the previous gene is <= cutoff AND upstream genes may be unbiased, leave the gene as downstream.
		# otherwise, the gene is an orphan, so classify as monocistronic
		d <- which(chrs$Status=="downstream" & chrs$Distance<=cutoff)
		m <- d[c("dummy", chrs$Status)[d]!="downstream" & 
				c(chrs$Status, "dummy")[d+1]!="downstream" & (upstream.bias=="yes" | 
			(upstream.bias=="no" & c(Inf, chrs$Distance)[d]>cutoff))]
		chrs$Status[m] <- "monocistronic"

		# define operon annotations from runs of downstream genes
		# upstream genes can then be added on, or be part of the run if SL2-bias required
		rl <- rle(chrs$Status)
		for(rn in which(rl$values=="downstream")){
		  dstr.end <- sum(rl$lengths[0:rn])
		  dstr.start <- dstr.end-(rl$lengths[rn]-1)
		  # unbiased upstream gene is right before dstr.start (and must conform to distance constraint)
		  if(upstream.bias=="no" & dstr.start>1) { 
			# don't assign a dummy gene (from breaking multiple operons above) as upstream!!!
			if(chrs$Distance[dstr.start-1]<=cutoff & chrs$Status[dstr.start-1]!="dummy") {
				# assign previous gene as upstream
				first <- dstr.start-1
			} else {
				# need an upstream gene, but none available
				# so, designate first downstream gene as upstream
				first <- dstr.start
				# designate gene quality as "adhoc"
				chrs$Quality[chrs$Geneid==chrs$Geneid[first]] <- "adhoc"
			}
		  } else if(upstream.bias=="yes") {
		    # upstream gene requires SL2 bias
			# so, designate first downstream gene as upstream
			first <- dstr.start
		  } else {
			# need an upstream gene, but none available
			# so, designate first downstream gene as upstream
			first <- dstr.start
			# designate gene quality as "adhoc"
			chrs$Quality[chrs$Geneid==chrs$Geneid[first]] <- "adhoc"
		  }
		  chrs$Status[chrs$Geneid==chrs$Geneid[first]] <- "upstream"
		  # provisional operon ID
		  op.id <- paste0("predicted:Operon_", chrs$Chr[1], "-", chrs$Start[first],
						"-", chrs$End[dstr.end],
						 "(", chrs$Strand[1], ")")
		  chrs$Operon[chrs$Geneid %in% chrs$Geneid[first:dstr.end]] <- op.id
		  chrs$Distance[chrs$Geneid %in% chrs$Geneid[dstr.end]] <- NA 
		}  
		# revert if negative strand
		if(any(chrs$Strand=="-")) { chrs <- chrs[nrow(chrs):1,] }
		#print(chrs)
		return(chrs)
	}, mc.cores=threads)
	cts <- do.call(rbind, cts)
	cts <- subset(cts, Status != "dummy")
	cts$Status <- factor(cts$Status, levels=c("upstream", "downstream", "monocistronic", "not trans-spliced"))
	
	# replace provisional operon ID with proper ID
	cts <- cts[order(cts$Chr, cts$Start),]
	op.id <- cts$Operon[!is.na(cts$Operon)]
	op.id <- factor(op.id, levels=unique(op.id))
	op.id <- paste0(op.prefix, as.numeric(op.id))
	cts$Operon[!is.na(cts$Operon)] <- op.id
	return(cts)
}
#cts.op <- lapply(cts.op, inferOperons)
# first, infer cutoffs
#dist.cutoff="x"
if(dist.cutoff=="x") { #infer
	tmp <- lapply(cts.op, inferOperons)
	# bisect the distance distribution (K-means clustering)
	cutoffs <- sapply(tmp, function(x){
		x <- subset(x, !is.na(Operon))$Distance
		x[x<0] <- 0
		x <- na.exclude(log10(x+1))
		if(length(x)>=3 & length(unique(x))>=2) {
			km <- kmeans(x, 2)
			cu <- round(10^min(tapply(x, km$cluster, max)), digits=0)-1
		} else {
			cu = Inf
		}
		cu
	})
	cat("\nInferred intercistronic distance cutoffs:")
	print(unname(data.frame(c("SL2", "Cluster1", "Cluster2"), cutoffs)), row.names = FALSE)
} else { # supplied
	cutoffs <- rep(as.numeric(dist.cutoff), 3)
}
# then, run operon inference
cts.op <- mapply(inferOperons, cts.op, cutoffs, SIMPLIFY=FALSE)
names(cts.op) <- c("SL2", "Cluster1", "Cluster2")

	# within each operon, correct distant genes to monocistronic
	#m <- which(!is.na(cts$Operon) & cts$Distance>cutoff)
	#cts[m,"Operon"] <- NA
	#cts[m,"Status"] <- "monocistronic"
	# we also cannot have operons comprising only one gene (distance would be senseless)
	#z <- which(y>50)
	#which(c(Inf, y)[z]>50)

# summarise intergenic distances
sink(file.path(wdir, paste0(outpref, ".intergenic_distances.txt")))
for(x in c("SL2", "Cluster1", "Cluster2")){
	cat(paste0(x, "\n"))
	di <- rbind(summary(cts.op[[x]]$Distance[!is.na(cts.op[[x]]$Operon)]),
		summary(cts.op[[x]]$Distance[is.na(cts.op[[x]]$Operon)]))
	rownames(di) <- c("operonic", "not operonic")
	print(di)
	cat("\n")
}
sink()
cat("\nSummary of intergenic distances:\n\n")
cat(readLines(file.path(wdir, paste0(outpref, ".intergenic_distances.txt"))), sep="\n")
# plot
ddf <- lapply(cts.op, function(x){
	x <- x[,c("Operon", "Distance")]
	x$Operonic <- "not operonic"
	x$Operonic[!is.na(x$Operon)] <- "operonic"
	x$Operonic <- factor(x$Operonic, levels=c("operonic", "not operonic"))
	x
})
ddf <- melt(ddf)
ddf$value[ddf$value<0] <- 0
ddf$L1 <- factor(ddf$L1, rev(names(cts.op)))
pdf(file.path(wdir, paste0(outpref, ".intergenic_distances.pdf")), width=8, height=8)
	gg.d <- ggplot(ddf, aes(x=L1, y=value, fill=Operonic)) +
	  geom_violin(trim=F, bw=0.3, alpha=1, color="white", position=position_dodge(width = 0.7)) +
	  #geom_point(shape=21, alpha=0.2, color="white", position=position_dodge(width = 0.7)) + 
	  geom_boxplot(color="white", width=0.2, position=position_dodge(width = 0.7), outlier.shape=NA) +
	  scale_fill_manual(name="Gene class", values=colours.cl[levels(ddf$Operonic)]) +
	  scale_y_log10(labels=label_number(accuracy = 1, big.mark = ","), limits=c(1, 1e6)) +
	  scale_x_discrete(position="top") +
	  #scale_x_discrete(labels=c("Cluster1", "Cluster2", "SL2")) +
	  guides(fill=guide_legend(ncol=1)) +
	  xlab("Polycistron resolvers") + ylab("Intergenic distance (bp)") +
	  coord_flip () +
	  theme_bw() +
	  theme(legend.position="right")
	gg.d
invisible(dev.off())

pdf(file.path(wdir, paste0(outpref, ".combined.pdf")), width=12, height=22)
	#grid.arrange(fc$MDSplots, gg.rr, gg.d, layout_matrix=matrix(c(1,2,1,3), ncol=2), 
		#widths=c(0.4, 0.6), heights=c(0.6,0.4))
	grid.arrange(fc$MDSplots, gg.prop, gg.tpm, gg.rr, gg.d, 
					layout_matrix=matrix(c(1,2,3,4,1,2,3,5), ncol=2), heights=c(0.35,0.25,0.2,0.2)) 
		#widths=c(0.4, 0.6), heights=c(0.6,0.4))
invisible(dev.off())

# print summaries

# operon sizes
sz <- lapply(cts.op, function(x) { 
	x <- data.frame(table(table(x$Operon)))
	# if no operons exist, make a dummy dataframe instead
	if(nrow(x)==0) x <- data.frame(Var1=factor(2), Freq=NA)
	x
})
sz <- Reduce(function(x,y) merge(x = x, y = y, by="Var1", all=T), sz)
sz <- sz[order(as.numeric(as.character(sz$Var1))),]
rownames(sz) <- paste0("n=", sz[,1])
sz <- t(sz[,-1])
n.ops <- t(sapply(cts.op, function(x) { c(sum(table(x$Operon)), sum(x$Quality=="adhoc"), length(unique(na.exclude(x$Operon)))) }))
colnames(n.ops) <- c("operonic", "adhoc", "operons")
rownames(n.ops) <- c("SL2", "Cluster1", "Cluster2")

cat("Predicted operonic genes, operons and operon sizes (n = genes in operon):\n")
print(cbind(n.ops, sz))
cat("\nFinal gene classes:\n")
gcl <- do.call(cbind, lapply(cts.op, function(x) cbind(table(x$Status))))
colnames(gcl) <- names(cts.op)
print(gcl)
cat("\n")

# write operon summaries to file
write.table(cbind(t(gcl), n.ops, sz), file.path(wdir, paste0(outpref, ".operons_summary.txt")), quote=F, col.names=NA, sep="\t")

#
# 4) write GFF3
# 
write.gff3 <- function(pref){
	pref.inv <- gsub("[12]", 2/as.integer(substr(pref, nchar(pref), nchar(pref))), pref)
	dd <- cts.op[[pref]]
	ops <- subset(dd, !is.na(Operon))
	operons <- split(ops, ops$Operon)
	if(length(operons)==0) { # no operons predicted
		cat(paste("No operons predicted for", pref, "\n"))
		return(NULL)
	}
	operons <- operons[order(sapply(operons, function(x) x$Chr[1]), sapply(operons, function(x) min(x$Start)))]
	gff3 <- lapply(1:length(operons), function(i){
		opid <- unique(operons[[i]]$Operon) #opid <- paste0(op.prefix, i)     # operon name
		strand <- unique(operons[[i]]$Strand)[1]
		opg <- nrow(operons[[i]])
		op.qual <- ifelse(any(operons[[i]]$Quality=="adhoc"), "provisional", "pass")
		gff.1 <- paste(operons[[i]]$Chr[1], ".", "operon", 
				  min(operons[[i]]$Start), max(operons[[i]]$End), ".", operons[[i]]$Strand[1], ".", 
				  paste0("ID=", opid, ";Name=", opid, ";Note=", 
						"size:", nrow(operons[[i]]), 
						", genes:", paste(operons[[i]]$Geneid, collapse=","), 
						", quality:", op.qual), sep="\t")
		gff.2 <- sapply(1:opg, function(g){
		  opgene <- operons[[i]][g,]
		  gid <- paste(opid, ".", ifelse(strand=="+", g, opg+1-g), sep="")
		  paste(opgene[,"Chr"], ".", "gene", 
				opgene[,"Start"], opgene[,"End"], ".", opgene[,"Strand"], ".", 
				paste0("ID=", gid, 
					  ";Name=", opgene[,"Geneid"], 
					  ";Parent=", opid, 
					  ";Note=", "meanTPM:", round(opgene[,"meanTPM"], 2),
					  ", SL1:", round(opgene[,grep(pref.inv, colnames(dd))], 2),
					  ", SL2:", round(opgene[,grep(pref, colnames(dd))], 2),
					  ", SL2:SL1-ratio:", opgene[,"Ratio"], 
					  ", intercistronic distance:", opgene[,"Distance"],
					  ", ", opgene[,"Status"], " in operon", ", quality:", opgene[,"Quality"])
				, sep="\t")
		})
		c(gff.1, gff.2)
	})
	gff3 <- do.call(c, gff3)
	
	gff3 <- c("##gff-version 3", 
					paste("# predicted using SLOPPR", sloppr_version),
					paste("# Libraries:", paste(unique(fc$Meta$Library), collapse=", ")), 
					paste("# SL1-type SLs:", paste(unique(subset(fc$Meta, SL.type==pref.inv | SL.cluster==pref.inv)$SL), collapse=", ")),
					paste("# SL2-type SLs:", paste(unique(subset(fc$Meta, SL.type==pref | SL.cluster==pref)$SL), collapse=", ")),
					gff3)
	gff.f <- file.path(wdir, paste0(outpref, ".operons.", pref, ".gff3"))
	writeLines(gff3, gff.f)
	#write.table(gff3, gff.f, sep="\t", col.names=F, row.names=F, quote=F)
	# write resolved SL count matrix for the record
	write.table(dd, file.path(wdir, paste0(outpref, ".operons.", pref, ".txt")), row.names=F, col.names=T, sep="\t", quote=F)
	cat(paste0("Operon annotations written to ", gff.f, ".\n"))
}
for(a in names(cts.op)){ write.gff3(a) }

q()





