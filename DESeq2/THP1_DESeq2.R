# load libraries
# R version 3.4.1 (2017-06-30)
library(DESeq2) # DESeq2 1.16.1
library(pvclust) # pvclust 2.0.0

# set options
options(stringsAsFactors = F)

######################################################

# source code
analyze.tissue("2019-05-13_THP1_kallisto_mappings.txt", "ALL", 0.05)
analyze.tissue("2019-05-13_THP1_kallisto_mappings.txt", "PMA_vs_MOTScPMA", 0.05)
analyze.tissue("2019-05-13_THP1_kallisto_mappings.txt", "CTL_vs_PMA", 0.05)

######################
###### FUNCTION ######
######################
# This function takes in the THP1 kallisto-mapped gene counts file and processes it through DESeq2 modeling to find 
# differentially expressed genes between control, PMA treated, and MOTSc + PMA treated samples. Using DESeq2 results, the function will
# also create a volcano plot (for PMA_vs_MOTScPMA) and multidimensional scaling (MDS) analysis plot (for all).
# INPUT: counts.file = 2019-04-25_C2C12_RNAseq_kallisto_mapping.txt
#        comparison = "ALL", "PMA_vs_MOTScPMA", "CTL_vs_PMA", "CTL_vs_MOTScPMA"
#        FDR = false discovery rate (0.05)

analyze.tissue <- function(counts.file, comparison, FDR) {

  # read kallisto mappings
  my.data <- read.csv(counts.file, sep = "\t", header = T)
  
  # reorder columns (ENS_GID and tGeneSymbol are switched) and rename with correct colnames
  my.data <- my.data[c(2,1,3:19)]
  colnames(my.data)[1:2] <- c("tGeneSymbol", "ENS_GID") 
  
  # sum read over genes (to not have results over transcripts for DEseq2)
  my.data.per.gene <- aggregate(my.data[,5:19],by=list(my.data$tGeneSymbol),FUN=sum)

  # round counts (DESeq needs integers)
  my.data.per.gene[,2:16] <- round(my.data.per.gene[,2:16])
  rownames(my.data.per.gene) <- my.data.per.gene$Group.1
  
  # get the genes with no reads out
  my.null <- which(apply(my.data.per.gene[,2:16], 1, sum) <= 5) # see deseq2 vignetter
  my.filtered.matrix <- my.data.per.gene[-my.null,2:16]
  
  # get number of mapped reads
  print(apply(my.filtered.matrix,2,sum))
  
  # remove genes that start with "ERCC"
  remove <- grep("ERCC", row.names(my.filtered.matrix), ignore.case = F, fixed = T)
  my.filtered.matrix <- my.filtered.matrix[-remove,]
  
  # add prefix to colnames to reflect samples
  colnames(my.filtered.matrix) <- paste(c(rep("CTL",5), rep("PMA", 5), rep ("MOTScPMA",5)), colnames(my.filtered.matrix), sep = "_")
  
  if (comparison == "ALL") {
    # prepare matrix for DESeq2
    my.status <- rep("ALL",dim(my.filtered.matrix)[2])
    my.status[grep(("CTL"),colnames(my.filtered.matrix))] <- "CTL"
    my.status[grep(("PMA"),colnames(my.filtered.matrix))] <- "PMA"
    my.status[grep(("MOTScPMA"),colnames(my.filtered.matrix))] <- "MOTScPMA"
    
    dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                             status = my.status )
    
    # get matrix using treatment as a modeling covariate
    dds <- DESeqDataSetFromMatrix(countData = my.filtered.matrix,
                                  colData = dataDesign,
                                  design = ~ status)
    
    # run DESeq normalizations and export results
    dds.deseq <- DESeq(dds) # no outliers reported
    
    # parse sample names
    my.sample.names <- unlist(strsplit(colnames( my.filtered.matrix ), c("_abundance.tsv")))
    
    # assign colors
    my.colors <- c(rep("orange",5), rep("dodgerblue",5), rep("dodgerblue4",5))
    
    # determine normalized expression value
    tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
    
    colnames(tissue.cts) <- c(paste("CTL_S",c(1:5),sep=""), paste("PMA_S",c(6:10),sep=""), paste("MOTScPMA_S",c(11:15),sep=""))
    
    ### MDS Analysis Plot
    mds.result <- cmdscale(1-cor(my.filtered.matrix,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
    x <- mds.result[, 1]
    y <- mds.result[, 2]
    
    # make plot with sample names labelled
    my.mds.out <- paste(Sys.Date(), comparison, "THP1_analysis_MDS_plot_with_sample_names.pdf", sep ="_")
    
    pdf(my.mds.out)
    plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling")
    points(x, y, pch=16,col=my.colors,cex=2)
    text(x, y,my.sample.names ,col="black",cex=0.5, pos  = 1)
    legend("topleft",c("CTL","PMA", "MOTScPMA"),col=c("orange","dodgerblue", "dodgerblue4"),pch=16,bty='n',pt.cex=2)
    dev.off()
    
    # make plot with no sample names labelled
    my.mds.out <- paste(Sys.Date(), comparison, "THP1_analysis_MDS_plot_without_sample_names.pdf", sep ="_")
    
    pdf(my.mds.out)
    plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",main="Multi-dimensional Scaling")
    points(x, y, pch=16,col=my.colors,cex=2)
    legend("topleft",c("CTL","PMA", "MOTScPMA"),col=c("orange","dodgerblue", "dodgerblue4"),pch=16,bty='n',pt.cex=2)
    dev.off()
    
    ## Clustering Plot
    my.pv <- pvclust(tissue.cts,nboot=100)
    my.heatmap.out <- paste(Sys.Date(), comparison, "THP1_PVCLUST_result.pdf", sep = "_")
    
    pdf(my.heatmap.out)
    plot(my.pv)
    dev.off()
    
    # output result tables to files
    write.table(tissue.cts, file = paste(Sys.Date(),"ALL_TISSUES_RNAseq_log2_counts.txt", sep = "_") , sep = "\t" , row.names = T, quote=F)
  }
  
  if (comparison == "PMA_vs_MOTScPMA") {
    # prepare matrix for DESeq2
    a.subset <- subset.data.frame(my.filtered.matrix, subset =TRUE, select =6:15, drop = FALSE)
    my.a <- rep("PMA_vs_MOTScPMA",dim(a.subset)[2])
    my.a[grep(("PMA"),colnames(a.subset))] <- "PMA"
    my.a[grep(("MOTSc"),colnames(a.subset))] <- "MOTScPMA"
    
    dataDesign = data.frame( row.names = colnames( a.subset ), 
                             PMA_vs_MOTScPMA = my.a )
    
    # get matrix using treatment as a modeling covariate
    dds <- DESeqDataSetFromMatrix(countData = a.subset,
                                  colData = dataDesign,
                                  design = ~ PMA_vs_MOTScPMA)
    
    # run DESeq normalizations and export results
    dds.deseq <- DESeq(dds) # no outliers reported
    
    res <- results(dds.deseq, contrast=c("PMA_vs_MOTScPMA","MOTScPMA", "PMA")) # added the name of the tested variable: doesn't seem to be taken correctly by default for FC
    
    # exclude NA in res
    res <- res[!is.na(res$padj),]
    
    # find number of significant genes from res at FDR
    my.genes <- rownames(res)[res$padj < FDR]
    my.num <- length(my.genes)
    print(my.num) # 945 significant genes
    
    # parse sample names
    my.sample.names <- unlist(strsplit(colnames( a.subset ), c("_abundance.tsv")))
    
    # assign colors
    my.colors <- c(rep("dodgerblue",5), rep("dodgerblue4",5))
    
    # determine normalized expression value
    tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
    
    colnames(tissue.cts) <- c(paste("PMA_S",c(6:10),sep=""), paste("MOTScPMA_S",c(11:15),sep=""))
    
    ### Volcano Plot
    # with 10 most significant genes annotated
    my.sig <- res$padj < 0.05
    my.volcano.out <- paste(Sys.Date(), comparison, "THP1_Volcano_plot.pdf", sep = "_")

    my.dwn.x <- subset(res$log2FoldChange[my.sig], res$log2FoldChange[my.sig] < 0)
    my.dwn.y <- subset(res$padj[my.sig], res$log2FoldChange[my.sig] < 0)

    my.up.x <- subset(res$log2FoldChange[my.sig], res$log2FoldChange[my.sig] > 0)
    my.up.y <- subset(res$padj[my.sig], res$log2FoldChange[my.sig] > 0)

    res1 <- res[with(res, order(res$padj)),]
    res2 <- res1[c(1:10),]
    
    pdf(my.volcano.out)
    smoothScatter(res$log2FoldChange,-log10(res$padj), col = "black", xlim=c(-1.5,1.5))
    points(my.dwn.x ,-log10(my.dwn.y), cex= 0.6, col = "red")
    points(my.up.x ,-log10(my.up.y), cex= 0.6, col = "orange")
    text(res2$log2FoldChange, -log10(res2$padj), rownames(res2), col="black",cex=0.5, pos  = 1)
    dev.off()
    
    # output result tables to files
    my.outprefix <- paste(Sys.Date(), comparison, "THP1_RNAseq", sep = "_")
    my.out.ct.mat <- paste(my.outprefix,"log2_counts_matrix.txt", sep = "_")
    
    my.out.stats <- paste(my.outprefix,"all_genes_statistics.txt", sep = "_")
    my.out.fdr <- paste(my.outprefix,"FDR", (FDR*100), "genes_statistics.txt", sep = "_")
    my.out.rdata <- paste(my.outprefix,"statistics.RData", sep = "_")
    
    write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
    write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
    write.table(res[my.genes,], file = my.out.fdr, sep = "\t" , row.names = T, quote=F)
    save(res,file=my.out.rdata)
    
  }
    
  if (comparison == "CTL_vs_PMA") {
    # prepare matrix for DESeq2
    b.subset <- subset.data.frame(my.filtered.matrix, subset =TRUE, select =1:10, drop = FALSE)
    my.b <- rep("CTL_vs_PMA",dim(b.subset)[2])
    my.b[grep(("CTL"),colnames(b.subset))] <- "CTL"
    my.b[grep(("PMA"),colnames(b.subset))] <- "PMA"
    
    dataDesign = data.frame( row.names = colnames( b.subset ), 
                             CTL_vs_PMA = my.b )
    
    # get matrix using treatment as a modeling covariate
    dds <- DESeqDataSetFromMatrix(countData = b.subset,
                                  colData = dataDesign,
                                  design = ~ CTL_vs_PMA)
    
    # run DESeq normalizations and export results
    dds.deseq <- DESeq(dds) # no outliers reported
    
    res <- results(dds.deseq, contrast=c("CTL_vs_PMA","PMA", "CTL")) # added the name of the tested variable: doesn't seem to be taken correctly by default for FC
    
    # exclude NA in res
    res <- res[!is.na(res$padj),]
    
    # find number of significant genes from res at FDR
    my.genes <- rownames(res)[res$padj < FDR]
    my.num <- length(my.genes)
    print(my.num) # 6752 significant genes
    
    # parse sample names
    my.sample.names <- unlist(strsplit(colnames( b.subset ), c("_abundance.tsv")))
    
    # assign colors
    my.colors <- c(rep("dodgerblue",5), rep("dodgerblue4",5))
    
    # determine normalized expression value
    tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
    
    colnames(tissue.cts) <- c(paste("CTL_S",c(1:5),sep=""), paste("PMA_S",c(6:10),sep=""))
    
    # output result tables to files
    my.outprefix <- paste(Sys.Date(), comparison, "THP1_RNAseq", sep = "_")
    my.out.ct.mat <- paste(my.outprefix,"log2_counts_matrix.txt", sep = "_")
    
    my.out.stats <- paste(my.outprefix,"all_genes_statistics.txt", sep = "_")
    my.out.fdr <- paste(my.outprefix,"FDR", (FDR*100), "genes_statistics.txt", sep = "_")
    my.out.rdata <- paste(my.outprefix,"statistics.RData", sep = "_")
    
    write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
    write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
    write.table(res[my.genes,], file = my.out.fdr, sep = "\t" , row.names = T, quote=F)
    save(res,file=my.out.rdata)
  }
  
  if (comparison == "CTL_vs_MOTScPMA") {
    # prepare matrix for DESeq2
    c.subset <- subset.data.frame(my.filtered.matrix, subset =TRUE, select = c(1:5, 11:15), drop = FALSE)
    my.c <- rep("CTL_vs_MOTScPMA",dim(c.subset)[2])
    my.c[grep(("CTL"),colnames(c.subset))] <- "CTL"
    my.c[grep(("MOTScPMA"),colnames(c.subset))] <- "MOTScPMA"
    
    dataDesign = data.frame( row.names = colnames( c.subset ), 
                             CTL_vs_MOTScPMA = my.c )
    
    # get matrix
    dds <- DESeqDataSetFromMatrix(countData = c.subset,
                                  colData = dataDesign,
                                  design = ~ CTL_vs_MOTScPMA)
    
    # run DESeq normalizations and export results
    dds.deseq <- DESeq(dds) # no outliers reported
    
    res <- results(dds.deseq, contrast=c("CTL_vs_MOTScPMA","MOTScPMA", "CTL")) # added the name of the tested variable: doesn't seem to be taken correctly by default for FC
    
    # exclude NA in res
    res <- res[!is.na(res$padj),]
    
    # find number of significant genes from res at FDR
    my.genes <- rownames(res)[res$padj < FDR]
    my.num <- length(my.genes)
    print(my.num) # 8148 significant genes
    
    # parse sample names
    my.sample.names <- unlist(strsplit(colnames( c.subset ), c("_abundance.tsv")))
    
    # colors for CTL_vs_MOTScPMA
    my.colors <- c(rep("orange",5), rep("dodgerblue4",5))
    
    # determine normalized expression value
    tissue.cts <- log2( counts(dds.deseq, normalize = TRUE) + 0.01)
    
    colnames(tissue.cts) <- c(paste("CTL_S",c(1:5),sep=""), paste("MOTScPMA_S",c(11:15),sep=""))
  
    # output result tables to files
    my.outprefix <- paste(Sys.Date(), comparison, "THP1_RNAseq", sep = "_")
    my.out.ct.mat <- paste(my.outprefix,"log2_counts_matrix.txt", sep = "_")
    
    my.out.stats <- paste(my.outprefix,"all_genes_statistics.txt", sep = "_")
    my.out.fdr <- paste(my.outprefix,"FDR", (FDR*100), "genes_statistics.txt", sep = "_")
    my.out.rdata <- paste(my.outprefix,"statistics.RData", sep = "_")
    
    write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)
    write.table(res, file = my.out.stats , sep = "\t" , row.names = T, quote=F)
    write.table(res[my.genes,], file = my.out.fdr, sep = "\t" , row.names = T, quote=F)
    save(res,file=my.out.rdata)
  }
}
