# load libraries
# R version 3.4.1 (2017-06-30)
library(phenoTest) # phenoTest 1.24.0
library(qusage) # qusage 2.10.0

# set options
options(stringsAsFactors = F)

######################################################

# PART I: RETRIEVING GENE SETS
# retrieve C5 (GO) gmt file from Molecular Signature Database (MSigDB) 
# http://software.broadinstitute.org/gsea/msigdb/genesets.jsp?collection=C5

# source code 
write.genesets("c5.all.v6.2.symbols.gmt")

####################
##### FUNCTION #####
####################
# This function creates separate txt files of GeneIDs for each pathway in a gmt file
# INPUT: gmt.file = name of gmt file

write.genesets <- function(gmt.file) {
  
  # read gmt file
  pathway <- read.gmt(gmt.file)
  
  # create a table of GeneIDs for each pathway set (txt file)
  for (x in c(1:length(pathway))) {
    write.table(pathway[x], file = names(pathway[x]), row.names=FALSE, quote=FALSE, col.names=FALSE)
  }
  
}

######################################################

# PART II: CONDUCT GSEA

# create command to call on all gene set txt files within the gene set collection
go_files = list.files()[startsWith(list.files(), "GO")]

# source code
call.gsea("2019-08-21_PMA_vs_MOTScPMA_THP1_RNAseq_all_genes_statistics.txt", go_files, "GO", "THP1_PMA_vs_MOTScPMA")
call.gsea("2019-08-21_CTL_vs_PMA_THP1_RNAseq_all_genes_statistics.txt", go_files, "GO", "THP1_CTL_vs_PMA")
call.gsea("2019-08-21_CTL_vs_MOTScPMA_THP1_RNAseq_all_genes_statistics.txt", go_files, "GO", "THP1_CTL_vs_MOTScPMA")

######################
###### FUNCTION ######
######################
# my.gsea and call.gsea functions work together to conduct gene set enrichment analysis (GSEA) using 
# DESeq2 output files and the C5 (GO) gene set from MSigDB
# INPUT: RNA.data = all gene statistics file from DESeq2
#        pathway.file = txt file of GeneIDs from a pathway
#        gene.set = name of gene set collection that pathways are from MSigDB
#        sample = sample name from DEseq2 all gene statistics file
#        nperm = number of permutations
#        nmingenes = number of minimum genes
#        nmaxgenes = number of maximum genes

my.gsea <- function(RNA.data, pathway.file, gene.set, sample, nperm = 1000, nmingenes = 0, nmaxgenes = 10000000000) {
  
  # get RNAseq statistics
  my.data1 <- read.csv(RNA.data, sep = "\t", header = T)
  
  # read RNAseq statistics and get GeneID and logfold2change
  my.data <- data.frame("GeneID" = as.character(rownames(my.data1)),
                        "log2FoldChange" = my.data1$log2FoldChange)
  my.data$GeneID <- sapply(my.data$GeneID, tolower)
  
  # read pathway list (contains GeneIDs)
  pathway <- read.csv(pathway.file, sep = "\t", header = T)
  pathway1 <- sapply(pathway[1], tolower)
  
  # create signature (a matrix that contains the GeneIDs of pathway and RNAseq statistics)
  my.signature = list(pathway1)
  names(my.signature) = pathway.file
  
  # reformat my.data for GSEA
  my.matrix = my.data$log2FoldChange
  names(my.matrix) = my.data$GeneID
  
  # run GSEA and get enrichment score
  tryCatch({
    print(pathway.file)
    gsea.data <- gsea(x = my.matrix, gsets = my.signature, mc.cores=2, logScale = F, B = nperm, minGenes = nmingenes, maxGenes = nmaxgenes)
    my.summary <- summary(gsea.data)
    return(my.summary)
  }, error = function(cond) {return("NoMatch")}) # for genes in the RNAseq data that are not listed as a gene in the gene set, return with "NoMatch"
}

call.gsea <- function(RNA.data, my.pathways, gene.set, sample) {
  my.current.pathway <- data.frame(matrix(ncol = 6, nrow = 0))
  colnames(my.current.pathway) <- c("n",	"es", "nes",	"pval.es",	"pval.nes",	"fdr")
  
  for (filename in my.pathways) {
    gsea.scores <- my.gsea(RNA.data, filename, gene.set, sample)
    if (gsea.scores == "NoMatch" || anyNA(gsea.scores)) {next;}
    my.current.pathway <- rbind(my.current.pathway,gsea.scores)
  }
  
  my.prefix <- paste(Sys.Date(),"GSEA",sample, gene.set,sep = "_")
  
  # run FDR correction
  my.current.pathway$fdr <- p.adjust(my.current.pathway$pval.nes, method = "fdr")
  
  write.table(my.current.pathway, append = F, file = paste(my.prefix,"all_pathways.txt",sep = "_"), sep = "\t" , row.names = T, quote=F)   
  write.table(my.current.pathway[my.current.pathway$fdr < 0.05,], append = F, file = paste(my.prefix,"FDR5_pathways.txt",sep = "_"), sep = "\t" , row.names = T, quote=F)   
}
