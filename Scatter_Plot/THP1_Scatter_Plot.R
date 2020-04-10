# load libraries
library(ggplot2) # ggplot2 3.1.0
library(bitops) # bitops 1.0.6
library(data.table ) # data.table 1.11.8 

# set options
options(stringsAsFactors = F)

######################################################

# set working directory
setwd('/Volumes/RLAI/David_THP1/Paper_Submission/Scatter_Plot')

# read data
my.res.CTLvsPMA.all   <- read.table(list.files(".", pattern = "*CTL_vs_PMA_THP1_RNAseq_all_genes_statistics.txt"), header = T, sep ="\t")
my.res.PMA_vs_MOTScPMA.all   <- read.table(list.files(".", pattern = "*PMA_vs_MOTScPMA_THP1_RNAseq_all_genes_statistics"), header = T, sep  ="\t", row.names = 1)

my.res.CTLvsPMA.all$names <- rownames(my.res.CTLvsPMA.all)
my.res.PMA_vs_MOTScPMA.all$names <- rownames(my.res.PMA_vs_MOTScPMA.all)

# merge data into one dataframe
my.merged <- merge(my.res.CTLvsPMA.all, my.res.PMA_vs_MOTScPMA.all, by = "names", suffixes = c(".CTLvsPMA.all",".PMA_vs_MOTScPMA.all"))
rownames(my.merged) <- my.merged$names

# correlation test
my.test <- cor.test(my.merged$log2FoldChange.CTLvsPMA.all, my.merged$log2FoldChange.PMA_vs_MOTScPMA.all, method = "spearman")

# run linear regression
my.motscfit <- lm(log2FoldChange.PMA_vs_MOTScPMA.all ~ log2FoldChange.CTLvsPMA.all, data= my.merged)

# determine the ideal fit, and the standard deviation around it to identify outliers
my.perfect.ys <- predict(my.motscfit, my.merged)
my.real.y.sd <- sd(my.merged$log2FoldChange.PMA_vs_MOTScPMA.all)

my.perfect.ys.up  <- my.perfect.ys + my.real.y.sd
my.perfect.ys.dwn <- my.perfect.ys - my.real.y.sd

# extract genes significant between PMA vs MOTscPMA, not significant with just CTl vs PMA
# and that are more than one standard deviation from the fit
my.motsc.only.fdr1.filter <- rownames(my.merged)[bitAnd(bitAnd(my.merged$padj.PMA_vs_MOTScPMA.all < 0.01, my.merged$padj.CTLvsPMA.all > 0.1) > 0, 
                                                        !between(my.merged$log2FoldChange.PMA_vs_MOTScPMA.all,my.perfect.ys.dwn, my.perfect.ys.up) ) >0]

my.motsc.only.fdr5.filter <- rownames(my.merged)[bitAnd(bitAnd(my.merged$padj.PMA_vs_MOTScPMA.all < 0.05, my.merged$padj.CTLvsPMA.all > 0.5) > 0, 
                                                        !between(my.merged$log2FoldChange.PMA_vs_MOTScPMA.all,my.perfect.ys.dwn, my.perfect.ys.up) ) >0]

# extract top 20 genes with biggest fold gene satisfying the above criteria
my.sort <- sort(abs(my.merged[my.motsc.only.fdr5.filter,]$log2FoldChange.PMA_vs_MOTScPMA.all), index.return = T, decreasing = T)
my.highlight <- rownames(my.merged[my.motsc.only.fdr5.filter,])[my.sort$ix][1:20]

### Scatter Plot
pdf(paste(Sys.Date(),"_Plot_THP1_CTLvsPMA_vs_PMAvsMOTScPMA_all_DEseq2_annotated_specific_genes.pdf"))
smoothScatter(my.merged$log2FoldChange.CTLvsPMA.all, my.merged$log2FoldChange.PMA_vs_MOTScPMA.all, ylim = c(-1.5,1.5), xlim = c(-3,3))
abline(h = 0, col = "brown", lty = "dashed")
abline(v = 0, col = "brown", lty = "dashed")

# plot regression line
abline(my.motscfit, col = "red", lty = "dashed") 

# plot SD around regression line
points(my.merged$log2FoldChange.CTLvsPMA.all, my.perfect.ys.up, type = 'l', col = "grey") 
points(my.merged$log2FoldChange.CTLvsPMA.all, my.perfect.ys.dwn, type = 'l', col = "grey") 

# plot significant genes
points(my.merged[my.motsc.only.fdr5.filter,]$log2FoldChange.CTLvsPMA.all,         
       my.merged[my.motsc.only.fdr5.filter,]$log2FoldChange.PMA_vs_MOTScPMA.all,
       col = "deeppink", pch = 16, cex = 0.5)

# plot top 20 significant gene names
text(my.merged[my.highlight,]$log2FoldChange.CTLvsPMA.all,                 
       my.merged[my.highlight,]$log2FoldChange.PMA_vs_MOTScPMA.all,
     rownames(my.merged[my.highlight,]),
       col = "black", cex = 0.5, pos = 4)

# plot Rho and p-value
text(-2.75,1.5,paste0("Rho = ", signif(my.test$estimate,3)), pos = 4)
text(-2.75,1.3,paste0("p = ", signif(my.test$p.value,3)), pos = 4)

dev.off()

# export table of significant genes to file
write.table(my.merged[my.motsc.only.fdr5.filter,], file = paste(Sys.Date(),"_Specific_MOTSc_Reg_THP1_CTLvsPMA_vs_PMAvsMOTScPMA_all_DEseq2.txt"), sep = "\t")

