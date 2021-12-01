setwd('/Users/benayoun/Desktop/Work/Single_cell/2020-04-18_BMDM_Single_Cell/M3_analysis')
options(stringsAsFactors = F)
options(connectionObserver = NULL)

library(clusterProfiler)
library(org.Mm.eg.db)

# 2021-07-08
# get marker files and ratios

# 2021-12-01
# use updated list on rerun

cluster5   <- read.table('2021-12-01_cluster5_MarkersWilcoxon_FDR5.txt', header = T, sep = "\t")
cluster6   <- read.table('2021-12-01_cluster6_MarkersWilcoxon_FDR5.txt', header = T, sep = "\t")
background <- read.table('2021-12-01_All_BMDMs_scGenes_Background.txt', header = T, sep = "\t")

cluster5.pos <- rownames(cluster5)[cluster5$avg_logFC >0]   # 111
cluster6.pos <- rownames(cluster6)[cluster6$avg_logFC >0]   # 181

# keytypes(org.Mm.eg.db)
entrezID.cluster5.pos <- bitr(cluster5.pos , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.cluster6.pos <- bitr(cluster6.pos , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.background   <- bitr(background$x   , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

######################################################################################################################################################
### CLUSTER 5
######## A. GO over-representation test
ego.clust5 <- enrichGO(gene          = entrezID.cluster5.pos$ENTREZID,
                       universe      = entrezID.background$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

# write results to file
write.table(ego.clust5@result,  file = paste(Sys.Date(),"Cluster5_UP_GO_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# make some plots
pdf(paste(Sys.Date(),"Top20_Cluster5_UP_GO_dotplot_FDR5.pdf", sep = "_"), width = 12, height = 8)
dotplot(ego.clust5,  x = "Count", showCategory=20, title = "Upregulated GO Terms")
dev.off()

pdf(paste(Sys.Date(),"Cluster5_UP_GO_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 10, height = 10)
cnetplot(ego.clust5, colorEdge = TRUE, categorySize="pvalue")
dev.off()

######## B. KEGG over-representation test
kk.clust5 <- enrichKEGG(gene          = entrezID.cluster5.pos$ENTREZID,
                        universe      = entrezID.background$ENTREZID,
                        keyType       = 'ncbi-geneid',
                        organism      = 'mmu',
                        pvalueCutoff  = 0.05)

# write results to file
write.table(kk.clust5@result,  file = paste(Sys.Date(),"Cluster5_UP_KEGG_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# make some plots
pdf(paste(Sys.Date(),"Top20_Cluster5_UP_KEGG_dotplot_FDR5.pdf", sep = "_"), width = 12, height = 8)
dotplot(kk.clust5,  x = "Count", showCategory=20, title = "Upregulated KEGG Terms")
dev.off()

pdf(paste(Sys.Date(),"Cluster5_UP_KEGG_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 20, height = 20)
cnetplot(kk.clust5, colorEdge = TRUE, categorySize="pvalue")
dev.off()





###################6##################################################################################################################################
### CLUSTER 6
######## A. GO over-representation test
ego.clust6 <- enrichGO(gene          = entrezID.cluster6.pos$ENTREZID,
                       universe      = entrezID.background$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

# write results to file
write.table(ego.clust6@result,  file = paste(Sys.Date(),"cluster6_UP_GO_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")

# make some plots
pdf(paste(Sys.Date(),"Top20_cluster6_UP_GO_dotplot_FDR5.pdf", sep = "_"), width = 12, height = 8)
dotplot(ego.clust6,  x = "Count", showCategory=20, title = "Upregulated GO Terms")
dev.off()

pdf(paste(Sys.Date(),"cluster6_UP_GO_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 10, height = 10)
cnetplot(ego.clust6, colorEdge = TRUE, categorySize="pvalue")
dev.off()

######## B. KEGG over-representation test
kk.clust6 <- enrichKEGG(gene          = entrezID.cluster6.pos$ENTREZID,
                        universe      = entrezID.background$ENTREZID,
                        keyType       = 'ncbi-geneid',
                        organism      = 'mmu',
                        pvalueCutoff  = 0.05)

# write results to file
write.table(kk.clust6@result,  file = paste(Sys.Date(),"cluster6_UP_KEGG_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# make some plots
pdf(paste(Sys.Date(),"Top20_cluster6_UP_KEGG_dotplot_FDR5.pdf", sep = "_"), width = 12, height = 8)
dotplot(kk.clust6,  x = "Count", showCategory=20, title = "Upregulated KEGG Terms")
dev.off()

pdf(paste(Sys.Date(),"cluster6_UP_KEGG_GeneConcept_Network_FDR5.pdf", sep = "_"), width = 20, height = 20)
cnetplot(kk.clust6, colorEdge = TRUE, categorySize="pvalue")
dev.off()


#######################
sink(file = paste(Sys.Date(),"_ClusterProfiler_session_Info.txt", sep =""))
sessionInfo()
sink()

