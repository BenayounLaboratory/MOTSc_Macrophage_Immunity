setwd('/Users/benayoun/Desktop/Work/Single_cell/2020-04-18_BMDM_Single_Cell/M3_analysis')
options(stringsAsFactors = F)

library('Seurat')
library(bitops)
library(MAST)
library(DoubletFinder)
library(sctransform)
library(SingleR)
library(pheatmap)
library(clustree)
library(scales)
library(EnhancedVolcano)
library(dplyr)
library(Scillus) # https://scillus.netlify.app/vignettes/plotting.html; advanced plotting

# 2021-06-24
# analyze BMDM single cell, with/without MOTSc treatment for first 3 days of differentation

# 2021-12-01
# small bug led to not discarding doublets (found by Michelle) - rerun after doublet removal

#####################################################################################################################
#### 1. Load Phantom purge object
load("../Phantom_Purge_Preprocessing/2020-06-30_PhantomPurged_BMDM_CTL_M3_Seurat.RData")

my.bmdm.all
# An object of class Seurat
# 31053 features across 15675 samples (instead of 7775 for CTL) within 1 assay 
# Active assay: RNA (31053 features, 0 variable features)

my.good.ix <-  apply(my.bmdm.all@assays$RNA@counts>0,1,sum) > 50 # expressed in at least 50 cells 
summary(my.good.ix)
#    Mode   FALSE    TRUE 
# logical   19431   11622

# remove low/null genes
subsample <- subset(my.bmdm.all, features = rownames(my.bmdm.all@assays$RNA@counts)[my.good.ix])
my.bmdm <- subsample

my.bmdm
# An object of class Seurat
# 11622 featuresacross 15675 samples within 1 assay 
# Active assay: RNA (11622 features, 0 variable features)

################################################################################################################################################################
#### 2. QC on mitochondrial reads
# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
my.bmdm[["percent.mito"]] <- PercentageFeatureSet(my.bmdm, pattern = "^mt-")
head(my.bmdm@meta.data)
tail(my.bmdm@meta.data)

pdf(paste(Sys.Date(),"BMDM_aging_violinPlots_QC_gene_UMI_mito.pdf", sep = "_"), height = 5, width = 10)
VlnPlot(object = my.bmdm, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(my.bmdm, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(my.bmdm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste(Sys.Date(),"BMDM_aging_QC_scatter.pdf", sep = "_"), height = 5, width = 10)
CombinePlots(plots = list(plot1, plot2))
dev.off()

# filter dead/low Q cells
my.bmdm <- subset(my.bmdm, subset = nFeature_RNA > 1000 & percent.mito < 10)
my.bmdm

# An object of class Seurat
# 11622 features across 15415 samples within 1 assay 
# Active assay: RNA (11622 features, 0 variable features)


################################################################################################################################################################
#### 3. Store info on biological origin of sample

head(my.bmdm@meta.data) # data shown below from CTL. TODO update head/tail

######                          barcode   sum detected Condition Treatment nCount_RNA nFeature_RNA percent.mito
######4m_F_CTL_cell_1 AAACCTGCATTATCTC_5 10752     2800  4m_F_CTL       CTL      10747         2796     2.586768
######4m_F_CTL_cell_2 AAACCTGGTACGCTGC_5  6250     1744  4m_F_CTL       CTL       6247         1741     2.337122
######4m_F_CTL_cell_3 AAACCTGGTGCTAGCC_5 10802     2412  4m_F_CTL       CTL      10797         2407     2.426600
######4m_F_CTL_cell_4 AAACCTGTCACCTTAT_5  9883     2530  4m_F_CTL       CTL       9878         2525     3.927921
######4m_F_CTL_cell_5 AAACCTGTCTCGGACG_5  8528     2094  4m_F_CTL       CTL       8526         2092     3.002580
######4m_F_CTL_cell_6 AAACGGGCAGGGAGAG_5  7039     2093  4m_F_CTL       CTL       7035         2089     3.198294

################################################################################################################################################################
#### 4. Normalizing the data
# global-scaling normalization method ???LogNormalize??? that normalizes the gene expression measurements for each cell by the total expression,
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
my.bmdm <- NormalizeData(object = my.bmdm, normalization.method = "LogNormalize",  scale.factor = 10000)

################################################################################################################################################################
#### 5. Cell cycle regression
# Read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "/Users/benayoun/Desktop/Work/Single_cell/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")

# make into mouse gene names
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

cc.genes.mouse <- firstup(tolower(cc.genes))


# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes   <- cc.genes.mouse[1:43]
g2m.genes <- cc.genes.mouse[44:97]

# Assign Cell-Cycle Scores
my.bmdm <- CellCycleScoring(object = my.bmdm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = my.bmdm@meta.data)

# write predictions to file
write.table(my.bmdm@meta.data, file = paste(Sys.Date(),"BMDM_aging_MOTSc_CellCycle_predictions.txt", sep = "_"), sep = "\t", quote = F)

pdf(paste0(Sys.Date(), "BMDM_aging_MOTSc__cell_cycle_pie_charts.pdf"), height = 13, width = 8)
par(mfrow=c(4,2))
par(oma=c(0.1,0.1,0.1,0.1))
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "4m_F_CTL"]) , main = "4m_F_CTL")
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "4m_M_CTL"]) , main = "4m_M_CTL")
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "20m_F_CTL"]), main = "20m_F_CTL")
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "20m_M_CTL"]), main = "20m_M_CTL")
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "4m_F_M3"]) , main = "4m_F_M3")
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "4m_M_M3"]) , main = "4m_M_M3")
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "20m_F_M3"]), main = "20m_F_M3")
pie(table(my.bmdm@meta.data$Phase[my.bmdm@meta.data$Condition == "20m_M_M3"]), main = "20m_M_M3")
par(mfrow=c(1,1))
dev.off()

################################################################################################################################################################
##### 6.Find and remove doublets using doublet finder workflow
# https://github.com/chris-mcginnis-ucsf/DoubletFinder

my.bmdm <- SCTransform(object = my.bmdm, vars.to.regress = c("nFeature_RNA", "percent.mito", "Phase"))
my.bmdm <- RunPCA(my.bmdm)

my.bmdm <- RunUMAP(my.bmdm, dims = 1:30)

my.bmdm <- FindNeighbors(my.bmdm, dims = 1:30)
my.bmdm <- FindClusters(object = my.bmdm)

## pK Identification (no ground-truth)
sweep.res.list_bmdm <- paramSweep_v3(my.bmdm, PCs = 1:30, sct = TRUE) # for Mac computers may include: , num.cores	 = 4
sweep.stats_bmdm <- summarizeSweep(sweep.res.list_bmdm , GT = FALSE)
bcmvn_bmdm <- find.pK(sweep.stats_bmdm)

pk.bmdm <- as.numeric(bcmvn_bmdm$pK[bcmvn_bmdm$BCmetric == max(bcmvn_bmdm$BCmetric)])
# 22
bcmvn_bmdm[pk.bmdm,]
#0.2

## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(my.bmdm@meta.data$seurat_clusters)     ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi       <- round(0.023*length(my.bmdm@meta.data$barcode))        ## Assuming 2.3% doublet formation rate, based on expectation for 3000 cells in 10X genomics manual
nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
my.bmdm <- doubletFinder_v3(my.bmdm, PCs = 1:15, pN = 0.25, pK = 0.2, nExp = nExp_poi,     reuse.pANN = FALSE, sct = T)
my.bmdm <- doubletFinder_v3(my.bmdm, PCs = 1:15, pN = 0.25, pK = 0.2, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.2_355", sct = T)

pdf(paste(Sys.Date(),"BMDM_aging_MOTSc_Doublet_Finder_UMAP.pdf", sep = "_"), height = 5, width = 10)
DimPlot(my.bmdm, reduction = "umap", group.by = "DF.classifications_0.25_0.2_355")
# DimPlot(my.bmdm, reduction = "umap", group.by = "DF.classifications_0.25_0.2_320")
dev.off()

pdf(paste(Sys.Date(),"BMDM_aging_MOTSc_Adgre1_UMAP.pdf", sep = "_"), height = 5, width = 10)
FeaturePlot(my.bmdm, features = c("Adgre1"))
dev.off()

# save data for singlets
my.bmdm.singlets    <- subset(my.bmdm, subset = DF.classifications_0.25_0.2_355 %in% "Singlet")  # only keep singlets

my.bmdm.singlets
# An object of class Seurat
# 23244 features across 15060 samples within 2 assays
# Active assay: SCT 11622 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

pdf(paste(Sys.Date(),"BMDM_aging_MOTSc_Singlets_UMAP_by_Group.pdf", sep = "_"), height = 5, width = 6)
DimPlot(my.bmdm.singlets, reduction = "umap", group.by = "Condition")
dev.off()


################################################################################################################################################################
##### 7.Scaling the data and removing unwanted sources of variation
# In Seurat v2 we also use the ScaleData function to remove unwanted sources of variation from a single-cell dataset.
# For example, we could ‘regress out’ heterogeneity associated with (for example) cell cycle stage, or mitochondrial contamination.
# hese features are still supported in ScaleData in Seurat v3, i.e.:
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
# However, particularly for advanced users who would like to use this functionality, we strongly recommend the use of our new normalization workflow,
# sctransform.
# The method is described in our recent preprint, with a separate vignette using Seurat v3 here. As with ScaleData, the function SCTransform also includes
# a vars.to.regress parameter.

my.bmdm.singlets <- SCTransform(object = my.bmdm.singlets, vars.to.regress = c("nFeature_RNA", "percent.mito", "S.Score", "G2M.Score"))


################################################################################################################################################################
##### 8. Cluster the cells

# Determine the ‘dimensionality’ of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
ElbowPlot(my.bmdm.singlets, ndims = 50)
# pick 30, as it is where the slope goes down

# created an issue because not on singlets before
my.bmdm.singlets <- FindNeighbors(my.bmdm.singlets, dims = 1:30)

### test different resolutions for Seurat to find optimal
my.bmdm.singlets <- FindClusters(object = my.bmdm.singlets,
                                 reduction.type = "pca",
                                 k.param = 100,
                                 dims.use = 1:30,
                                 resolution = 0.5,
                                 print.output = 1,
                                 save.SNN = TRUE,
                                 plot.SNN = TRUE,
                                 n.start = 100)
# Number of communities: 13

my.bmdm.singlets <- FindClusters(object = my.bmdm.singlets,
                                 reduction.type = "pca",
                                 k.param = 100,
                                 dims.use = 1:26,
                                 resolution = 0.4,
                                 print.output = 1,
                                 save.SNN = TRUE,
                                 plot.SNN = TRUE,
                                 n.start = 100)
# Number of communities: 11


my.bmdm.singlets <- FindClusters(object = my.bmdm.singlets,
                                 reduction.type = "pca",
                                 k.param = 100,
                                 dims.use = 1:26,
                                 resolution = 0.3,
                                 print.output = 1,
                                 save.SNN = TRUE,
                                 plot.SNN = TRUE,
                                 n.start = 100)
# Number of communities: 10

my.bmdm.singlets <- FindClusters(object = my.bmdm.singlets,
                                 reduction.type = "pca",
                                 k.param = 100,
                                 dims.use = 1:26,
                                 resolution = 0.2,
                                 print.output = 1,
                                 save.SNN = TRUE,
                                 plot.SNN = TRUE,
                                 n.start = 100)
# Number of communities: 8


my.bmdm.singlets <- FindClusters(object = my.bmdm.singlets,
                                 reduction.type = "pca",
                                 k.param = 100,
                                 dims.use = 1:26,
                                 resolution = 0.1,
                                 print.output = 1,
                                 save.SNN = TRUE,
                                 plot.SNN = TRUE,
                                 n.start = 100)
# Number of communities: 5


my.bmdm.singlets <- FindClusters(object = my.bmdm.singlets,
                                 reduction.type = "pca",
                                 k.param = 100,
                                 dims.use = 1:26,
                                 resolution = 0,
                                 print.output = 1,
                                 save.SNN = TRUE,
                                 plot.SNN = TRUE,
                                 n.start = 100)
# Number of communities: 1


# clustree: Visualise Clusterings at Different Resolutions
# https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
pdf(paste(Sys.Date(),"BMDM_aging_MOTSc_Singlets_clustree_res.pdf", sep = "_"), height = 5, width = 6)
clustree(my.bmdm.singlets, prefix = "SCT_snn_res.")
dev.off()


# 
pdf(paste(Sys.Date(),"BMDM_aging_UMAP_clusters_0.2_res_Singlets.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(my.bmdm.singlets, label = TRUE, group.by = "SCT_snn_res.0.2")
dev.off()

pdf(paste(Sys.Date(),"BMDM_aging_UMAP_clusters_0.2_res_Singlets_SPLIT_by_GROUP.pdf", sep = "_"), height = 5, width = 25)
DimPlot(my.bmdm.singlets, label = TRUE, group.by = "SCT_snn_res.0.2", split.by = "Condition")
dev.off()


write.table(my.bmdm.singlets@meta.data, file = paste(Sys.Date(),"Seurat_PCA_SNN_clusters_and_Metadata_Singlets.txt", sep = "_"), sep = "\t", quote = F, col.names = F, row.names = T)


# 
my.freq.table <- prop.table(x = table(my.bmdm.singlets@meta.data$SCT_snn_res.0.2, my.bmdm.singlets@meta.data$Condition), margin = 2)
my.freq.table.av <- apply(my.freq.table,1,mean)
my.freq.table.av.sort <- sort(my.freq.table.av, decreasing = T,index.return = T)

#  Check proportions
pdf(paste(Sys.Date(),"BMDM_Aging_10X_clusters_res_0.2_barplot2.pdf",sep = "_"), height = 10, width = 15)
barplot(my.freq.table[my.freq.table.av.sort$ix,], las = 2,
        legend.text = rownames(my.freq.table)[my.freq.table.av.sort$ix],
        col =  c("royalblue","deeppink","orange","cyan","purple", "lightgreen","dodgerblue","yellow")[my.freq.table.av.sort$ix],
        ylab = "Cluster frequency in 10xGenomics single cell library",
        main = "BMDM clusters (res = 0.2)")
box()
dev.off()

#### get M3/EV ratios
my.freq.table.2 <- as.data.frame.matrix(my.freq.table)
my.ratio.table <- data.frame("YM_ratio" = my.freq.table.2$`4m_M_M3`/my.freq.table.2$`4m_M_CTL`,
                             "YF_ratio" = my.freq.table.2$`4m_F_M3`/my.freq.table.2$`4m_F_CTL`,
                             "OM_ratio" = my.freq.table.2$`20m_M_M3`/my.freq.table.2$`20m_M_CTL`,
                             "OF_ratio" = my.freq.table.2$`20m_F_M3`/my.freq.table.2$`20m_F_CTL`,
                             row.names = rownames(my.freq.table.2))
my.ratios <- data.frame(t(my.ratio.table))
colnames(my.ratios) <- paste0("Cluster_",0:7)

pdf(paste(Sys.Date(),"BMDM_Aging_10X_clusters_res_0.2_BOXPLOT_M3_over_CTL_frequency_change.pdf",sep = "_"), height = 10, width = 15)
boxplot(my.ratios, log = 'y', 
        col = c("royalblue","deeppink","orange","cyan","purple", "lightgreen","dodgerblue","yellow"),
        las = 2, ylab = "Ratio of cell proportion in M3/CTL")
beeswarm::beeswarm(my.ratios, add = T, pch = 16, cex= 1.5)
abline(h = 1, col = "red", lty = "dashed")
dev.off()


################################################################################################################################################################
##### 9. Visualization of potential marker gene expression

pdf(paste(Sys.Date(),"BMDM_aging_MOTSc_Feature_plot_Singlets_Mph_markers.pdf", sep = "_"))
FeaturePlot(object = my.bmdm.singlets, features = c("Adgre1","Itgam","Spi1","Cd14")  )
dev.off()

pdf(paste(Sys.Date(),"BMDM_aging_MOTSc_Xist_ridge_plot_Singlets.pdf", sep = "_"))
RidgePlot(object = my.bmdm.singlets, features = "Xist", ncol = 1, group.by = "Condition")
dev.off()

pdf(paste(Sys.Date(),"BMDM_aging_MOTSc_Y_genes_ridge_plot_Singlets.pdf", sep = "_"))
RidgePlot(object = my.bmdm.singlets, features = c("Ddx3y", "Eif2s3y"), ncol = 2, group.by = "Condition")
dev.off()

save(my.bmdm.singlets, file = paste(Sys.Date(),"BMDM_aging_MOTSc_10XGenomics_Singlets_Seurat_object.RData",sep = "_"))

library('scales')
my.umap.colors <- c(alpha("deeppink4"   , alpha = 0.3 ) ,
                    alpha("deepskyblue4", alpha = 0.3 ) ,
                    alpha("deeppink"    , alpha = 0.3 ) ,
                    alpha("deepskyblue" , alpha = 0.3 ) ,
                    alpha("salmon1"   , alpha = 0.3 ) ,
                    alpha("CornflowerBlue" , alpha = 0.3 ) ,
                    alpha("IndianRed2" , alpha = 0.3 ) ,
                    alpha("SteelBlue1" , alpha = 0.3 ) )


pdf(paste(Sys.Date(),"BMDM_Aging_Singlets_UMAP_color_by_Group.pdf", sep = "_"), height = 5, width = 6.5)
DimPlot(my.bmdm.singlets, reduction = "umap", group.by = "Condition", cols = my.umap.colors, pt.size	= 1)
dev.off()



################################################################################################################################################################
##### 10. DE-gene analysis


### Checking markers for clusters 5 and 6, (which seemed enriched upon MOTS-ctreatment) using resolution 0.2


my.bmdm.singlets <- SetIdent(my.bmdm.singlets, value = "SCT_snn_res.0.2")

###################################
clust5.marks <- FindMarkers(my.bmdm.singlets, ident.1 = "5", test.use = "wilcox")
sum(clust5.marks$p_val_adj < 0.05) #110
head(clust5.marks, n = 10)
#               p_val  avg_logFC pct.1 pct.2     p_val_adj
# Cd74      0.000000e+00  3.6420401 1.000 0.095  0.000000e+00
# H2-Aa     0.000000e+00  2.5411843 0.840 0.013  0.000000e+00
# H2-Eb1    0.000000e+00  2.1468638 0.670 0.006  0.000000e+00
# H2-Ab1    0.000000e+00  2.0737934 0.837 0.034  0.000000e+00
# Ciita     0.000000e+00  0.4367810 0.307 0.008  0.000000e+00
# H2-DMa   8.892002e-156  0.8223038 0.650 0.158 1.033428e-151
# Ccr2     1.012847e-134  0.4929264 0.221 0.019 1.177131e-130
# Itgb7     3.982868e-97  0.2601720 0.172 0.016  4.628889e-93
# Tmem176b  2.642708e-88  0.8305346 0.622 0.219  3.071355e-84
# Spp1      6.987731e-88 -1.3809487 0.943 0.989  8.121141e-84
write.table(clust5.marks, file = paste0(Sys.Date(),"_cluster5_MarkersWilcoxon_FDR5.txt"), sep = "\t", quote = F)


pdf(paste(Sys.Date(),"BMDM_aging_Volcano_Cluster5_Markers_Wilcoxon.pdf", sep = "_"))
EnhancedVolcano(clust5.marks, lab = rownames(clust5.marks),
                x = 'avg_logFC', y = 'p_val_adj', 
                pointSize = 3.0, labSize = 3.0, FCcutoff = 0,
                pCutoff = 0.05,
                pLabellingCutoff = 1e-150,
                xlim = c(-4,4), ylim = c(0,300))
dev.off()

pdf(paste(Sys.Date(),"BMDM_aging_Feature_plot_Singlets_Cluster5_Markers_Wilcoxon.pdf", sep = "_"))
FeaturePlot(object = my.bmdm.singlets, features = c("Cd74", "H2-Aa", "Ciita", "H2-Eb1")  )
dev.off()

###################################
clust6.marks <- FindMarkers(my.bmdm.singlets, ident.1 = "6", test.use = "wilcox")
sum(clust6.marks$p_val_adj < 0.05) #176
head(clust6.marks, n = 10)
# p_val avg_logFC pct.1 pct.2 p_val_adj
# Rsad2      0 1.9595793 0.761 0.089         0
# Ifit1      0 1.7779875 0.797 0.112         0
# Cxcl10     0 1.7325930 0.385 0.016         0
# Ifit3      0 1.4510152 0.821 0.069         0
# Ifit2      0 1.1982929 0.715 0.046         0
# Cmpk2      0 1.0787236 0.582 0.053         0
# Ifit3b     0 0.7958888 0.509 0.020         0
# Gbp2       0 0.7723744 0.500 0.027         0
# Ifi211     0 0.7410827 0.615 0.061         0
# Ms4a4c     0 0.6584047 0.364 0.021         0
write.table(clust6.marks, file = paste0(Sys.Date(),"_cluster6_MarkersWilcoxon_FDR5.txt"), sep = "\t", quote = F)

pdf(paste(Sys.Date(),"BMDM_aging_Volcano_Cluster6_Markers_Wilcoxon.pdf", sep = "_"))
EnhancedVolcano(clust6.marks, lab = rownames(clust6.marks),
                x = 'avg_logFC', y = 'p_val_adj', 
                pointSize = 3.0, labSize = 3.0, FCcutoff = 0,
                pCutoff = 0.05,
                pLabellingCutoff = 1e-150,
                xlim = c(-4,4), ylim = c(0,300))
dev.off()

pdf(paste(Sys.Date(),"BMDM_aging_Feature_plot_Singlets_Cluster6_Markers_Wilcoxon.pdf", sep = "_"))
FeaturePlot(object = my.bmdm.singlets, features = c("Rsad2", "Ifit1", "Cxcl10", "Cmpk2")  )
dev.off()

# get background
all.genes <- rownames(my.bmdm.singlets)
write.table(all.genes, file = paste0(Sys.Date(),"_All_BMDMs_scGenes_Background.txt"), sep = "\t", quote = F)


################################################################
#### get marker heatmap
bmdm.markers <- FindAllMarkers(my.bmdm.singlets, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10.bmdm.markers <- bmdm.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

png(paste0(Sys.Date(),"BMDM_aging_MOSTc_marker_0.2_heatmap.png"),width = 20, height = 18, units = "cm", res = 300)
DoHeatmap(my.bmdm.singlets, features = top10.bmdm.markers$gene)  + theme(axis.text.y = element_text(size = 5))
dev.off()

pdf(paste0(Sys.Date(),"BMDM_aging_MOSTc_marker_0.2_heatmap.pdf"), width = 10, height = 8)
DoHeatmap(my.bmdm.singlets, features = top10.bmdm.markers$gene)  + theme(axis.text.y = element_text(size = 8))
dev.off()


top10.5_6 <- top10.bmdm.markers[top10.bmdm.markers$cluster == c(5,6),]

pdf(paste0(Sys.Date(),"BMDM_aging_MOSTc_marker_0.2_Dtoplot_5_6.pdf"), width = 7, height = 5)
DotPlot(my.bmdm.singlets, idents = c(5,6), features = top10.5_6$gene, cols = c("blue", "red"), dot.scale = 8, split.by = "Treatment") + RotatedAxis()
dev.off()

pdf(paste0(Sys.Date(),"BMDM_aging_MOSTc_marker_0.2_Dtoplot_ALL.pdf"), width = 30, height = 10)
DotPlot(my.bmdm.singlets, features = top10.bmdm.markers$gene, cols = c("blue", "red"), dot.scale = 8, split.by = "Treatment") + RotatedAxis()
dev.off()

save(my.bmdm.singlets, file = paste(Sys.Date(),"BMDM_aging_MOTSc_10XGenomics_Singlets_Seurat_object.RData",sep = "_"))



# ################################################################################################################################################################
# ##### 11. CEMI tools
# library('CEMiTool')
# library(ggplot2)
# 
# # scale.data being pearson residuals; sctransform::vst intermediate results are saved in misc slot of new assay.
# my.bmdm.data <- as.data.frame(my.bmdm.singlets@assays$SCT@data)
# dim(my.bmdm.data)
# # [1] 11622 15415
# 
# # Adding sample annotation
# # More information can be included in CEMiTool to build a more complete object 
# # and generate richer reports about the expression data. Sample annotation can 
# # be supplied in a data.frame that specifies a class for each sample. 
# # Classes can represent different conditions, phenotypes, cell lines, time points, etc. 
# # load your sample annotation data
# sample_annot.1 <- my.bmdm.singlets@meta.data
# sample_annot.1$SampleName <- rownames(sample_annot.1)
# sample_annot <- sample_annot.1[,c("SampleName","Condition")]
# 
# # For ORA analysis
# gmt_in <- read_gmt("/Volumes/BB_Home_HQ/PATHWAY_ANNOT/ENSEMBL/2020-04-10_mouse_Ens99_GO_BP.gmt")
# 
# # Adding interactions
# int_df <- read.delim("/Volumes/BB_Home_HQ/PATHWAY_ANNOT/Interactions/STRING_2021-05-03/2021-05-03_STRING_Interaction_mouse.txt")
# 
# # run cemitool
# cem.bmdm <- cemitool(expr          = my.bmdm.data,
#                      annot         = sample_annot   ,
#                      gmt           = gmt_in         ,
#                      interactions  = int_df         ,
#                      cor_method    = "pearson"      ,
#                      class_column  = "Condition"    ,
#                      gsea_max_size = 5000           ,
#                      filter        = TRUE           ,
#                      plot          = TRUE           ,
#                      verbose       = TRUE           )
# 
# # create report as html document
# generate_report(cem.bmdm, directory = "./Report")
# 
# # write analysis results into files
# write_files(cem.bmdm   , directory = "./Tables")
# 
# # save all plots
# save_plots(cem.bmdm    , "all", directory = "./Plots")
# 



#######################
sink(file = paste(Sys.Date(),"_Seurat_session_Info.txt", sep =""))
sessionInfo()
sink()




