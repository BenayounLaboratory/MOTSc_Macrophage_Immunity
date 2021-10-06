setwd('/Users/berenice/Desktop/Single_cell/2020-04-18_BMDM_Single_Cell/Phantom_Purge_Preprocessing')
options(stringsAsFactors = F)

# purge index hopping
library("PhantomPurgeR")
library("furrr")
library(AnnotationHub)
library("stringr")
library(SingleCellExperiment)
library(tibble)
library(dplyr)
library(purrr)
library(scater)
library(Seurat)

# 2020-06-26
# analyze BMDM single cell 
# clean up barcode hopping

#####################################################################################################################
#### 0. Remove and purge phantom molecules/barcode hopping

# https://csglab.github.io/PhantomPurgeR/assets/notebooks/hiseq4000_vignette.html
# need the molecule_info.h5 files from all samples that were multiplexed

input_dir                <- "/Users/berenice/Desktop/Single_cell/2020-04-18_BMDM_Single_Cell/Cell_Ranger/Cell_Ranger_count/h5_input"
samples_filepaths        <- get_h5_filenames(input_dir)
names(samples_filepaths) <- unlist(strsplit(names(samples_filepaths) , ".molecule_info"))
samples_filepaths

### (1) Read data
out <- read10xMolInfoSamples(samples_filepaths)
# Loading samples data from molecule_info.h5 files produced by 10X Genomics CellRanger software: 105.926 sec elapsed

### (2) Join and merge data
out <- join_data(out)
# Joining data. Step 1: creating read counts datatables from lists: 0.177 sec elapsed
# Joining data. Step 2: joining and merging datatables for all samples keyed by cell, umi, and gene: 741.658 sec elapsed
# Joining data. Step 3: creating an outcome variable for each row: 1002.102 sec elapsed

# Read counts datatable
out$read_counts


### (3) Estimate sample index hopping rate
out <- estimate_hopping_rate(out)
# Estimating SIHR. Step 1: creating outcome counts datatable: 22.842 sec elapsed
# Estimating SIHR. Step 2: creating chimera counts datatable: 0.471 sec elapsed
# Estimating SIHR. Step 3: fitting GLM: 0.622 sec elapsed
# Estimating SIHR. Step 4: computing summary statistics: 1.123 sec elapsed

# Estimates
out$glm_estimates
# max_r  phat phat_low phat_high    SIHR   SBIHR
# <int> <dbl>    <dbl>     <dbl>   <dbl>   <dbl>
#   1   875 0.996    0.996     0.996 0.00367 0.00406

out$summary_stats$summary_estimates
# n_reads    n_cugs    n_pm      u        g   RMR p_chimeras
# <int>     <dbl>   <dbl>  <dbl>    <dbl> <dbl>      <dbl>
#   1 1090433345 151449310 3994683 0.0258 0.000539  7.20     0.0254
# p_chimeras is the proportion CUGs that are chimeric. 
# g is the estimated proportion of fugue molecules and u is the molecule inflation factor such that 
# n_cugs x u would give the number of non-fugue phantom molecules. 
# The estimated total number of phantom molecules present in the dataset is given by n_pm=n_cugs x (u+g).

# Marginal summary statistics
out$summary_stats$marginal
# A tibble: 8 x 6
# sample         n_reads prop_reads n_molecs p_molecs   RMR
# <chr>            <dbl>      <dbl>    <int>    <dbl> <dbl>
#   1 OF_CTL_BMDMs 142545359      0.131 21136136   0.136   6.74
# 2 OF_M3_BMDMs  125810722      0.115 17348039   0.112   7.25
# 3 OM_CTL_BMDMs 129026167      0.118 19838129   0.128   6.50
# 4 OM_M3_BMDMs  138537657      0.127 23380017   0.150   5.93
# 5 YF_CTL_BMDMs 134834499      0.124 15468086   0.0996  8.72
# 6 YF_M3_BMDMs  132411945      0.121 16145543   0.104   8.20
# 7 YM_CTL_BMDMs 144133212      0.132 18412417   0.119   7.83
# 8 YM_M3_BMDMs  143133784      0.131 23634059   0.152   6.06

### (4) Purge phantom molecules
# Call garbage collection first
gc()
# used    (Mb) gc trigger    (Mb) limit (Mb)   max used    (Mb)
# Ncells    9277217   495.5   17250137   921.3         NA   15384634   821.7
# Vcells 1833385240 13987.7 5162916586 39390.0     102400 5162916586 39390.0

# Set the trade-off ratio cost cutoff (torc). 
# The parameter torc represents the number of real molecules one is willing to incorrectly discard in order to correctly purge one phantom molecule. 
# Since discarding a large proportion of the data is undesirable, reasonable values of torc are expected to be within the range of 1-5.

out <- purge_phantoms(out, torc=2)
# Purging phantom molecules. Step 1: reassigning hopped reads: 17.159 sec elapsed
# Purging phantom molecules. Step 2: getting observed tor threshold below user-provided cutoff: 0.075 sec elapsed
# Purging phantom molecules. Step 3: marking retained observations in read counts table: 87.081 sec elapsed
# Purging phantom molecules. Step 4: creating sparse count matrices of cleaned and discarded data: 375.074 sec elapsed

out$summary_stats$cutoff_dt
# approach      outcome               s       qr       n    tor     o      FP      FN        TP      TN   FPm     FNm
# <chr>         <chr>             <int>    <dbl>   <int>  <dbl> <dbl>   <dbl>   <dbl>     <dbl>   <dbl> <dbl>   <dbl>
# 1 torc_before   8,0,0,11,0,0,10,0     4  0.00495       1  38.1  0.986   70937 2155208 149212535 3923746 55352 2110487
# 2 discard_torc  0,0,0,0,0,1,0,0       6  0.00512 2061428   1.33 0.999   81485  104328 151263415 3913198 44804   59607
# 3 torc_after    0,14,0,2,13,5,0,2     2  0.00594       1   1.33 0.999   81485  104327 151263416 3913198 44804   59606
# 4 no_discarding 0,1,0,0,1,0,1,0       7  0.649        19 NaN    1.     126289   44721 151323021 3868395     0       0
# 5 no_purging    NA                   NA NA            NA  NA    1.03  3994683       0 151367743       0    NA      NA

out$summary_stats$purge_summary
# sample       umi_total retained purged_real purged_phantom phantom_prop
# <chr>            <int>    <int>       <int>          <int>        <dbl>
# 1 OF_CTL_BMDMs  21136136 20646445       17035         472656       0.0224
# 2 OF_M3_BMDMs   17348039 16900496       10158         437385       0.0252
# 3 OM_CTL_BMDMs  19838129 19401627       12673         423829       0.0214
# 4 OM_M3_BMDMs   23380017 22902601       26450         450966       0.0193
# 5 YF_CTL_BMDMs  15468086 14837355        2765         627966       0.0406
# 6 YF_M3_BMDMs   16145543 15713558          58         431927       0.0268
# 7 YM_CTL_BMDMs  18412417 17919374        8773         484270       0.0263
# 8 YM_M3_BMDMs   23634059 23023444       26498         584117       0.0247

# Make diagnostic plots
dataset_name <- "BMDMs_HiseqXten"
plots <- make_plots(out, dataset_name, x_lim = 150, legend_rel_width=0.2)
plots 

# save
save(out, file = paste0(Sys.Date(),"_phantompurge_output_object.RData"))

# get help from developper on github
# sink("output.txt")
# str(out)
# sink()

#####  convert to SCE/Seurat
# https://github.com/csglab/PhantomPurgeR/issues/1
# create_sce_objects_with_metadata
# https://csglab.github.io/PhantomPurgeR/assets/notebooks/downstream_analysis_umap.html

# identify empty drops
# OF_CTL_BMDMs.empty <- call_emptydrops(out$umi_counts$retained$OF_CTL_BMDMs)
my.empty.OF_CTL_BMDMs  <- call_emptydrops(out$umi_counts$retained$OF_CTL_BMDMs)
my.empty.OF_M3_BMDMs   <- call_emptydrops(out$umi_counts$retained$OF_M3_BMDMs )
my.empty.OM_CTL_BMDMs  <- call_emptydrops(out$umi_counts$retained$OM_CTL_BMDMs)
my.empty.OM_M3_BMDMs   <- call_emptydrops(out$umi_counts$retained$OM_M3_BMDMs )
my.empty.YF_CTL_BMDMs  <- call_emptydrops(out$umi_counts$retained$YF_CTL_BMDMs)
my.empty.YF_M3_BMDMs   <- call_emptydrops(out$umi_counts$retained$YF_M3_BMDMs )
my.empty.YM_CTL_BMDMs  <- call_emptydrops(out$umi_counts$retained$YM_CTL_BMDMs)
my.empty.YM_M3_BMDMs   <- call_emptydrops(out$umi_counts$retained$YM_M3_BMDMs )

# keep cells (remove empty droplets)
# filter_cells <- function(mat, cell) { mat[, cell$is_cell] }
out$umi_counts$retained$OF_CTL_BMDMs <- out$umi_counts$retained$OF_CTL_BMDMs[, my.empty.OF_CTL_BMDMs$is_cell ]
out$umi_counts$retained$OF_M3_BMDMs  <- out$umi_counts$retained$OF_M3_BMDMs[,  my.empty.OF_M3_BMDMs$is_cell  ]
out$umi_counts$retained$OM_CTL_BMDMs <- out$umi_counts$retained$OM_CTL_BMDMs[, my.empty.OM_CTL_BMDMs$is_cell ]
out$umi_counts$retained$OM_M3_BMDMs  <- out$umi_counts$retained$OM_M3_BMDMs[,  my.empty.OM_M3_BMDMs$is_cell  ]
out$umi_counts$retained$YF_CTL_BMDMs <- out$umi_counts$retained$YF_CTL_BMDMs[, my.empty.YF_CTL_BMDMs$is_cell ]
out$umi_counts$retained$YF_M3_BMDMs  <- out$umi_counts$retained$YF_M3_BMDMs[,  my.empty.YF_M3_BMDMs$is_cell  ]
out$umi_counts$retained$YM_CTL_BMDMs <- out$umi_counts$retained$YM_CTL_BMDMs[, my.empty.YM_CTL_BMDMs$is_cell ]
out$umi_counts$retained$YM_M3_BMDMs  <- out$umi_counts$retained$YM_M3_BMDMs[,  my.empty.YM_M3_BMDMs$is_cell  ]

# Load Transcriptome
# find annotation hub id for transcriptome info
# mm10 <- query(ah, c("GRCm38","GRanges", "Mus Musculus"))
# cbind(mm10$ah_id,mm10$title)
# # [84,] "AH79245" "Mus_musculus.GRCm38.99.gtf"                               

gene_metadata <- get_transcriptome("AH79245")  #Mus_musculus.GRCm38.99
save(gene_metadata, file = paste0(Sys.Date(),"_AH79245_object.RData"))


# make into SCE objects
# OF_CTL_BMDMs.sce <- to_sce(out$umi_counts$retained$OF_CTL_BMDMs, "1")
# my.sce.obj <- map2(out$umi_counts$retained, seq_along(out$umi_counts$retained), to_sce)
my.ann.sce.obj <- create_sce(out$umi_counts$retained, seq_along(out$umi_counts$retained), gene_metadata)


# convert to Seurat
seurat.OF_CTL_BMDMs <- as.Seurat(my.ann.sce.obj$OF_CTL_BMDMs, counts = "counts", data = NULL, project = "10x_BMDMs")
seurat.OF_M3_BMDMs  <- as.Seurat(my.ann.sce.obj$OF_M3_BMDMs , counts = "counts", data = NULL, project = "10x_BMDMs")
seurat.OM_CTL_BMDMs <- as.Seurat(my.ann.sce.obj$OM_CTL_BMDMs, counts = "counts", data = NULL, project = "10x_BMDMs")
seurat.OM_M3_BMDMs  <- as.Seurat(my.ann.sce.obj$OM_M3_BMDMs , counts = "counts", data = NULL, project = "10x_BMDMs")
seurat.YF_CTL_BMDMs <- as.Seurat(my.ann.sce.obj$YF_CTL_BMDMs, counts = "counts", data = NULL, project = "10x_BMDMs")
seurat.YF_M3_BMDMs  <- as.Seurat(my.ann.sce.obj$YF_M3_BMDMs , counts = "counts", data = NULL, project = "10x_BMDMs")
seurat.YM_CTL_BMDMs <- as.Seurat(my.ann.sce.obj$YM_CTL_BMDMs, counts = "counts", data = NULL, project = "10x_BMDMs")
seurat.YM_M3_BMDMs  <- as.Seurat(my.ann.sce.obj$YM_M3_BMDMs , counts = "counts", data = NULL, project = "10x_BMDMs")

# Merge Seurat objects
bmdm.combined <- merge(seurat.YF_CTL_BMDMs, 
                       y =  c(seurat.OF_CTL_BMDMs,
                              seurat.YM_CTL_BMDMs,
                              seurat.OM_CTL_BMDMs,
                              seurat.YF_M3_BMDMs ,
                              seurat.OF_M3_BMDMs ,
                              seurat.YM_M3_BMDMs ,
                              seurat.OM_M3_BMDMs   ), 
                       add.cell.ids = c("4m_F_CTL",
                                        "20m_F_CTL",
                                        "4m_M_CTL",
                                        "20m_M_CTL",
                                        "4m_F_M3",
                                        "20m_F_M3",
                                        "4m_M_M3",
                                        "20m_M_M3"), 
                       project = "10x_BMDMs")

bmdm.combined@meta.data <- bmdm.combined@meta.data[,-c(4:11)]

############################## metadata ############################## 
# create condition label
my.4mF.c  <- grep("4m_F_CTL", colnames(bmdm.combined@assays$RNA))
my.4mM.c  <- grep("4m_M_CTL", colnames(bmdm.combined@assays$RNA))
my.20mF.c <- grep("20m_F_CTL", colnames(bmdm.combined@assays$RNA))
my.20mM.c <- grep("20m_M_CTL", colnames(bmdm.combined@assays$RNA))
my.4mF.m  <- grep("4m_F_M3", colnames(bmdm.combined@assays$RNA))
my.4mM.m  <- grep("4m_M_M3", colnames(bmdm.combined@assays$RNA))
my.20mF.m <- grep("20m_F_M3", colnames(bmdm.combined@assays$RNA))
my.20mM.m <- grep("20m_M_M3", colnames(bmdm.combined@assays$RNA))


Condition <- rep("NA", length(colnames(bmdm.combined@assays$RNA)))
Condition[ my.4mF.c  ]   <- "4m_F_CTL"
Condition[ my.4mM.c  ]   <- "4m_M_CTL"
Condition[ my.20mF.c ]   <- "20m_F_CTL"
Condition[ my.20mM.c ]   <- "20m_M_CTL"
Condition[ my.4mF.m  ]   <- "4m_F_M3"
Condition[ my.4mM.m  ]   <- "4m_M_M3"
Condition[ my.20mF.m ]   <- "20m_F_M3"
Condition[ my.20mM.m ]   <- "20m_M_M3"

Condition <- data.frame(Condition)
rownames(Condition) <- colnames(bmdm.combined@assays$RNA)

bmdm.combined <- AddMetaData(object = bmdm.combined, metadata = as.vector(Condition), col.name = "Condition")


# create treatment label
my.ctl  <- grep("_CTL", colnames(bmdm.combined@assays$RNA))
my.m3   <- grep("_M3", colnames(bmdm.combined@assays$RNA))

Treatment <- rep("NA", length(colnames(bmdm.combined@assays$RNA)))
Treatment[ my.ctl ]   <- "CTL"
Treatment[ my.m3  ]   <- "M3"

Treatment <- data.frame(Treatment)
rownames(Treatment) <- colnames(bmdm.combined@assays$RNA)

bmdm.combined <- AddMetaData(object = bmdm.combined, metadata = as.vector(Treatment), col.name = "Treatment")


### make Seurat objects for the entire dataset or just the control
my.bmdm.all <- bmdm.combined
my.bmdm.ctl <- subset(bmdm.combined, subset = Treatment == "CTL")

save(my.bmdm.all, file = paste0(Sys.Date(),"_PhantomPurged_BMDM_CTL_M3_Seurat.RData"))
save(my.bmdm.ctl, file = paste0(Sys.Date(),"_PhantomPurged_BMDM_CTL_only_Seurat.RData"))


#######################
sink(file = paste(my.outprefix,"_PhantomPurge_session_Info.txt", sep =""))
sessionInfo()
sink()


