# phantom purge accessory function ### https://csglab.github.io/PhantomPurgeR/assets/notebooks/workflow_hiseq4000.nb.html

########################################################
# convert a *single* sample to Single Cell Experiment format
# modified from https://csglab.github.io/PhantomPurgeR/assets/notebooks/downstream_analysis_umap.html

create_sce <- function(data_list, suffixes, gene_metadata) {
  data_list <- map2(data_list, suffixes, to_sce)
  data_list <- map(data_list, add_metadata, gene_metadata)
  return(data_list)
}


########################################################
# convert a *single* sample to Single Cell Experiment format
# https://github.com/csglab/PhantomPurgeR/issues/2
to_sce <- function(data, suff) {
  col_names <- paste(colnames(data), suff, sep = "_")
  row_names <- rownames(data)
  
  colnames(data) <- NULL
  rownames(data) <- NULL
  
  data <- SingleCellExperiment(assays = list(counts = data),
                               rowData = DataFrame(gene_id = row_names),
                               colData = DataFrame(barcode = col_names)
  )
  return(data)
}


########################################################
# get metadata information
# modified from https://csglab.github.io/PhantomPurgeR/assets/notebooks/downstream_analysis_umap.html
add_metadata <- function(sce10x, gene_metadata) {
  rowData(sce10x) <- gene_metadata %>% left_join(as_tibble(rowData(sce10x)),
                                                 .,
                                                 by = "gene_id") %>% dplyr::rename(chr = seqnames)
  rownames(sce10x) <- uniquifyFeatureNames(rowData(sce10x)$gene_id, rowData(sce10x)$gene_name)
  sce10x <- addPerCellQC(sce10x, subsets = list(mito = which(rowData(sce10x)$chr == "chrM")), flatten = TRUE)
  sce10x <- addPerFeatureQC(sce10x)
  
  return(sce10x)
}

########################################################
# get transcriptome information
# modified from https://csglab.github.io/PhantomPurgeR/assets/notebooks/downstream_analysis_umap.html
get_transcriptome <- function(ah_id) {
  ah <- AnnotationHub()
  gene_metadata <- ah[[ah_id]]
  gene_metadata <-
    gene_metadata[mcols(gene_metadata)$type == "gene",c("gene_id", "gene_name", "gene_biotype")]
  seqlevelsStyle(gene_metadata) <- "UCSC"
  
  gene_metadata  <- gene_metadata %>% as.data.frame() %>% dplyr::select(gene_id,
                                                                        seqnames,
                                                                        gene_name,
                                                                        gene_biotype,
                                                                        gene_id)
  
  return(gene_metadata)
  
}


########################################################
# call empty droplets to remove 
#  from https://csglab.github.io/PhantomPurgeR/assets/notebooks/downstream_analysis_umap.html
call_emptydrops <- function(sample_umi_count, lower = 200, fdr_thresh = 0.005) {
  sample_emptydrops <- DropletUtils::emptyDrops(sample_umi_count,
                                                lower = lower) %>% as.data.frame() %>% rownames_to_column(var = "cell") %>% select(cell, FDR)
  
  sample_emptydrops[is.na(sample_emptydrops)] <- 1
  
  sample_emptydrops <- sample_emptydrops %>% mutate(is_cell = FDR <= fdr_thresh, FDR = NULL)
  
  return(sample_emptydrops)
}
