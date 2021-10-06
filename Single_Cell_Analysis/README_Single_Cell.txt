############   README   ############

Scripts for analysis of single-cell RNA-seq of mouse BMDM differentiated in the absence or presence of MOTSc

- Cell_Ranger: folder with cell-ranger mapping scripts
		* cell_ranger_count_BMDM.sh: running cellranger count on all 8 libraries

- PhantomPurge: removing chimeric molecules due to Illumina patterned flow cells
		* 2020-06-30_purge_phantoms_BMDM_aging_scRNAseq.R: script for phantompurge processing
		* Phantompurge_accessory_functions.R: accessory functions for phantompurge processing
		
- Seurat_and_enrichment: 
		* 2021-06-24_analyze_BMDM_aging_scRNAseq.R: script for Seurat data processing
		* run_GO_analysis.R: script for GO enrichment analysis on cluster markers
