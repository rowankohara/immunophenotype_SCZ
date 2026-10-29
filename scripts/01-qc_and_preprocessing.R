# Rowan K. O'Hara-Payne
# 4/8/2025
# conda activate /local/projects-t3/idea/amentlab_software/conda_envs/Seurat5
# Doublet detection for the D Kelly immunophentype of schizophrenia and healthy controls
# Base on cluster.R and biostatsquid code (https://biostatsquid.com/doubletfinder-tutorial/)


# Load libraries and data ----
library(Seurat)
library(harmony)
library(SeuratWrappers)
library(DoubletFinder)
library(dplyr)
library(tibble)
library(ggplot2)

# create Seurat object
mtx <- ReadParseBio('/local/projects/idea/sament2/dkelly_tcell/processed/output_combined/all-sample/DGE_filtered')
cell_meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/processed/output_combined/all-sample/DGE_filtered/cell_metadata.csv')
rownames(cell_meta) <- cell_meta$bc_wells
obj <- CreateSeuratObject(mtx, meta.data = cell_meta)


# Initial QC ----
obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = '^MT-')
obj <- subset(obj, nCount_RNA < 10000 & nCount_RNA > 500 & percent.mt < 15)

# remove participants with overall low quality cells
obj <- subset(obj, sample %in% c('M094746_M_40yo_HC',
				 'M095130_M_40yo_SCZ',
				 'M098087_F_36yo_SCZ'),
	      invert = T)


cat("sample\tmin_pc\toptimal_pK\tnum_cells\texp_doublets\thomotypic_prop\tadj_exp_doublets", "\n", file='output_doublet_detection.txt')

# Doublet detection using doubletFinder ---- 
# based on biostatsquid function (https://github.com/biostatsquid/scripts/blob/main/Analysis/scRNAseq/1_Preprocessing/scRNAseq_doublet_detection.R)
run_doubletfinder <- function(sample_subset) {
	

	cat(unique(sample_subset@meta.data$sample), "\t", file='output_doublet_detection.txt', sep="", append=T)

	# pre-processing Seurat object
	sample <- NormalizeData(sample_subset)
	sample <- FindVariableFeatures(sample)
	sample <- ScaleData(sample)
	sample <- RunPCA(sample, nfeatures.print = 10)

	# find significant PCs
	stdv <- sample[["pca"]]@stdev
	percent_stdv <- (stdv/sum(stdv)) * 100
	cumulative <- cumsum(percent_stdv)
	co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
	co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] -
			   percent_stdv[2:length(percent_stdv)]) > 0.1), decreasing = T)[1] + 1
	min_pc <- min(co1, co2)

	cat(min_pc, "\t", file='output_doublet_detection.txt', sep="", append=T)

	# finish pre-processing with min_pc
	sample <- RunUMAP(sample, dims = 1:min_pc)
	sample <- FindNeighbors(sample, dims = 1:min_pc)
	sample <- FindClusters(sample, resolution = 0.4)  # Satija lab suggest 0.4-1.2 for overall clustering analysis for ~3k cells

	# pK identification
	sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)
	sweep_stats <- summarizeSweep(sweep_list)
	bcmvn <- find.pK(sweep_stats)
	optimal.pk <- bcmvn %>%
		dplyr::filter(BCmetric == max(BCmetric)) %>%
		dplyr::select(pK)
	optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))

	cat(optimal.pk, "\t", file='output_doublet_detection.txt', sep="", append=T)


	# estimate expected doublets
	num_cells <- nrow(sample@meta.data)
	expected_doublet_rate <- ifelse(
		num_cells <= 10000, 0.02,
		ifelse(num_cells <= 20000, 0.03,
		ifelse(num_cells <= 50000, 0.05,
		ifelse(num_cells <= 1000000, 0.07, 0.20)))
	)
	expected_doublets <- round(expected_doublet_rate * num_cells)

	cat(num_cells, "\t", expected_doublets, "\t", file='output_doublet_detection.txt', sep="", append=T)

	# homotypic doublet proportion estimate
	annotations <- sample@meta.data$seurat_clusters
	homotypic.prop <- modelHomotypic(annotations)
	nExp.poi.adj <- round(expected_doublets * (1 - homotypic.prop))

	cat(homotypic.prop, "\t", file='output_doublet_detection.txt', sep="", append=T)
	cat(nExp.poi.adj, "\n", file='output_doublet_detection.txt', sep="", append=T)

	# run DoubletFinder
	sample <- doubletFinder(seu = sample,
				PCs = 1:min_pc,
				pK = optimal.pk,
				nExp = nExp.poi.adj)
	colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"

	doublet_finder_res <- sample@meta.data['doublet_finder']
	doublet_finder_res <- rownames_to_column(doublet_finder_res, "row_names")
	return(doublet_finder_res)
}

# Run DoubletFinder -----
obj_split <- SplitObject(obj, split.by = "sample")
obj_split <- lapply(obj_split, run_doubletfinder)

sglt_dblt_metadata <- data.frame(bind_rows(obj_split))  # merge to a single dataframe
rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names  # assign cell IDs to row names to ensure match
sglt_dblt_metadata$row_names <- NULL
obj <- AddMetaData(obj, sglt_dblt_metadata, col.name = 'doublet_finder')

# visualize features between singlets and doublets per sample
pdf("doublet_detection_results.pdf")
VlnPlot(obj, group.by = 'sample', split.by = 'doublet_finder',
	features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
	combine = FALSE, pt.size = 0)
dev.off()


# get doublets per sample
doublets_summary <- obj@meta.data %>%
	group_by(sample, doublet_finder) %>%
	summarize(total_count = n(), .groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
	group_by(sample) %>%
	mutate(countT = sum(total_count)) %>%
	group_by(doublet_finder, .add = TRUE) %>%
	mutate(percent = paste0(round(100 * total_count/countT, 2), '%')) %>%
	dplyr::select(-countT)
write.table(doublets_summary, file = 'doubletFinder_doublets_summary.txt', quote = FALSE, row.names = FALSE, sep = '\t')


# save data
saveRDS(obj, file='dkelly_doubletFinder_results.rds')

singlets <- subset(obj, doublet_finder == 'Singlet')
saveRDS(singlets, file='dkelly_doubletFinder_singlets.rds')


# Preprocess on the integrated object ----
singlets <-  NormalizeData(singlets)
singlets <-  FindVariableFeatures(singlets,
                           selection.method = "vst",
                           nfeatures = 5000 )
singlets <-  ScaleData(singlets)

# run PCA and Harmony
singlets <-  RunPCA(singlets)
singlets <-  RunHarmony(singlets, group.by.vars = 'sample' )

# Run UMAP and clustering
singlets <- RunUMAP(singlets, reduction = "harmony", dims = 1:30)
singlets <- FindNeighbors(singlets, reduction = "harmony", dims = 1:30)
singlets <- FindClusters(singlets, resolution = 0.5)
singlets <- FindClusters(singlets, resolution = 1)
singlets <- FindClusters(singlets, resolution = 1.5)


# Save RDS
saveRDS(singlets, file='dkelly_singlets.postqc.harmony.30pcs.rds')


# plots
pdf('dkelly_singlets.all-sample.harmony.30pcs.umaps.pdf')
DimPlot(singlets, label = T , group.by = 'RNA_snn_res.0.5' , raster = T ) + NoLegend()
DimPlot(singlets, label = T , group.by = 'RNA_snn_res.1' , raster = T ) + NoLegend()
DimPlot(singlets, label = T , group.by = 'RNA_snn_res.1.5' , raster = T ) + NoLegend()
DimPlot(singlets, label = T , group.by = 'sample' , raster = T ) + NoLegend()
FeaturePlot(singlets, features = 'nCount_RNA' )
FeaturePlot(singlets, features = 'percent.mt' )
dev.off()

