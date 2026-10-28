library(miloR)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)

obj <- readRDS("GSE252762_kornberg_han_celiac_dataset_fig1-2.RDS")
obj_sce <- as.SingleCellExperiment(obj)
obj_milo <- Milo(obj_sce)


# get the UMAP from the single cell experiment object
reducedDim(obj_milo, "UMAP") <- reducedDim(obj_sce, "UMAP")
reducedDim(obj_milo, "PCA") <- Embeddings(obj, "pca")[, 1:30]

# plot variancce explained by PCs -- will help choose our d parameter
pdf("plots/elbow_plot_PCs_kornberg.pdf")
ElbowPlot(obj, ndims = 50)
dev.off()


obj_milo_conserv <- buildGraph(obj_milo, k = 200, d = 20)
	
obj_milo_conserv <- makeNhoods(obj_milo_conserv, prop = 0.2, k = 200, d = 20, refined = TRUE, refinement_scheme="graph")
obj_milo_conserv <- countCells(obj_milo_conserv, meta.data = data.frame(colData(obj_milo_conserv)), samples = "sample")

pdf("plots/NHoodSize_milo_kornberg_200.pdf")
p1 <- plotNhoodSizeHist(obj_milo_conserv)
print(p1)
dev.off()

obj_design <- data.frame(colData(obj_milo_conserv))[, c("sample", "general_condition", "sample_sex", "sample_age")]
obj_design <- distinct(obj_design)
rownames(obj_design) <- obj_design$sample
obj_design <- obj_design[colnames(nhoodCounts(obj_milo_conserv)), , drop=FALSE]


model <- model.matrix(~ 0 + general_condition + sample_sex + sample_age, data=obj_design)

obj_milo_conserv <- calcNhoodDistance(obj_milo_conserv, d=20)
contrastCvG <- testNhoods(obj_milo_conserv, design=model, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("general_conditionCELIAC - general_conditionCONTROL"))

obj_milo_conserv <- buildNhoodGraph(obj_milo_conserv)
pdf("plots/NhoodGraph_milo_kornberg_200.pdf")
p2 <- plotNhoodGraphDA(obj_milo_conserv, contrastCvG, alpha=0.15) + 
	labs(capton = "k = 200, d = 20, SpatialFDR < 0.15")
p22 <- plotNhoodGraphDA(obj_milo_conserv, contrastCvG, alpha=0.30) + 
	labs(caption = "k = 200, d = 20, SpatialFDR < 0.30")	
print(p2)
print(p22)
dev.off()

annot_df <- annotateNhoods(obj_milo_conserv, contrastCvG, "clusters")
write_csv(annot_df, paste0("milo_kornberg_results_summary_200.csv"))

annot_df$clusters <- ifelse(annot_df$clusters_fraction < 0.7, "Mixed", annot_df$clusters)
	
pdf("plots/DA_celltype_frac_neighborhoods_kornberg_200.pdf")
p3 <- ggplot(annot_df, aes(clusters_fraction)) + geom_histogram(bins=50)
p4 <- plotDAbeeswarm(annot_df, group.by = "clusters", alpha=0.30) + 
	labs(x = "Cell Type",
	     y = "logFC",
	     caption = "k = 200, d = 20, SpatialFDR < 0.30") +
	theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
p42 <- plotDAbeeswarm(annot_df, group.by = "clusters", alpha=0.15) + 
	labs(x = "Cell Type",
	     y = "logFC",
	     caption = "k = 200, d = 20, SpatialFDR < 0.15") +
	theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
print(p3)
print(p4)
print(p42)
dev.off()

pdf("plots/DA_results_milo_kornberg_200.pdf")
p5 <- ggplot(annot_df, aes(PValue)) + geom_histogram(bins=50)
p6 <- ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
  	geom_point(aes(color = clusters)) +
  	geom_hline(aes(yintercept = -log10(0.15), linetype = "spatialFDR < 0.15")) + ## Mark significance threshold (15% FDR)
  	geom_hline(aes(yintercept = -log10(0.3), linetype = "spatialFDR < 0.30")) +
	scale_linetype_manual(name = "threshold", values = c("spatialFDR < 0.15" = "solid", "spatialFDR < 0.30" = "dashed")) +
	theme_classic() +
        labs(color = "cell type") +
        theme(legend.text = element_text(size = 10),
              legend.title = element_text(size = 14),
              axis.title = element_text(size = 14))
print(p5)
print(p6)
dev.off()

nhood.dge <- findNhoodMarkers(obj_milo_conserv, contrastCvG, overlap=1, compute.new=TRUE)

contrastCvG <- groupNhoods(obj_milo_conserv, contrastCvG, max.lfc.delta = 1)
nhood.dge <- testDiffExp(x=obj_milo_conserv, da.res=contrastCvG, model.contrasts = c("general_conditionCELIAC - general_conditionCONTROL"), design=~ 0 + general_condition + sample_sex + sample_age, meta.data=obj_design)

saveRDS(nhood.dge, "milo_kornberg_200_markers.RDS")  # save object

