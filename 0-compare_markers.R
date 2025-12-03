# Rowan K. O'Hara-Payne
# 3/12/2025
# Compare seurat cluster marker genes to canonical markers


# Load libraries ----
library(Seurat)
library(tidyverse)


# Read in data ----
obj <- readRDS('/local/projects/idea/sament2/dkelly_tcell/clustering/dkelly_pilot.postqc.azimuth_celltypes.rds')

all_markers <- FindAllMarkers(obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

markers_sorted <- all_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group=T)

best_markers <- all_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5 & pct.1 > 0.5 & pct.2 < 0.2) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)  # Top 5 markers per cluster


# Visualize best markers ----
top_genes <- c("KLRF1", "CD28", "HLA-DQA1", "TSHZ2", "VCAN", "SLC4A10", "COL19A1", "BANK1", "NCAM1", "CCL5", "TOP2A", "TCF7L2", "RRM2", "CUX2", "IKZF2", "NKAIN2", "JCHAIN", "SSPN", "MIR924HG", "KIT", "PF4")

pdf("pbmc_violin_plots_dkelly.pdf", width = 10, height = 8)
VlnPlot(obj, features = top_genes, group.by = "predicted.celltype.l2", combine = FALSE)
dev.off()

for (gene in top_genes) {
  p <- FeaturePlot(obj, features = gene)
  ggsave(paste0(gene, "_FeaturePlot.pdf"), plot = p, width = 6, height = 5)
}

pdf("pbmc_dge_dkelly.pdf", width = 10, height = 8)
ggplot(best_markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.5) +
  theme_minimal()
dev.off()

pdf("pbmc_heatmap_dkelly.pdf", width = 10, height = 8)
DoHeatmap(obj, features = best_markers$gene, assay = "RNA", slot = "scale.data")
dev.off()


# Double Checking dnT annotations ----
obj <- readRDS("kelly_postqc_celltypes_meta.rds")
DefaultAssay(obj) <- 'RNA'
Idents(obj) <- obj$predicted.celltype.l2
t_cells_list <- c("CD4 Naive", "CD4 TCM", "CD4 Proliferating", "CD8 Naive", "CD8 TCM", "CD8 TEM", "dnT", "Treg", "NK", "NK Proliferating", "NK_CD56bright", "MAIT")
t_cells <- subset(obj, subset = predicted.celltype.l2 %in% t_cells_list)
table(t_cells$predicted.celltype.l2)

dnT_markers <- c("PTMA", "HSPH1", "SLC25A5", "HNRNPA1", "MTF2", "PSMA6", "BEX3", "MSI2", "UHRF1", "GAPDH", "JCHAIN", "TUBA1B", "TSC22D1", "HMGA1", "YWHAE", "MAP1A", "SIVA1", "SSBP2", "HSP90AB1", "PPIA", "ACTG1", "SMIM24", "STMN1", "MZB1", "CD99", "SOX4", "TFDP2", "GIHCG", "ADA", "HSP90AA1", "NUCB2", "NREP", "PTCRA", "TUBB", "HSPA1A", "HSPB1", "HSPA1B", "ARPP21", "HSPD1", "FXYD2", "RNASEH2B", "HMGB1", "CDK6", "HSPE1", "CD1E", "TKT", "CD28", "EEIG1", "ITGA6", "ENSG00000290067", "CCR7", "LCK", "TRABD2A")
markers <- c("CD3E", "CD3D", "TRAC", "TRBC1", "TRBC2", "IL7R", "CD2", "CD5", "CD4", "CD8A", "CD8B", "TRDC", "TRDV1", "TRGV9", "KLRB1", "ZBTB16", "GNLY", "NKG7", "KLRD1")


# dot plot of T cell markers
pdf("plots/dotplot_l2_Tcell_markers.pdf")
DotPlot(t_cells, features = markers) + 
	ggtitle("T Cell Markers") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# dot plot of dnT markers
pdf("plots/dotplot_l2_dnT_markers.pdf")
DotPlot(t_cells, features = dnT_markers) + 
	ggtitle("dnT Markers from TCell Atlas") + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



pdf("plots/dnT_markers_feature_dim.pdf")

DimPlot(t_cells, group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE) +
  ggtitle("Azimuth L2 annotations – T cell subsets") +
  theme(plot.title = element_text(size = 14, face = "bold"))

for (marker in markers) {
        p <- FeaturePlot(t_cells, features = marker, reduction = "umap", cols = c("lightgrey", "red"), order = TRUE) +
                theme_minimal() +
                labs(title = marker)
        print(p)
}


for (marker in dnT_markers) {
        p <- FeaturePlot(t_cells, features = marker, reduction = "umap", cols = c("lightgrey", "blue"), order = TRUE) +
                theme_minimal() +
                labs(title = marker)
        print(p)
}

for (marker in markers) {
        p <- VlnPlot(t_cells, features = marker, group.by = "predicted.celltype.l2", pt.size = 0) +
                ggtitle("TCell Markers: ", marker) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
}

for (marker in dnT_markers) {
        p <- VlnPlot(t_cells, features = marker, group.by = "predicted.celltype.l2", pt.size = 0) +
                ggtitle("dnT Markers: ", marker) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
}


dev.off()



## Trying Manual Annotation on the Seurat Clusters ----
obj <- readRDS("kelly_postqc_celltypes_meta.rds")
Idents(obj) <- obj$seurat_clusters


all_markers <- FindAllMarkers(obj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

markers_sorted <- all_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group=T)

best_markers <- all_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5 & pct.1 > 0.5 & pct.2 < 0.2) %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)  # Top 5 markers per cluster
write.csv(best_markers, "seurat_clusters_best_markers.csv")

markers <- c("CD3E", "CD3D", "TRAC", "TRBC1", "TRBC2", "IL7R", "CD2", "CD5", "CD4", "CD8A", "CD8B", "TRDC", "TRDV1", "TRGV9", "KLRB1", "ZBTB16", "GNLY", "NKG7", "KLRD1", "KLRF1", "CD28", "HLA-DQA1", "TSHZ2", "VCAN", "SLC4A10", "COL19A1", "BANK1", "NCAM1", "CCL5", "TOP2A", "TCF7L2", "RRM2", "CUX2", "IKZF2", "NKAIN2", "JCHAIN", "SSPN", "MIR924HG", "KIT", "PF4")

# UMAP of seurat clusters
pdf("plots/pbmc_umap_seurat_clusters.pdf", width=10, height=7)
DimPlot(obj,
	label = TRUE) +
        ggtitle("UMAP of Seurat Clusters")
dev.off()

# dot plot of PBMC markers
pdf("plots/dotplot_pbmc_markers_seurat_clusters.pdf", width=15, height=7)
DotPlot(obj, features = markers) +
        ggtitle("PBMC Markers") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# dot plot of FindAllMarkers
marker_list <- unique(best_markers$gene)

pdf("plots/dotplot_fam_seurat_clusters.pdf", width=20, height=7)
DotPlot(obj, features = marker_list) +
        ggtitle("FindAllMarkers() Markers") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# dot plot of dnT markers
pdf("plots/dotplot_dnT_markers_seurat_clusters.pdf", width=20, height=7)
DotPlot(obj, features = dnT_markers) +
        ggtitle("dnT Markers") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


