# Rowan K. O'Hara-Payne
# 3/12/2025
# Compare seurat cluster marker genes to canonical markers


# Load libraries ----
library(Seurat)
library(dplyr)


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


