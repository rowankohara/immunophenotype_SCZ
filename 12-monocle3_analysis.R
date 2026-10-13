# Rowan O'Hara-Payne
# 9/11/2026
# Monocle analysis on monocytes
# conda activate /home/sonia.malaiya/.conda/envs/r_4.3/

# Load libraries and read in data ----
library(monocle3)
library(Seurat)
library(dplyr)

pbmc <- readRDS('kelly_postqc_clusters_meta.rds')


# Preprocess the data ----
mono <- subset(pbmc, subset = seurat_clusters %in% c(
  "0",
  "1",
  "10",
  "14", 
  "25", 
  "12", 
  "24-1"
))
mono$mono_type <- dplyr::recode(
  as.character(mono$seurat_clusters),
  "0"   = "SOD2hi CD14+ Mono",
  "1"   = "IL1Bhi CD4+ Mono",
  "10"  = "MHChi CD14+ Mono",
  "14"  = "CD14+ Mono",
  "25"  = "ISG+ CD14+ Mono",
  "12"  = "C1Qhi CD16+ Mono",
  "24-1" = "CD16+ Mono"
)


expr_matrix <- GetAssayData(
  mono,
  assay = "RNA",
  layer = "counts"
)
cell_metadata <- mono@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expr_matrix),
  row.names = rownames(expr_matrix)
)

cds <- new_cell_data_set(
  expression_data = expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, resolution=1e-5)

pdf("plots/monocle3_monocytes_cluster_cells.pdf")
plot_cells(cds)	 # community detection by monocle used for clusters
plot_cells(cds, label_cell_groups = FALSE)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")	 # partitions by monocle used for clusters
plot_cells(cds, color_cells_by="partition", group_cells_by="partition", label_cell_groups = FALSE)
plot_cells(cds, color_cells_by="mono_type")  # seurat clusters colors
plot_cells(cds, color_cells_by="mono_type", label_cell_groups = FALSE)
dev.off()

cds <- learn_graph(cds)

pdf("plots/monocle3_monocyte_trajectories.pdf")
plot_cells(
  cds,
  color_cells_by = "mono_type",
  label_cell_groups = TRUE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)
plot_cells(
  cds,
  color_cells_by = "mono_type",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)
dev.off()


root_cells <- colnames(cds)[
  colData(cds)$mono_type == "SOD2hi CD14+ Mono"
]

cds <- order_cells(
  cds,
  root_cells = root_cells
)

pdf("plots/monocle3_monocytes_pseudotime.pdf")
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE
)
dev.off()


# Compare AGA- and AGA+ ----

pt <- pseudotime(cds)
mono$pseudotime <- pt[colnames(mono)]

pdf("plots/pseudotime_violin_AGA_status.pdf")
VlnPlot(
  mono,
  features = "pseudotime",
  group.by = "Class"
)
dev.off()

# use Wilcoxon test to compare
wilcox.test(
  pseudotime ~ Class,
  data = mono@meta.data
)

# get sample median and compare
sample_pt <- mono@meta.data %>%
  group_by(sample, Class) %>%
  summarise(
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    .groups = "drop"
  )
wilcox.test(
  median_pseudotime ~ Class,
  data = sample_pt
)

pdf("plots/monocle3_AGA_status_trajectory.pdf")
plot_cells(
  cds,
  color_cells_by = "Class",
  label_cell_groups = FALSE
)
plot_cells(
  cds,
  color_cells_by = "mono_type",
  label_cell_groups = TRUE
)
dev.off()



# DEG ----

deg <- graph_test(
  cds,
  neighbor_graph = "principal_graph",
  cores = 4
)

deg <- deg %>%
  arrange(q_value)
head(deg)


