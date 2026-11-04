# Deepak Salem
# 11/8/2025
# Subclustering Seurat clusters to further manually annotate and verify cell types
####################################################################################################################


# Load libraries and data ----

library(Seurat)
library(dplyr)
library(ggplot2)

s <- readRDS("dkelly_singlets.postqc.harmony.30pcs.rds")


# Subclustering for cell type annotation ----

## 1) Subset to global cluster 13 ----

DefaultAssay(s) <- "RNA"
Idents(s) <- "seurat_clusters"

cells13 <- WhichCells(s, idents = "13")
s_cl13  <- subset(s, cells = cells13)   # new object, won't overwrite s


## 2) Standard Seurat subclustering pipeline on cluster 13 ----

set.seed(1234)  # for reproducibility

# Normalize (within this subset)
s_cl13 <- NormalizeData(
  s_cl13,
  normalization.method = "LogNormalize",
  scale.factor         = 1e4
)

# Variable features
s_cl13 <- FindVariableFeatures(
  s_cl13,
  selection.method = "vst",
  nfeatures        = 2000
)

# Scale
s_cl13 <- ScaleData(
  s_cl13,
  features = VariableFeatures(s_cl13)
)

# PCA
s_cl13 <- RunPCA(
  s_cl13,
  features = VariableFeatures(s_cl13),
  npcs     = 20,
  verbose  = FALSE
)

# (Optional) elbow plot to decide how many PCs
ElbowPlot(s_cl13)

# Choose number of PCs to use (edit if you like)
dims_use <- 1:30

# Neighbors + clustering
s_cl13 <- FindNeighbors(s_cl13, dims = dims_use)
s_cl13 <- FindClusters(s_cl13, resolution = 0.2)  # try 0.2, 0.6, etc. as well

# UMAP
s_cl13 <- RunUMAP(s_cl13, dims = dims_use)

Idents(s_cl13) <- "seurat_clusters"


## 3) Visualize subclusters ----

# UMAP of cluster-13 subclusters
DimPlot(
  s_cl13,
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.5
) +
  ggtitle("Cluster 13 – subclusters 13-0 and 13-1")

# Optional: see how SCZ High / Low / HC sit in these subclusters
DimPlot(
  s_cl13,
  reduction = "umap",
  group.by  = "Group",
  pt.size   = 0.3
) +
  ggtitle("Cluster 24 – Group distribution (SCZ High / SCZ Low / HC)")


##
markers <- intersect(unique(chosen_markers), rownames(s_sub_combined))

dot_raw <- DotPlot(
  s_sub_combined,
  features = markers,
  group.by = "combined_cluster_13"
) + Seurat::RotatedAxis()

df_dot  <- dot_raw$data
rng     <- range(df_dot$avg.exp.scaled, na.rm = TRUE)
max_abs <- max(abs(rng))

dot_plot <- dot_raw +
  scale_color_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    limits   = c(-max_abs, max_abs),
    name     = "Scaled avg exp\n(z per gene)"
  ) +
  scale_size(range = c(0, 6), name = "% expressing") +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 8),
    panel.grid   = element_blank()
  ) +
  labs(
    x = "Markers",
    y = "Seurat Cluster",
    title = "T Cell Clusters" 
  )

print(dot_plot)


# <<< EDIT THIS VECTOR: choose and order the clusters as you want on the y-axis >>>
#clusters_to_plot <- c("3", "4", "16", "6", "7", "13-0", "8", "9", "17", "20", "22", "13-1", "31")

clusters_to_plot <- c("3", "5", "16", "6", "7", "20", "13-0", "31", "4", "8", "9", "17","22", "28", "13-1") 
# Subset s_sub_combined to just the clusters you want
cells_use <- colnames(s_sub_combined)[
  s_sub_combined$combined_cluster_13 %in% clusters_to_plot
]

s_dot <- subset(s_sub_combined, cells = cells_use)

# Ensure the cluster factor has the desired order for the dotplot y-axis
s_dot$combined_cluster_13 <- factor(
  s_dot$combined_cluster_13,
  levels = clusters_to_plot
)

# <<< EDIT THIS: your marker list >>>
chosen_markers <- c(
# --- CD4 / CD8 ---
"CD4", "CD8A",

"KLRB1", "ZBTB16", "SLC4A10", 

# --- Naive / memory / survival scaffold ---
"IL6ST", "FOXO1", "TCF7", "CD27", "IL7R", "CD28", "LEF1", "FAS", "CXCR4", "CD44", "IL2RA", "IL2RB", "KLF2", "DPP4",

#Migration
"CCR7", "SELL", "CD69", "ITGA1","ITGA4","ITGAE","ITGB1","ITGB7", "ABCB1",

# --- Function / polarization TFs ---
"CXCR3", "GATA3", "CCR4", "PTGDR2", "BCL6", "MAF", "ICOS", "CXCR5", "IKZF2", "FOXP3", "CTLA4", "TIGIT", "TNFRSF1B", "AHR", "CCR6", "RORC", "IL23R", "KIT", "IL10RA", "IL18R1", "IL18RAP", "NCAM1", "B3GAT1", "CD160", "EOMES", "TBX21", "IFNG", "CCR5", "CX3CR1", 

"NCAM1", "FCGR3A", "FCER1G", "KLRK1", "GZMK", "KLRC1", "KLRC2", "PRF1", "KLRG1", "GZMB", "GNLY", "FGFBP2", "ZEB2", "PRDM1",



# --- T cell lineage / identity (core) ---
"CD28","ZAP70",
"KLRB1","SLC4A10","ZBTB16",
"NCAM1","FCGR3A","KLRK1","CD160",
"KIR3DL1","KIR2DL3","KIR2DL4","KIR2DS4","KIR3DL2","KIR3DL3",
"KLRC1","KLRC2","B3GAT1","CX3CR1","ITGB2","FCER1G","TYROBP",


# --- TCR ---
"TRAC","TRBC1","TRBC2","TRDC","TRGC1",

# --- Regulatory / checkpoint / suppression ---
"FOXP3","IL2RA","CTLA4","TIGIT","TNFRSF1B","FGL2","IL10RA",

# --- TCR ---
"TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2",

# --- Cytotoxic / effector molecules ---
"GZMK","GZMB","PRF1","GNLY","FGFBP2","IFNG","FASLG","ZEB2",

# --- Migration / homing / tissue positioning ---
"CCR4","CCR5","CCR6","CCR7","CCR8","CCR10",
"CXCR3","CXCR4","CXCR5","CXCR6",
"S1PR1","ABCB1",

# --- Adhesion / integrins ---
"ITGA1","ITGA4","ITGAE","ITGB1","ITGB7",

# --- Naive / memory / survival scaffold ---
"IL7R","KLF2","BTG2",

# --- Activation / immediate-early / stress / antigen presentation ---
"NR4A2","DUSP1","HLA-DRA","CD38","ISG15",

# --- State markers (activation, memory, senescence) ---
"CD44","CD69","FAS",

# --- Helper/Tfh interaction markers ---
"ICOS","CXCL13","CD84","CD200",

# --- Cytokine / lipid receptors relevant to helper programs ---
"IL4R","IL2RB","IL18R1","IL18RAP","IL23R","DPP4","PTGDR2","KIT"
)

# Keep only genes present in this subset
markers <- intersect(unique(chosen_markers), rownames(s_dot))
markers

DefaultAssay(s_dot) <- "RNA"

dot_raw <- DotPlot(
  s_dot,
  features = markers,
  group.by = "combined_cluster_13",  # uses the combined cluster label
  scale    = TRUE                    # z-scales avg.exp per gene across these clusters
) + Seurat::RotatedAxis()

# Optionally polish the color scale a bit
df_dot  <- dot_raw$data
rng     <- range(df_dot$avg.exp.scaled, na.rm = TRUE)
max_abs <- max(abs(rng))

dot_plot <- dot_raw +
  scale_color_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0,
    limits   = c(-max_abs, max_abs),
    name     = "Scaled avg exp\n(z per gene)"
  ) +
  scale_size(range = c(0, 6), name = "% expressing") +
  theme(
    axis.text.x  = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 8),
    panel.grid   = element_blank()
  ) +
  labs(
    x = "Marker Genes",
    y = "Seurat Cluster",
    title = "T Cell Subsets"
  )

print(dot_plot)

##
# You can edit this marker list as you like
markers <- c(
  "CD4", "CD8A", "CD8B", "TRAC", "TRDC", "TRBC2", "TRGC2" 
  # CD4 vs CD8              # TCR alpha constant
)

# Keep only genes present in this subset
markers <- intersect(markers, rownames(s_cl22))
markers

FeaturePlot(
  s_cl22,
  features  = markers,
  reduction = "umap",
  ncol      = 3,
  cols      = c("lightgrey", "red")
)

DefaultAssay(s_cl14) <- "RNA"

## Repeat for other clusters ----
# New object just for cluster 21
s_cl24 <- subset(s, cells = cells24)
s_cl24
DefaultAssay(s_cl24) <- "RNA"
Idents(s_cl24) <- "seurat_clusters"   # should be 0,1,2 after subclustering

# UMAP of the 3 subclusters within cluster 21
DimPlot(
  s_cl24,
  reduction = "umap",
  label     = TRUE,
  repel     = TRUE,
  pt.size   = 0.6
) +
  ggtitle("Cluster 21 – Seurat subclusters (0,1,2)")

# T-cell markers
t_markers_19 <- c(
  "CD3D","CD3E","CD3G",   # TCR/CD3 complex
  "TRAC",                 # TCRα constant
  "CD4","CD8A",
  "MS4A1", "IGKC"
  # coreceptors
)

# Monocyte markers
mono_markers_19 <- c(
  "CD14","LYZ",           # classical mono
  "S100A8","S100A9",      # inflammatory mono
  "LST1","VCAN","CTSS",   # myeloid signature
  "FCGR3A"                # CD16
)
# Monocyte markers
mono_markers_19 <- c(
  "CD34"                # CD16
)

markers_19 <- c(t_markers_19, mono_markers_19)

# Keep only genes that actually exist in this subset
markers_19 <- intersect(markers_19, rownames(s_cl24))
markers_19

FeaturePlot(
  s_cl24,
  features  = markers_19,
  reduction = "umap",
  ncol      = 3,
  cols      = c("lightgrey", "red")
)


# Plot dot plots ----

cluster_labels_ <- c(
  "21-2" = "Dendritic Cell Type 1",
  "21-1" = "Dendritic Cell Type 2",
  "21-0" = "Dendritic Cell Type 3",
  "26"   = "Plasmacytoid Dendritic Cell",
  
  "0"    = "SOD2hi CD14+ Mono",
  "1"    = "IL-1Bhi CD14+ Mono",
  "10"   = "MHChi CD14+ CD16- Mono",
  "14"   = "CD14+ CD16- Mono",
  "25"   = "ISG+ CD14+ Mono",
  "12"   = "C1Qhi CD16+ Mono",
  "24-1" = "CD16+ Mono",
  
  "11"   = "Activated Naive B",
  "15"   = "Switched Memory B",
  "19"   = "CD5+ B",
  "24-0" = "Atypical Memory B",
  "29"   = "Plasmablasts",
  
  "2"    = "CD56dim NK",
  "18"   = "CD56bright NK",
  "27"   = "Platelets",
  "30"   = "Hematopoietic Stem and Progenitors",
  "32"   = "Neutrophils",
  
  "3"    = "CD26hi CD4+ Naïve T",
  "4"    = "CD8+ Naïve T",
  "5"    = "CD4+ Naïve T",
  "6"    = "CD4+ TH2 Effector Stem Memory",
  "7"    = "CD4+ TFH Stem Memory",
  "8"    = "CD8+ TC1 Effector Memory",
  "9"    = "CD8+ NKTC1 Terminal Effector",
  "13-0" = "CD4+ TH17 Effector Memory",
  "13-1" = "CD8+ MAIT-like TC17 Effector Memory",
  "16"   = "CD4+ T Effector Memory",
  "17"   = "CD161+ CD8+ NKTC1 Terminal Effector",
  "20"   = "CD4+ Treg Central Memory",
  "22"   = "CD8+ NKTC1 Effector Progenitor",
  "31"   = "CD4- CD8dim Treg Stem Memory",
  "23"   = "CD4+ TFH1 Central Memory",
  "28"   = "Ki67+ CD8+ NKTC1 Terminal Effector"
)

DefaultAssay(s) <- "RNA"


## Choose Cluster Order and Set ----

#All clusters
clusters_to_plot <- c("3", "5", "4", "16", "6", "7", "31", "20", "13-0", "23", "8", "9", "17","22", "13-1","28", "2", "18", "12","24-1", "25", "0", "1", "10", "14", "27", "21-0", "21-1", "21-2", "26", "11", "15", "19", "24-0", "29", "30", "32") #all clusters main markers

#T cells
clusters_to_plot <- c("3","4","5","6","7","8","9","13-0","13-1","16","17","20","22","31") #T cell clusters in increasing seurat cluster number order

#DCs
clusters_to_plot <- c("21-2", "21-1", "21-0", "26") 

#B cells + plasmablasts
clusters_to_plot <- c("11", "15", "19", "24-0", "29") 

#B cells without plasmablasts
clusters_to_plot <- c("11", "15", "19", "24-0") 

#Monocytes 
clusters_to_plot <- c("14", "0", "10", "1", "25", "24-1", "12") 

#NKs
clusters_to_plot <- c("2", "18") 


## Choose Marker Sets ----

#GENERAL MARKERS FOR ALL CLUSTERS
chosen_markers <- c("CD3D", "CD3G", "CD3E", "ZAP70", "NCAM1", "NCR1", "FCGR3A", "SYK", "CD14", "PF4", "CLEC10A", "CLEC9A", "LILRA4", "PAX5", "MS4A1", "CD19", "IGHM", "IGHD", "JCHAIN", "CD34", "LTF", "PAX5")

#DC MARKERS
chosen_markers <- c("CLEC9A", "XCR1", "CD81", "LAMP3", "ITGA1", "ITGB7", "ITGA4", "CD5", "CCR7", "CD14", "ITGAX", "CX3CR1", "FCER1A", "CLEC10A", "CD163", "LILRA4", "CLEC4C", "JCHAIN", "TCF4", "AXL", "SIGLEC6")

#B LINEAGE CELL MARKERS
#B CELLS MARKERS
chosen_markers <- c( "MS4A1", "TCL1A", "YBX3", "FCER2", "IL4R", "FCER2", "IL4R", "CD38", "CXCR5", "CD69", "CD83", "IGHM", "IGHD", "CD19", "CR2", "CD24", "CD27", "AIM2", "TNFRSF13B", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "FAS", "CD27", "CD5", "LEF1", "TBX21", "SPN", "ITGAX", "ZEB2", "MEF2B", "CD86", "BCL6", "LMO2", "FOS", "JUN", "JUNB")
#B CELL WITH PLASMABLAST MARKERS
chosen_markers <- c( "MS4A1", "TCL1A", "YBX3", "FCER2", "IL4R", "FCER2", "IL4R", "CXCR5", "CD69", "CD83", "IGHM", "IGHD", "CD19", "CR2", "CD24", "CD27", "AIM2", "TNFRSF13B", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "FAS", "CD27", "CD5", "LEF1", "TBX21", "SPN", "ITGAX", "ZEB2", "MEF2B", "CD86", "BCL6", "LMO2", "FOS", "JUN", "JUNB","TXNDC5", "MYDGF", "JCHAIN", "MZB1", "XBP1", "PRDM1", "CD38", "MKI67", "SDC1")

#MONOCYTE MARKERS
chosen_markers <- c("CD14", "S100A8", "S100A9", "S100A12", "FCGR3A","SOD2", "MGST1", "NAMPT", "CST3", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1","IL1B", "IRF8", "ICAM1", "JUNB", "EGR3", "CCL3", "CCR2","ISG15", "IFI44L", "MX1", "MX2", "MS4A7", "CX3CR1", "C1QA", "C1QB", "C1QC", "CD68")

#T CELL MARKERS
chosen_markers <- c(
  #LINEAGE (kind of, outside of gamma-delts)
  "CD4", "CD8A", "CD8B", "KLRB1", "ZBTB16", "SLC4A10", 
  
  #ACTIVATION STATUS 
  "IL6ST", "FOXO1", "TCF7", "CD27", "IL7R", "CD28", "LEF1", "FAS", "CXCR4", "CD44", "IL2RA", "IL2RB", "KLF2", "DPP4",
  
  #MIGRATION
  "CCR7", "SELL", "S1PR1", "CD69", "ITGA1","ITGA4","ITGAE","ITGB1","ITGB7", "ABCB1",
  
  #FUNCTION/POLARIZATION
  "CXCR3", "GATA3", "CCR4", "PTGDR2", "BCL6", "MAF", "ICOS", "CXCR5", "IKZF2", "FOXP3", "CCR8", "CTLA4", "TIGIT", "HAVCR2", "TNFRSF1B", "CCR10", "AHR", "GPR35", "CCR6", "RORC", "IL23R", "KIT", "IL10RA", "IL18R1", "IL18RAP", "NCAM1", "B3GAT1", "CD160", "EOMES", "TBX21", "IFNG", "CCR5", "CX3CR1", "NCAM1", "FCGR3A", "FCER1G", "KLRK1", "GZMK", "KLRC1", "KLRC2", "PRF1", "KLRG1", "GZMB", "GNLY", "FGFBP2", "ZEB2", "PRDM1", "LAG3", "PDCD1", "FASLG", "FGL2", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2", "CD38", "MKI67"
) 

#NK MARKERS 
chosen_markers <- c("NCAM1", "FCGR3A", "KLRC1", "KLRC2", "FCER1G", "GZMB", "GZMK", "B3GAT1", "MX1", "MKI56")

#PLOTTING THE MANUALLY CHOSEN CLUSTERS WITH MANUALLY CHOSEN MARKER SETS
DefaultAssay(s) <- "RNA"

# Preserve manually specified cluster and gene order
clusters <- as.character(clusters_to_plot)
markers <- unique(unlist(chosen_markers, use.names = FALSE))


## Plot ----

# Create dot
p <- DotPlot(
  s,
  features = markers,
  idents = clusters,
  group.by = "seurat_clusters"
) + RotatedAxis()

# Set cluster order (bottom to top)
p$data$id <- factor(p$data$id, levels = clusters)

# Set gene order (left to right)
p$data$features.plot <- factor(
  p$data$features.plot,
  levels = markers
)
print(p)
