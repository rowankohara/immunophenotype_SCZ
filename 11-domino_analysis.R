# Rowan O'Hara-Payne
# 06/17/2026
# Script to prep data for dominoSignal


# Downsample Cells for PySCENIC ----

library(Seurat)
library(tidyverse)

obj <- readRDS("kelly_postqc_clusters_meta.rds")

celltype_col <- "seurat_clusters"
condition_col <- "Group"

set.seed(123)

max_cells <- 150  # maximum number of cells for each cluster x group strata

meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell")

meta$strata <- paste(
  meta[[celltype_col]],
  meta[[condition_col]],
  sep = "_"
)

cell_counts <- meta %>%
  count(strata)

# remove strata with too few cells
keep_strata <- cell_counts %>%
  filter(n >= 30) %>%
  pull(strata)
meta_filt <- meta %>%
  filter(strata %in% keep_strata)

# downsample using max-cells
cells_keep <- meta_filt %>%
  group_by(strata) %>%
  group_modify(~ {
    
    n_keep <- min(nrow(.x), max_cells)
    
    slice_sample(.x, n = n_keep)
    
  }) %>%
  pull(cell)
obj_down <- subset(obj, cells = cells_keep)

# reduce the number of genes by selecting for TF genes and the most variable features
counts <- GetAssayData(obj_down, assay = "RNA", layer = "counts")
keep_genes <- rowSums(counts > 0) >= 100  # remove lowly expressed genes
obj_down <- subset(obj_down, features = rownames(counts)[keep_genes])
obj_down <- FindVariableFeatures(
  obj_down,
  nfeatures = 5000
)
tf_list <- readLines("scenic/allTFs_hg38.txt")
genes_use <- unique(c(
  VariableFeatures(obj_down),
  tf_list
))
# remove mitochondrial genes and immunoglobin genes
genes_use <- genes_use[
  !grepl("^MT-", genes_use) &
  !grepl("^RPL|^RPS", genes_use) &
  !grepl("^IG[HKL]", genes_use)
]
obj_down <- subset(obj_down, features = genes_use)

table(
  obj_down@meta.data[[celltype_col]],
  obj_down@meta.data[[condition_col]]
)

ncol(obj_down)

# save downsampled Seurat object as matrix for pySCENIC 
pbmc <- as.SingleCellExperiment(obj_down)
pbmc_counts <- SummarizedExperiment::assay(pbmc, "counts")
write.table(t(as.matrix(pbmc_counts)), "scenic/pbmc13k_counts.tsv", sep = "\t",
    col.names = NA)

# run pySCENIC using the bash script


# Download CellPhone DB ----
# URL for desired version of CellPhoneDB
cellphone_url <- "https://github.com/ventolab/cellphonedb-data/archive/refs/tags/v4.0.0.tar.gz"

# download compressed database
cellphone_tar <- "~/cellphoneDB_v4.tar.gz"
download.file(url = cellphone_url, destfile = cellphone_tar)

# move contents of the compressed file to a new directory
untar(tarfile = cellphone_tar, exdir = "cellphone")
cellphone_data <- "/local/projects/idea/rohara/immunophenotype_SCZ/cellphone/cellphonedb-data-4.0.0/data"



# Load all input files ----
library(dominoSignal)

regulons <- read.csv("scenic/regulons_pbmc_13k.csv")
auc <- read.table("scenic/auc_pbmc_60k.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE,
    sep = ",")

regulons <- regulons[-1:-2, ]
colnames(regulons) <- c("TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity",
    "Annotation", "Context", "TargetGenes", "RankAtMax")
regulon_list <- create_regulon_list_scenic(regulons = regulons)

auc_in <- as.data.frame(t(auc))
# Remove pattern '...' from the end of all rownames:
rownames(auc_in) <- gsub("\\.\\.\\.$", "", rownames(auc_in))

complexes <- read.csv(paste0(cellphone_data, "/complex_input.csv"), stringsAsFactors = FALSE)
genes <- read.csv(paste0(cellphone_data, "/gene_input.csv"), stringsAsFactors = FALSE)
interactions <- read.csv(paste0(cellphone_data, "/interaction_input.csv"), stringsAsFactors = FALSE)
proteins <- read.csv(paste0(cellphone_data, "/protein_input.csv"), stringsAsFactors = FALSE)

rl_map <- create_rl_map_cellphonedb(genes = genes, proteins = proteins, interactions = interactions,
    complexes = complexes, database_name = "CellPhoneDB_v4.0"  # database version used
)

knitr::kable(head(rl_map))


# Start domino analysis ----
obj <- readRDS("kelly_postqc_clusters_meta.rds")
pbmc <- as.SingleCellExperiment(obj)
counts <- SummarizedExperiment::assay(pbmc, "counts")
logcounts <- SummarizedExperiment::assay(pbmc, "logcounts")
logcounts <- logcounts[rowSums(logcounts) > 0, ]
z_scores <- t(scale(t(logcounts)))
clusters <- factor(pbmc$seurat_clusters)
names(clusters) = colnames(pbmc)

pbmc_dom <- create_domino(rl_map = rl_map, features = auc_in, counts = counts, z_scores = z_scores,
    clusters = clusters, tf_targets = regulon_list, use_clusters = TRUE, use_complexes = TRUE,
    remove_rec_dropout = FALSE)
# rl_map: receptor - ligand map data frame
# features: TF activation scores (AUC matrix) 
# counts: counts matrix 
# z_scores: scaled expression data 
# clusters: named vector of cell cluster assignments 
# tf_targets: list of TFs and their regulons 
# use_clusters: assess receptor activation and ligand expression on a per-cluster basis 
# use_complexes: include receptors and genes that function as a complex in results 
# remove_rec_dropout: whether to remove zeroes from correlation calculations

pbmc_dom <- build_domino(dom = pbmc_dom, min_tf_pval = 0.001, max_tf_per_clust = 25,
    max_rec_per_tf = 25, rec_tf_cor_threshold = 0.25, min_rec_percentage = 0.1)
# min_tf_pval: Threshold for p-value of DE for TFs 
# rec_tf_cor_threshold: Minimum correlation between receptor and TF 
# min_rec_percentage: Minimum percent of cells that must express receptor


# plotting results
celltype_pal <- c(
  "3"    = "coral3",
  "5"    = "orangered",
  "4"    = "chocolate1",
  "16"   = "goldenrod",
  "6"    = "darkgoldenrod",
  "7"    = "sandybrown",
  "31"   = "darkorange3",
  "20"   = "indianred1",
  "13-0" = "red2",
  "23"   = "orange",
  "8"    = "rosybrown1",
  "9"    = "tomato",
  "17"   = "tan1",
  "22"   = "salmon",
  "13-1" = "palevioletred2",
  "28"   = "deeppink3",
  "2"    = "gold",
  "18"   = "yellow",
  "12"   = "dodgerblue2",
  "24-1" = "cyan",
  "25"   = "aquamarine3",
  "0"    = "darkturquoise",
  "1"    = "steelblue2",
  "10"   = "cadetblue2",
  "14"   = "darkcyan",
  "27"   = "orchid3",
  "21-0" = "dodgerblue3",
  "21-1" = "deepskyblue4",
  "21-2" = "darkblue",
  "26"   = "royalblue3",
  "11"   = "springgreen2",
  "15"   = "limegreen",
  "19"   = "olivedrab",
  "24-0" = "seagreen",
  "29"   = "mediumseagreen",
  "30"   = "navy",
  "32"   = "mediumpurple3"
)
celltype_labels <- c(
  "3"    = "CD26^'hi'*CD4^'+'*TS['B']*N",
  "5"    = "CD4^'+'*TS['B']*N",
  "4"    = "CD8^'+'*TS['B']*N",
  "16"   = "CD4^'+'*TD['W']*M",
  "6"    = "CD4^'+'*T['H2']*D['W']*M['P']",
  "7"    = "CD4^'+'*T['FH']*S['B']*M['P']",
  "31"   = "CD4^'-'*CD8^'dim'*T['Reg']*SM",
  "20"   = "CD4^'+'*T['Reg']*S['W']*M",
  "13-0" = "CD4^'+'*T['H17']*D['W']*M",
  "23"   = "CD4^'+'*T['FH1']*S['W']*M",
  "8"    = "CD8^'+'*T['C1']*D['W']*M",
  "9"    = "CD8^'+'*NKT['C1']*D['W']*A['T']",
  "17"   = "CD161^'+'*CD8^'+'*NKT['C1']*D['W']*A['T']",
  "22"   = "CD8^'+'*NKT['C1']*D['W']*A['P']",
  "13-1" = "CD8^'+'*'MAIT-like'*T['C1']*D['W']*M",
  "28"   = "Ki67^'+'*CD8^'+'*NKT['C1']*D['W']*A['T']",
  "2"    = "CD56^'dim'~NK",
  "18"   = "CD56^'bright'~NK",
  "12"   = "C1Q^'hi'*CD16^'+'~Mono",
  "24-1" = "CD16^'+'~Mono",
  "25"   = "ISG^'+'*CD14^'+'~Mono",
  "0"    = "SOD2^'hi'*CD14^'+'~Mono",
  "1"    = "'IL-1B'^'hi'*CD14^'+'~Mono",
  "10"   = "MHC^'hi'*CD14^'+'~Mono",
  "14"   = "CD14^'+'~Mono",
  "27"   = "Platelets",
  "21-0" = "DC3",
  "21-1" = "DC2",
  "21-2" = "DC1",
  "26"   = "pDC",
  "11"   = "Activated~Naive~B",
  "15"   = "Switched~Memory~B",
  "19"   = "CD5^'+'~B",
  "24-0" = "Atypical~Memory~B",
  "29"   = "Plasmablasts",
  "30"   = "HSPCs",
  "32"   = "Neutrophils"
)

pdf("plots/domino_feature_heatmap.pdf")
feat_heatmap(pbmc_dom, norm = TRUE, bool = FALSE, use_raster = FALSE, cols = celltype_pal, row_names_gp = grid::gpar(fontsize = 4))
dev.off()

pdf("plots/domino_signaling_diagram.pdf")
signaling_network(pbmc_dom, cols = celltype_pal, edge_weight = 0.5, max_thresh = 3)
signaling_network(pbmc_dom, cols = celltype_pal, showOutgoingSignalingClusts = c("13-0"), edge_weight = 0.5, max_thresh = 3)
signaling_network(pbmc_dom, cols = celltype_pal, showIncomingSignalingClusts = c("13-0"), edge_weight = 0.5, max_thresh = 3)
dev.off()

pdf("plots/domino_gene_network_th17.pdf", width = 10, height = 10)
gene_network(pbmc_dom, clust = "13-0", layout = "sphere")
gene_network(pbmc_dom, clust = "13-0", OutgoingSignalingClust = "21-1", layout = "sphere")
gene_network(pbmc_dom, clust = "13-0", OutgoingSignalingClust = "26", layout = "sphere")
gene_network(pbmc_dom, clust = "13-0", OutgoingSignalingClust = "19", layout = "sphere")
dev.off()

pdf("plots/domino_gene_network_from_th17.pdf", width = 10, height = 10)
gene_network(pbmc_dom, clust = "24-0", OutgoingSignalingClust = "13-0", layout = "grid")
gene_network(pbmc_dom, clust = "21-1", OutgoingSignalingClust = "13-0", layout = "grid")
gene_network(pbmc_dom, clust = "26", OutgoingSignalingClust = "13-0", layout = "grid")
gene_network(pbmc_dom, clust = "1", OutgoingSignalingClust = "13-0", layout = "grid")
gene_network(pbmc_dom, clust = "0", OutgoingSignalingClust = "13-0", layout = "grid")
gene_network(pbmc_dom, clust = "24-1", OutgoingSignalingClust = "13-0", layout = "grid")
dev.off()

pdf("plots/domino_incoming_signaling_heatmap.pdf", width = 10, height = 12)
incoming_signaling_heatmap(pbmc_dom, rec_clust = "13-0", max_thresh = 2.5, use_raster = FALSE)
dev.off()

pdf("plots/domino_circos_CD74.pdf")
circos_ligand_receptor(pbmc_dom, receptor = "CD74")
dev.off()

pdf("plots/domino_corr_heatmap.pdf", width = 14, height = 12)
cor_heatmap(pbmc_dom, title = "PBMC R-TF Correlations", column_names_gp = grid::gpar(fontsize = 8))
cor_heatmap(pbmc_dom, bool = TRUE, bool_thresh = 0.25)
cor_heatmap(pbmc_dom, bool = FALSE, mark_connections = TRUE)
dev.off()

pdf("plots/domino_signaling_heatmap.pdf")
signaling_heatmap(pbmc_dom)
signaling_heatmap(pbmc_dom, scale = "sqrt")
signaling_heatmap(pbmc_dom, normalize = "rec_norm")
dev.off()


# Differential signaling testing ----
obj <- readRDS('kelly_postqc_clusters_meta.rds')
meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
for (i in 1:15) {
        if (startsWith(unique(obj@meta.data$sample)[i], meta$MPRCID[i])) {
                meta$sample[i] <- unique(obj@meta.data$sample)[i]
}}
meta <- meta[, c(17, (1:ncol(meta))[-17])]
samples <- meta$sample

domino_list <- lapply(samples, function(s) {

    sub <- subset(obj, sample == s)
    pbmc <- as.SingleCellExperiment(sub)
    counts <- SummarizedExperiment::assay(pbmc, "counts")
    logcounts <- SummarizedExperiment::assay(pbmc, "logcounts")
    logcounts <- logcounts[rowSums(logcounts) > 0, ]
    z_scores <- t(scale(t(logcounts)))
    clusters <- factor(pbmc$seurat_clusters)
    names(clusters) = colnames(pbmc)

    create_domino(rl_map = rl_map, 
		  features = auc_in, 
		  counts = counts, 
		  z_scores = z_scores,
		  clusters = clusters, 
		  tf_targets = regulon_list, 
		  use_clusters = TRUE, 
		  use_complexes = TRUE,
		  remove_rec_dropout = FALSE)
})

names(domino_list) <- samples


dom_summary <- summarize_linkages(domino_list, meta, subject_names = NULL)


th17_signaling_diff <- test_differential_linkages(
  dom_summary,
  "13-0",
  "Group",
  linkage = "rec_lig",
  subject_names = samples,
  test_name = "fishers.exact"
)


