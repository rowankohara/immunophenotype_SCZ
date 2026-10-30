# Rowan O'Hara-Payne
# 06/17/2026
# Cell-cell communication analysis and differential signaling using dominoSignal (PMID: 41746282)
# used the available vignettes: https://fertiglab.github.io/dominoSignal/articles/dominoSignal.html
####################################################################################################################


# Downsample Cells for PySCENIC ----
# running pySCENIC is the first step of this workflow
# I needed to downsample my cells to get it to complete a run

library(Seurat)
library(tidyverse)

obj <- readRDS("kelly_postqc_clusters_meta.rds")

celltype_col <- "seurat_clusters"
condition_col <- "Class"

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

# run pySCENIC using the bash script (08.1-pySCENIC.sh)
# based on https://fertiglab.github.io/dominoSignal/articles/tf_scenic_vignette.html


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
library(Seurat)

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
levels(pbmc$seurat_clusters) <- celltype_labels[levels(pbmc$seurat_clusters)]
clusters <- factor(pbmc$seurat_clusters)
names(clusters) <- colnames(pbmc)

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

# save domino object
saveRDS(pbmc_dom, "kelly_domino_obj.RDS")


# plotting results
celltype_pal <- c(
  "CD26hiCD4TSbN"    = "coral3",
  "CD4TSbN"    = "orangered",
  "CD8TSbN"    = "chocolate1",
  "CD4TDwM"   = "goldenrod",
  "CD4Th2DwMp"    = "darkgoldenrod",
  "CD4TfhSbMp"    = "sandybrown",
  "CD4negCD8dimTregSM"   = "darkorange3",
  "CD4TregSwM"   = "indianred1",
  "CD4Th17DwM" = "red2",
  "CD4Tfh1SwM"   = "orange",
  "CD8Tc1DwM"    = "rosybrown1",
  "CD8NKTc1DwAt"    = "tomato",
  "CD161CD8NKTc1DwAt"   = "tan1",
  "CD8NKTc1DwAp"   = "salmon",
  "CD8MAITlikeTc1DwM" = "palevioletred2",
  "Ki67CD8NKTc1DwAt"   = "deeppink3",
  "CD56dimNK"    = "gold",
  "CD56brightNK"   = "yellow",
  "C1QhiCD16Mono"   = "dodgerblue2",
  "CD16Mono" = "cyan",
  "ISGCD14Mono"   = "aquamarine3",
  "SOD2hiCD14Mono"    = "darkturquoise",
  "IL1BhiCD14Mono"    = "steelblue2",
  "MHChiCD14Mono"   = "cadetblue2",
  "CD14Mono"   = "darkcyan",
  "Platelets"   = "orchid3",
  "DC3" = "dodgerblue3",
  "DC2" = "deepskyblue4",
  "DC1" = "darkblue",
  "pDC"   = "royalblue3",
  "ActivatedNaiveB"   = "springgreen2",
  "SwitchedMemoryB"   = "limegreen",
  "CD5B"   = "olivedrab",
  "AtypicalMemoryB" = "seagreen",
  "Plasmablasts"   = "mediumseagreen",
  "HSPCs"   = "navy",
  "Neutrophils"   = "mediumpurple3"
)

celltype_labels <- c(
  "3"    = "CD26hiCD4TSbN",
  "5"    = "CD4TSbN",
  "4"    = "CD8TSbN",
  "16"   = "CD4TDwM",
  "6"    = "CD4Th2DwMp",
  "7"    = "CD4TfhSbMp",
  "31"   = "CD4negCD8dimTregSM",
  "20"   = "CD4TregSwM",
  "13-0" = "CD4Th17DwM",
  "23"   = "CD4Tfh1SwM",
  "8"    = "CD8Tc1DwM",
  "9"    = "CD8NKTc1DwAt",
  "17"   = "CD161CD8NKTc1DwAt",
  "22"   = "CD8NKTc1DwAp",
  "13-1" = "CD8MAITlikeTc1DwM",
  "28"   = "Ki67CD8NKTc1DwAt",
  "2"    = "CD56dimNK",
  "18"   = "CD56brightNK",
  "12"   = "C1QhiCD16Mono",
  "24-1" = "CD16Mono",
  "25"   = "ISGCD14Mono",
  "0"    = "SOD2hiCD14Mono",
  "1"    = "IL1BhiCD14Mono",
  "10"   = "MHChiCD14Mono",
  "14"   = "CD14Mono",
  "27"   = "Platelets",
  "21-0" = "DC3",
  "21-1" = "DC2",
  "21-2" = "DC1",
  "26"   = "pDC",
  "11"   = "ActivatedNaiveB",
  "15"   = "SwitchedMemoryB",
  "19"   = "CD5B",
  "24-0" = "AtypicalMemoryB",
  "29"   = "Plasmablasts",
  "30"   = "HSPCs",
  "32"   = "Neutrophils"
)

pdf("plots/domino_feature_heatmap.pdf", width = 10, height = 10)
feat_heatmap(pbmc_dom, norm = TRUE, bool = FALSE, use_raster = FALSE, cols = celltype_pal, row_names_gp = grid::gpar(fontsize = 4))
dev.off()

pdf("plots/domino_incoming_signaling_heatmap.pdf", width = 10, height = 12)
for (type in celltypes) {
  incoming_signaling_heatmap(pbmc_dom, rec_clust = type, max_thresh = 2.5, use_raster = FALSE)
}
dev.off()

pdf("plots/domino_corr_heatmap.pdf", width = 14, height = 12)
cor_heatmap(pbmc_dom, title = "PBMC R-TF Correlations", column_names_gp = grid::gpar(fontsize = 8))
cor_heatmap(pbmc_dom, bool = TRUE, bool_thresh = 0.25)
cor_heatmap(pbmc_dom, bool = FALSE, mark_connections = TRUE)
dev.off()

pdf("plots/domino_signaling_heatmap.pdf", width = 10, height = 10)
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

# need to run per sample
domino_list <- lapply(samples, function(s) {
  sub <- subset(obj, sample == s)
  pbmc <- as.SingleCellExperiment(sub)
  counts <- SummarizedExperiment::assay(pbmc, "counts")
  logcounts <- SummarizedExperiment::assay(pbmc, "logcounts")
  logcounts <- logcounts[rowSums(logcounts) > 0, ]
  z_scores <- t(scale(t(logcounts)))
  clusters <- factor(pbmc$seurat_clusters)
  names(clusters) = colnames(pbmc)

  pbmc_dom <- create_domino(rl_map = rl_map, 
                            features = auc_in[, colnames(pbmc)], 
                            counts = counts, 
                            z_scores = z_scores,
                            clusters = clusters, 
                            tf_targets = regulon_list, 
                            use_clusters = TRUE, 
                            use_complexes = TRUE,
                            remove_rec_dropout = FALSE)
  pbmc_dom <- build_domino(dom = pbmc_dom, 
                           min_tf_pval = 0.001, 
                           max_tf_per_clust = 25,
                           max_rec_per_tf = 25, 
                           rec_tf_cor_threshold = 0.25, 
                           min_rec_percentage = 0.1)
  pbmc_dom
})

names(domino_list) <- samples

dom_summary <- summarize_linkages(domino_list, meta, subject_names = NULL)


# Run differential signaling test for all celltypes and combine results into a dataframe
celltypes <- unique(obj$seurat_clusters)
results_list <- list()
for (type in celltypes) {

  cat("Testing:", type, "\n")

  res <- tryCatch(
    test_differential_linkages(
      dom_summary,
      type,
      "Class",
      linkage = "rec_lig",
      subject_names = samples,
      test_name = "fishers.exact"
    ),
    error = function(e) {
      message(type, " failed: ", e$message)
      return(NULL)
    }
  )

  if (!is.null(res)) {
    res$cluster <- type
    results_list[[type]] <- res
  }
}

diff_link <- do.call(rbind, results_list)
rownames(diff_link) <- NULL


# Summary plots -----

# recode cluster labels to their spelled out form
diff_link$cluster <- celltype_labels[diff_link$cluster]

# a function to get the top 3 senders for a particular ligand
get_top_n_senders <- function(dom, receiver, receptor, ligand, n = 3) {
  
  if (!receiver %in% names(dom@cl_signaling_matrices)) {
    return(NULL)
  }
  
  mat <- dom@cl_signaling_matrices[[receiver]]
  
  if (!ligand %in% rownames(mat)) {
    return(NULL)
  }
  
  scores <- mat[ligand, ]
  
  # remove zero scores
  scores <- scores[scores > 0]
  
  # sort and take top n
  scores <- sort(scores, decreasing = TRUE)
  scores <- head(scores, n)
  
  if (length(scores) == 0) {
    return(NULL)
  }
  
  data.frame(
    receiver = receiver,
    receptor = receptor,
    ligand = ligand,
    sender = sub("^L_", "", names(scores)),
    sender_score = as.numeric(scores),
    sender_rank = seq_along(scores),
    stringsAsFactors = FALSE
  )
}

sig_diff <- diff_link[diff_link$p.value < 0.05, ]
sig_diff$receptor <- sub(" <- .*", "", sig_diff$feature)
sig_diff$ligand <- sub(".* <- ", "", sig_diff$feature)

### Supplemental Data Table: differential_signaling ----
write.csv(sig_diff, "sig_diff_linkages_all.csv")


top_senders_list <- vector("list", nrow(sig_diff))

# iterate through the rows of significantly different signaling linkages
for (i in seq_len(nrow(sig_diff))) {

  receiver <- sig_diff$cluster[i]
  receptor <- sig_diff$receptor[i]
  ligand <- sig_diff$ligand[i]

  senders <- get_top_n_senders(
    dom = pbmc_dom,
    receiver = receiver,
    receptor = receptor,
    ligand = ligand,
    n = 3
  )

  if (!is.null(senders)) {

    # add differential signaling statistics
    result_info <- sig_diff[rep(i, nrow(senders)), ]

    top_senders_list[[i]] <- cbind(
      result_info,
      senders[, c("sender", "sender_score", "sender_rank")]
    )
  }
}

top_senders <- do.call(rbind, top_senders_list)
rownames(top_senders) <- NULL

### Supplemental Data Table: incoming_ligand_senders ----
write.csv(top_senders, "incoming_ligand_senders.csv")
