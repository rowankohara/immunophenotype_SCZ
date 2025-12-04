# Rowan K. O'Hara-Payne
# 5/25/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Do differential abundance analysis


# Read in libraries and data ----
library(Seurat)
library(dplyr)
library(edgeR)
library(ggplot2)
library(pheatmap)

# read in Seurat object
obj <- readRDS('dkelly_singlets.postqc.azimuth_celltypes.rds')

# read in pseudobulk data
elist <- readRDS('dkelly_singlets.pseudobulk.rds')

# read in meta data
sample <- colnames(elist[[11]])
MPRCID <- gsub( '\\_(.*)' , '' , sample )
MPRCID <- gsub( '\\-(.*)' , '' , MPRCID )

meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
rownames(meta) <- meta$MPRCID
meta <- meta[ MPRCID , ]


# Differential Abundance Analysis ----

# create an empty data frame to store cell-level metadata
cell_meta <- data.frame(row.names = colnames(obj))

# loop through each metadata column (excluding MPRCID)
meta_filtered <- meta %>% select(MPRCID, Sex, Group, Age)
for (col in colnames(meta_filtered[-1])) {
  vec <- rep(NA, ncol(obj))
  names(vec) <- colnames(obj)
  
  for (i in 1:nrow(meta_filtered)) {
    sample_id <- meta_filtered$MPRCID[i]
    
    # find cells whose sample value starts with the sample_id
    matched_idx <- grep(paste0("^", sample_id), obj@meta.data$sample)
    matched_cells <- colnames(obj)[matched_idx]
    
    # assign sample-level metadata to those cells
    vec[matched_cells] <- meta_filtered[[col]][i]
  }
  
  cell_meta[[col]] <- vec
}

# add the metadata to the Seurat object
obj <- Seurat::AddMetaData(obj, metadata = cell_meta)

# summarize to one row per sample
sample_metadata <- obj@meta.data %>%
  select(sample, Group, Sex, Age) %>%
  distinct()

# cell counts per sample per cell type
counts_table <- table(obj$sample, obj$predicted.celltype.l2)
counts_matrix <- as.matrix(counts_table)  # samples x cell types

# match sample order
sample_metadata <- sample_metadata %>%
  filter(sample %in% rownames(counts_matrix)) %>%
  arrange(match(sample, rownames(counts_matrix)))

# confirm alignment
stopifnot(all(sample_metadata$sample == rownames(counts_matrix)))

# transpose so rows = cell types, cols = samples
dge <- DGEList(counts = t(counts_matrix))
dge$sample <- cbind(dge$sample, sample_metadata)
dge <- calcNormFactors(dge)

# convert categorical variables to factors
sample_metadata$Group <- factor(sample_metadata$Group)
sample_metadata$Group <- forcats::fct_recode(sample_metadata$Group, "SCZLow" = "SCZ Low", "SCZHigh" = "SCZ High")
sample_metadata$Sex <- factor(sample_metadata$Sex)

# create design matrix
design <- model.matrix(~ Group + Age + Sex, data = sample_metadata)
colnames(design)

# model fitting and testing
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)

# test the effect of group (e.g., SCZ High vs SCZ Low vs Control)
lrt_low <- glmLRT(fit, coef = "GroupSCZLow")
res_low <- topTags(lrt_low, n = Inf)$table
res_low$comparison <- "SCZ Low vs HC"

lrt_high <- glmLRT(fit, coef = "GroupSCZHigh")
res_high <- topTags(lrt_high, n = Inf)$table
res_high$comparison <- "SCZ High vs HC"

res_all <- bind_rows(res_low, res_high)
res_all$cell_type <- rownames(res_all)

# make contrast for SCZ HIgh vs SCZ Low
contrast_matrix <- limma::makeContrasts(SCZHigh_vs_SCZLow = GroupSCZHigh - GroupSCZLow, levels = design)

# fit model
fit2 <- glmFit(dge, design)
lrt_contrast <- glmLRT(fit2, contrast = contrast_matrix)
res_contrast <- topTags(lrt_contrast, n = Inf)$table
res_contrast$comparison <- "SCZ High vs SCZ Low"
res_contrast$cell_type <- rownames(res_contrast)

# combine results
res_all <- bind_rows(res_all, res_contrast)

# fix celltype names
res_all$cell_type <- stringr::str_remove(res_all$cell_type, "\\.\\.\\.[0-9]+")
write.csv(res_all, "DAseq_edgeR_results_allcomparisons.csv", row.names=F)


# Plots ----

## volcano plot ----

# classify significance and direction
volcano_df <- res_all %>%
  mutate(
    direction = case_when(
      logFC > 0 ~ "Higher in Group 1",
      logFC < 0 ~ "Lower in Group 1",
      TRUE ~ "No Change"
    )
  )

pdf("plots/DAseq_volcano.pdf")
ggplot(volcano_df, aes(x = logFC, y = -log10(FDR), color = direction)) +
  geom_jitter(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("Higher in Group 1" = "firebrick",
                                "Lower in Group 1" = "dodgerblue",
                                "No Change" = "gray80")) +
  #geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  facet_wrap(~comparison) +
  theme_bw(base_size = 12) +
  labs(
    title = "Volcano Plots of Cell Type DA Across Comparisons",
    x = "log2 Fold Change",
    y = "-log10(FDR)",
    color = "Direction"
  )
dev.off()


## heatmap ----

# pivot to wide format to have a matrix of logFCs; row = cell type, col = contrast
lfc_matrix <- res_all %>%
  select(cell_type, logFC, comparison) %>%
  tidyr::pivot_wider(names_from = comparison, values_from = logFC) %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()

pdf("plots/DAseq_heatmap.pdf")
pheatmap(lfc_matrix,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "logFC of Cell Types")
dev.off()

# none at FDR < 0.05
sig_celltypes <- res_all %>%
  filter(FDR < 0.05) %>%
  pull(cell_type) %>%
  unique()

lfc_matrix_sig <- lfc_matrix[rownames(lfc_matrix) %in% sig_celltypes, ]


## boxplots ----

# convert to long format and compute proportions
prop_df <- counts_table %>%
  as.data.frame() %>%
  rename(sample = Var1, cell_type = Var2, count = Freq) %>%
  group_by(sample) %>%
  mutate(total = sum(count),
         prop = count / total) %>%
  ungroup()

# add sample meta data
prop_df <- left_join(prop_df, sample_metadata, by = "sample")

pdf("plots/DAseq_boxplots.pdf")
ggplot(prop_df, aes(x = Group, y = prop, fill = Group)) +
  geom_boxplot() +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_bw() +
  labs(y = "Proportion", title = "Cell Type Proportions by Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("plots/DAseq_Treg_boxplot.pdf")
prop_treg <- prop_df %>% filter(cell_type == "Treg")
ggplot(prop_treg, aes(x = Group, y = prop, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Proportion", title = "Treg Proportions by Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("plots/celltype_prop_DAseq.pdf")
ggplot(prop_df, aes(x = cell_type, y = prop, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Cell Type Proportion", x = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


## bar plot ----

# bar plot of cell type proportions for each sample
pdf("plots/celltype_sample_DAseq.pdf")
ggplot(prop_df, aes(x = sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Sample", y = "Proportion of Cells", fill = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x")
dev.off()



# Milo Strategy ----
# using miloR package for differential abundance testing using a K-nearest neighbors algorithm
# using this vignette for contrasts https://bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_contrasts.html
library(miloR)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)

obj <- readRDS("kelly_postqc_celltypes_meta.rds")
obj_sce <- as.SingleCellExperiment(obj)
obj_milo <- Milo(obj_sce)

# get the UMAP from the single cell experiment object
reducedDim(obj_milo, "UMAP") <- reducedDim(obj_sce, "UMAP")

pdf("plots/milo_umap_celltype_group.pdf", width=15, height=10)
DimPlot(obj, reduction="umap", group.by="seurat_clusters", label=TRUE) + ggtitle("Seurat Clusters")
DimPlot(obj, reduction="umap", group.by="predicted.celltype.l2", label=TRUE) + ggtitle("L2 Azimuth Labels")
DimPlot(obj, reduction="umap", group.by="Group") + ggtitle("Group")
dev.off()

obj_milo <- buildGraph(obj_milo, k=80, d=30)
obj_milo <- makeNhoods(obj_milo, prop = 0.2, k = 80, d=30, refined = TRUE, refinement_scheme="graph")

obj_milo <- countCells(obj_milo, meta.data = data.frame(colData(obj_milo)), samples = "sample")


# how many cells are in each neighborhood
pdf("plots/milo_Nhood_sizes_hist.pdf")
plotNhoodSizeHist(obj_milo)
dev.off()

obj_milo$Group <- gsub(" ", "", obj_milo$Group)
obj_milo$Group <- factor(obj_milo$Group)

obj_design <- data.frame(colData(obj_milo))[, c("sample", "Group")]
obj_design <- distinct(obj_design)
rownames(obj_design) <- obj_design$sample
obj_design <- obj_design[colnames(nhoodCounts(obj_milo)), , drop=FALSE]
table(obj_design$Group)

# making contrasts between groups
obj_milo$Group <- gsub(" ", "", obj_milo$Group)
obj_milo$Group <- factor(obj_milo$Group)
contrast_all <- c(
    "GroupSCZHigh - GroupSCZLow",
    "GroupSCZHigh - GroupHC",
    "GroupSCZLow - GroupHC"
)
model <- model.matrix(~ 0 + Group, data=obj_design)
mod.contrast <- makeContrasts(contrasts=contrast_all, levels=model)

contrast1.res <- testNhoods(obj_milo, design=~0+ Group, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = contrast_all)
table(contrast1.res$SpatialFDR < 0.1)

contrastHvL <- testNhoods(obj_milo, design=~0+ Group, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("GroupSCZHigh - GroupSCZLow"))

obj_milo <- buildNhoodGraph(obj_milo)

pdf("plots/milo_results_umap_nhood_graph.pdf", width=30, height=15)
DimPlot(obj, reduction="umap", group.by="predicted.celltype.l2", label=TRUE) + plotNhoodGraphDA(obj_milo, contrastHvL, alpha=0.3)
dev.off()



