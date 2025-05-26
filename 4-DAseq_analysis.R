# Rowan K. O'Hara-Payne
# 5/25/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Do differential abundance analysis


# Read in libraries and data ----
library(Seurat)
library(dplyr)
library(edgeR)
library(ggplot2)

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
#res_all$cell_type <- gsub("...\\d", "", res_all$cell_type)

# Plots ----

## volcano plot ----
pdf("plots/DAseq_volcano.pdf")
ggplot(res_all, aes(x = logFC, y = -log10(FDR), label = cell_type, color = comparison)) +
  geom_point() +
  ggrepel::geom_text_repel(data = subset(res_all, FDR < 0.05)) +
  facet_wrap(~comparison) +
  theme_bw() +
  labs(title = "DA of Cell Types Across Groups",
       x = "logFC", y = "-log10(FDR)")
dev.off()


## heatmap ----
lfc_matrix <- res_all %>%
  select(cell_type, logFC, comparison) %>%
  pivot_wider(names_from = comparison, values_from = logFC) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

# filling NAs with zero to allow visualization
lfc_matrix_filled <- lfc_matrix
lfc_matrix_filled[is.na(lfc_matrix_filled)] <- 0

pdf("plots/DAseq_heatmap.pdf")
pheatmap(lfc_matrix_filled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "logFC of Cell Types (All Comparisons, NA=0)")
dev.off()


## boxplots ----

# convert to long format and compute proportions
prop_df <- counts_table %>%
  as.data.frame() %>%
  rename(sample = Var1, cell_type = Var2, count = Freq) %>%
  group_by(sample) %>%
  mutate(total = sum(count),
         prop = count / total) %>%
  ungroup()

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


# compute proportions manually for plotting
counts_table$sample <- gsub("_*", "", counts_table$sample)
plot_df <- as.data.frame(counts_table) %>%
  rename(sample = Var1, cell_type = Var2, count = Freq) %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count))
plot_df$sample <- substring(plot_df$sample, 1, 7)
plot_df <- left_join(plot_df, meta_filtered, by = join_by("sample" == "MPRCID"))

plot_df <- as.data.frame(counts_table) %>%
  rename(sample = Var1, cell_type = Var2, count = Freq) %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count)) %>%
  left_join(meta_filtered, by = join_by("sample" == "MPRCID"))

# plot using ggplot2
pdf("plots/celltype_sample_DAseq.pdf")
ggplot(plot_df, aes(x = sample, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Sample", y = "Proportion of Cells", fill = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x")
dev.off()

pdf("plots/celltype_prop_DAseq.pdf")
ggplot(plot_df, aes(x = cell_type, y = prop, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Cell Type Proportion", x = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


