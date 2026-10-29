# Rowan K. O'Hara-Payne
# 5/25/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Do differential abundance analysis


# Read in libraries and data ----
library(Seurat)
library(tidyverse)
library(edgeR)
library(pheatmap)
library(emmeans)
library(ggpubr)

# read in Seurat object
obj <- readRDS('kelly_postqc_clusters_meta.rds')

# read in pseudobulk data
elist <- readRDS('dkelly_singlets.pseudobulk.rds')

# read in meta data
sample <- colnames(elist[[11]])
MPRCID <- gsub( '\\_(.*)' , '' , sample )
MPRCID <- gsub( '\\-(.*)' , '' , MPRCID )

meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
rownames(meta) <- meta$MPRCID
meta <- meta[ MPRCID , ]


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
counts_table <- table(obj$sample, obj$seurat_clusters)
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
design <- model.matrix(~ 0 + Group + Age + Sex, data = sample_metadata)
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

# 19 at FDR < 0.05
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


### ANOVA ----
run_binom_model <- function(df_ct) {

  ct_name <- as.character(unique(df_ct$cell_type)[1])

  group_levels <- c("SCZHigh", "SCZLow", "HC")
  df_ct$Group <- factor(df_ct$Group, levels = group_levels)

  # Fit model
  fit <- glm(
    cbind(count, total - count) ~ Group + Sex + Age,
    family = binomial,
    data = df_ct
  )

  # Overall test
  anova_res <- car::Anova(fit, type = 3)
  p_group <- anova_res["Group", "Pr(>Chisq)"]

  # Emmeans (response scale)
  emm <- emmeans(fit, ~ Group, type = "response")
  emm_df <- as.data.frame(emm) %>%
    mutate(cell_type = ct_name)


  # Split contrast safely
  contr <- emmeans::contrast(emm, method = "pairwise", adjust = "none")

  pairs_df <- as.data.frame(contr) %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " / ") %>%
    dplyr::mutate(
      group1 = trimws(group1),
      group2 = trimws(group2),
      p_label = sprintf("p=%.3g", p.value),
      cell_type = ct_name
    )
  

  # y position for significance bars
  y_max <- max(c(df_ct$prop, emm_df$prob), na.rm = TRUE)
  y_range <- diff(range(c(df_ct$prop, emm_df$prob), na.rm = TRUE))

  pairs_df <- pairs_df %>%
    dplyr::arrange(group1, group2) %>%
    dplyr::mutate(
      y.position = y_max + 0.08 * y_range * seq_len(n())
    )

  # Plot
  p <- ggplot(df_ct, aes(x = Group, y = prop)) +
	  geom_col(data = emm_df,
                aes(x = Group, y = prob, fill = Group), color = "black", lwd = 1.25, width = 0.75) +
	  geom_jitter(width = 0.15, size = 4, alpha = 0.8, color = "black") +
	  geom_errorbar(
    		data = emm_df,
    		aes(x = Group, ymin = asymp.LCL, ymax = asymp.UCL),
    		width = 0.35,
		lwd = 1.25,
    		inherit.aes = FALSE) +
	  scale_x_discrete(
			   labels = c("SCZHigh" = "AGA+ SRD",
				      "SCZLow" = "AGA- SRD",
				      "HC" = "AGA- Control")
			   ) +
	  scale_fill_manual(values = c("#BB5566", "#DDAA33", "#004488")) +
	  labs(title = parse(text = celltype_labels[[ct_name]]),
	       y = "Proportion",
	       fill = "group") +
	  theme_classic() +
	  scale_y_continuous(limits = c(0, NA),
			     expand = expansion(mult = c(0, 0.05))
			     ) +
	  theme(plot.title = element_text(hjust = 0.5,
					  size = 30),
		axis.text = element_text(size = 20),
		axis.title = element_text(size = 17),
		axis.line = element_line(linewidth = 1.25),
		axis.ticks = element_line(linewidth = 1.25))

  # Add significance bars (AFTER p is created)
  p <- p + ggpubr::stat_pvalue_manual(
				      pairs_df,
    				      label = "p_label",
				      xmin = "group1",
				      xmax = "group2",
				      y.position = "y.position",
				      tip.length = 0.015,
				      inherit.aes = FALSE,
				      size = 5,
				      bracket.size = 0.7) +
			coord_cartesian(clip = "off")

  return(list(
    summary = data.frame(cell_type = ct_name, p_value = p_group),
    emm = emm_df,
    pairs = pairs_df,
    plot = p
  ))
}

results_list <- prop_df %>%
  group_split(cell_type) %>%
  map(run_binom_model)

summary_df <- map_dfr(results_list, "summary") %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

emm_all <- map_dfr(results_list, "emm")

pairs_all <- map_dfr(results_list, "pairs") %>%
  mutate(FDR = p.adjust(p.value, method = "fdr"))

pdf("plots/celltype_proportion_plots.pdf", width = 12, height = 8)
walk(results_list, ~ print(.x$plot))
dev.off()


### t-test ----
comparisons <- list(
  c("HC", "SCZHigh"),
  c("SCZHigh", "SCZLow"),
  c("HC", "SCZLow")
)
pdf("plots/DAseq_boxplots.pdf")
for (ct in levels(obj@meta.data$seurat_clusters)) {
	prop_ct <- filter(prop_df, cell_type == ct)
	p <- ggplot(prop_ct, aes(x = Group, y = prop, fill = Group)) +
  		geom_boxplot() +
		geom_point() +
		stat_compare_means(comparisons = comparisons,
				   method = "t.test", 
				   label = "p.signif") +
  		theme_bw() +
  		labs(y = "Proportion", title = paste(ct)) +
  		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	print(p)
	 }
dev.off()

pdf("plots/DAseq_scatterplots.pdf")
for (ct in levels(obj@meta.data$seurat_clusters)) {
        prop_ct <- filter(prop_df, cell_type == ct)
        p <- ggplot(prop_ct, aes(x = AGA_level, y = prop)) +
                geom_point(aes(color = Group)) +
		geom_smooth(method = "lm", formula = y ~ x, color = "black") +
        	stat_cor(aes(label = after_stat(rr.label)),
                 label.x.npc = 0.75, label.y.npc = 0.00001) +
                coord_cartesian(clip = "off") +
		theme_bw() +
                labs(y = "Proportion", x = "AGA Level", title = paste(ct)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
         }
dev.off()


pdf("plots/celltype_prop_DAseq.pdf", width = 15, height = 6)
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
  scale_fill_manual(values = col_pal) +
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
reducedDim(obj_milo, "PCA") <- Embeddings(obj, "pca")[, 1:30]

# plot variancce explained by PCs -- will help choose our d parameter
pdf("plots/elbow_plot_PCs_kelly.pdf")
ElbowPlot(obj, ndims = 50)
dev.off()

pca_stdev <- obj[["pca"]]@stdev
var_explained <- pca_stdev^2 / sum(pca_stdev^2)
pdf("plots/var_explained_PCs_kelly.pdf")
plot(var_explained[1:50], type = "b",
   ylab = "Proportion variance explained",
   xlab = "PC")
dev.off()


pdf("plots/milo_umap_celltype_group.pdf", width=15, height=10)
DimPlot(obj, reduction="umap", group.by="seurat_clusters", label=TRUE) + ggtitle("Seurat Clusters")
DimPlot(obj, reduction="umap", group.by="predicted.celltype.l2", label=TRUE) + ggtitle("L2 Azimuth Labels")
DimPlot(obj, reduction="umap", group.by="Group") + ggtitle("Group")
dev.off()

k_list <- c(25, 30, 35, 45, 50, 55, 70, 75)
d_list <- c(10, 20, 30, 40, 50)

summary_list <- c()
plot_list <- c()
results_obj <- c()

for (k in k_list) {
	for (d in d_list) {
		obj_milo_iteration <- obj_milo
		label <- paste(k, d, sep='_')
		obj_milo_iteration <- buildGraph(obj_milo_iteration, k = k, d = d)
		obj_milo_iteration <- makeNhoods(obj_milo_iteration, prop = 0.2, k = k, d = d, refined = TRUE, refinement_scheme="graph")

		obj_milo_iteration <- countCells(obj_milo_iteration, meta.data = data.frame(colData(obj_milo_iteration)), samples = "sample")


		# how many cells are in each neighborhood
		#p1 <- plotNhoodSizeHist(obj_milo_iteration)


		obj_design <- data.frame(colData(obj_milo_iteration))[, c("sample", "Group")]
		obj_design <- distinct(obj_design)
		rownames(obj_design) <- obj_design$sample
		obj_design <- obj_design[colnames(nhoodCounts(obj_milo_iteration)), , drop=FALSE]
		table(obj_design$Group)

		# making contrasts between groups
		contrast_all <- c(
    			"GroupSCZHigh - GroupSCZLow",
    			"GroupSCZHigh - GroupHC",
    			"GroupSCZLow - GroupHC"
		)
		model <- model.matrix(~ 0 + Group, data=obj_design)
		mod.contrast <- makeContrasts(contrasts=contrast_all, levels=model)

		contrast1.res <- testNhoods(obj_milo_iteration, design=~0+ Group, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = contrast_all)
		table(contrast1.res$SpatialFDR < 0.1)

		contrastHvL <- testNhoods(obj_milo_iteration, design=~0+ Group, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("GroupSCZHigh - GroupSCZLow"))
		table(contrastHvL$SpatialFDR < 0.1)

		obj_milo_iteration <- buildNhoodGraph(obj_milo_iteration)

		p2 <- plotNhoodGraphDA(obj_milo_iteration, contrastHvL, alpha=0.3) + ggtitle(label)
		#plot_list[[label]] <- p1 / p2
		plot_list[[label]] <- p2

		results_obj[[label]] <- obj_milo_iteration
		summary_res <- cbind(k, d, annot_df)
                summary_list[[label]] <- summary_res

	}
}


# save table
results_all <- bind_rows(summary_list)
write_csv(results_all, "milo_results_nhood_summary.csv")

# save objects
saveRDS(results_obj, "milo_obj_list.RDS")


# plot results
plot_info <- data.frame(
  name = names(plot_list),
  k = as.numeric(str_extract(names(plot_list), "^[0-9]+")),
  d = as.numeric(str_extract(names(plot_list), "(?<=_)[0-9]+"))
)

plot_info$plot <- plot_list

plot_info <- plot_info[order(plot_info$k, plot_info$d), ]

k_vals <- sort(unique(plot_info$k))
d_vals <- sort(unique(plot_info$d))

plot_matrix <- lapply(k_vals, function(k) {
  row_plots <- plot_info$plot[plot_info$k == k]
  # reorder within row by d
  row_plots <- row_plots[order(plot_info$d[plot_info$k == k])]
  row_plots
})

# plot a grid of all of the different graphs I made for comparisons
pdf("plots/milo_k_by_d_grid3.pdf", width = 5 * length(d_vals), height = 5 * length(k_vals))
cowplot::plot_grid(
  plotlist = unlist(plot_matrix, recursive = FALSE),
  nrow = length(k_vals),
  ncol = length(d_vals),
  labels = names(plot_info$plot),
  label_size = 8
)
dev.off()

# plot an individual milo graph next to the Azimuth labeled UMAP
pdf("plots/milo_results_umap_nhood_graph.pdf", width=30, height=15)
DimPlot(obj, reduction="umap", group.by="predicted.celltype.l2", label=TRUE) + plotNhoodGraphDA(obj_milo, contrastHvL, alpha=0.3)
dev.off()



## check if these look right -----

### check neighborhood sizes ----
nhood_df <- do.call(rbind, lapply(names(results_obj), function(nm) {
 	split_matrix <- str_split_fixed(nm, "_", n = 2)
	k_val <- split_matrix[1, 1]
	d_val <- split_matrix[1, 2]
	sizes <- rowSums(results_obj[[nm]]@nhoodCounts)
  data.frame(k = as.numeric(k_val), d = as.numeric(d_val), size = sizes)
}))


pdf("plots/nhood_size_by_kd.pdf", width=25, height=30)
ggplot(nhood_df, aes(x = factor(1), y = size)) +
  geom_violin(fill = "grey80") +
  geom_boxplot(width = 0.1) +
  facet_grid(k ~ d, switch = "both") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("Neighborhood size")
dev.off()

### Using conservative parameters to reduce false positives ----
sig_table <- results_all %>% group_by(k, d) %>% filter(predicted.celltype.l2 == "CD8 TEM") %>% summarize(mean = mean(FDR < 0.3, na.rm = TRUE))

pdf("plots/frac_sig_CD8TEM_by_k.pdf")
plot(as.numeric(sig_table$k), sig_table$mean, type='b',
     xlab="k", ylab="fraction significant in CD8 TEM pop")
dev.off()


### chose k = 45, d = 20 ----
obj_milo_conserv <- buildGraph(obj_milo, k = 45, d = 20)
obj_milo_conserv <- makeNhoods(obj_milo_conserv, prop = 0.2, k = 45, d = 20, refined = TRUE, refinement_scheme="graph")
obj_milo_conserv <- countCells(obj_milo_conserv, meta.data = data.frame(colData(obj_milo_conserv)), samples = "sample")

pdf("plots/NHoodSize_milo_kelly_45.pdf")
plotNhoodSizeHist(obj_milo_conserv)
dev.off()

obj_design <- data.frame(colData(obj_milo_conserv))[, c("sample", "Group", "Sex", "Age")]
obj_design <- distinct(obj_design)
rownames(obj_design) <- obj_design$sample
obj_design <- obj_design[colnames(nhoodCounts(obj_milo_conserv)), , drop=FALSE]
table(obj_design$Group)


model <- model.matrix(~ 0 + Group + Sex + Age, data=obj_design)
contrast_all <- c(
                        "GroupSCZHigh - GroupSCZLow",
                        "GroupSCZHigh - GroupHC",
                        "GroupSCZLow - GroupHC"
                )
mod.contrast <- makeContrasts(contrasts=contrast_all, levels=model)
obj_milo_conserv <- calcNhoodDistance(obj_milo_conserv, d=20)
contrast1.res <- testNhoods(obj_milo_conserv, design=model, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = contrast_all)
contrastHvL <- testNhoods(obj_milo_conserv, design=model, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("GroupSCZHigh - GroupSCZLow"))

obj_milo_conserv <- buildNhoodGraph(obj_milo_conserv)
pdf("plots/NhoodGraph_milo_kelly.pdf")
plotNhoodGraphDA(obj_milo_conserv, contrastHvL, alpha=0.30) + ggtitle("k = 45, d = 20, SpatialFDR < 0.30")
plotNhoodGraphDA(obj_milo_conserv, contrastHvL, alpha=0.15) + ggtitle("k = 45, d = 20, SpatialFDR < 0.15")
dev.off()

annot_df <- annotateNhoods(obj_milo_conserv, contrastHvL, "seurat_clusters")
write_csv(annot_df, "milo_kelly_results_summary.csv")

annot_df$seurat_clusters <- ifelse(annot_df$seurat_clusters_fraction < 0.7, "Mixed", annot_df$seurat_clusters)
pdf("plots/DA_kelly_frac_neighborhoods_45.pdf")
ggplot(annot_df, aes(seurat_clusters_fraction)) + geom_histogram(bins=50)
plotDAbeeswarm(annot_df, group.by = "seurat_clusters")
dev.off()

pdf("plots/DA_results_milo_kelly.pdf")
col_pal <- c("darkturquoise", "darkslategrey", "cadetblue2", "blue3", "darkcyan", "cyan", "aquamarine3", "cornflowerblue", "dodgerblue3", "deepskyblue4", "darkblue", "royalblue3", "lightblue", "slateblue1", "navy", "powderblue", "seagreen", "limegreen", "olivedrab", "springgreen4", "yellow", "violetred1", "tomato", "tan1", "salmon", "rosybrown1", "red3", "palevioletred2", "orange", "maroon", "magenta3", "coral3", "gold", "firebrick", "deeppink", "darkorange3", "brown2")
ggplot(annot_df, aes(PValue)) + geom_histogram(bins=50)
ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
  geom_point(aes(color = seurat_clusters)) +
  scale_color_manual(values = col_pal) +
  geom_hline(yintercept = 1) + ## Mark significance threshold (10% FDR)
  geom_hline(yintercept = 3, linetype = "dashed")
dev.off()

#contrastHvL <- groupNhoods(obj_milo_conserv, contrastHvL, max.lfc.delta = 1)
#nhood.dge <- testDiffExp(x=obj_milo_conserv, da.res=contrastHvL, model.contrasts = c("GroupSCZHigh - GroupSCZLow"), design=~ 0 + Group + Sex + Age, meta.data=obj@meta.data)


# Create CIBERSORTx input ----
library(Matrix)
obj <- readRDS("kelly_postqc_clusters_meta.rds")
Idents(obj) <- "seurat_clusters"

# filter cells to the top markers for each cell type
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 100)
marker_genes <- unique(top_markers$gene)  # 1871 genes

# average expression per cell type in our data set
avg_exp <- AggregateExpression(
  obj,
  group.by = "seurat_clusters",
  assays = "RNA",
  slot = "data"  # use log-normalized data
)$RNA

# remove mitochondrial and ribosomal genes
avg_exp <- avg_exp[marker_genes, ]  # 1871 genes
avg_exp <- avg_exp[!grepl("^MT-", rownames(avg_exp)), ]  # 1864 genes
avg_exp <- avg_exp[!grepl("^RPS|^RPL", rownames(avg_exp)), ]  # 1794 genes

# wrangle Glahn data set into matrix
bulk <- read_csv("salmon_merged_gene_counts.csv")
gene_id <- bulk$gene_id
bulk_matrix <- as.matrix(bulk[ , -1])
rownames(bulk_matrix) <- gene_id

# check that the same genes are present in both
common_genes <- intersect(rownames(avg_exp), rownames(bulk_matrix))
final_signature <- final_signature[common_genes, ]  # 1595 genes
bulk_matrix <- bulk_matrix[common_genes, ]

# export for CIBERSORTx
final_export <- cbind(gene = rownames(final_signature), as.data.frame(final_signature))
write.table(
  final_export,
  file = "kelly_signature_matrix.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# normalize Glahn counts to prep for CIBERSORTx
dge <- DGEList(counts = bulk_matrix)
dge <- calcNormFactors(dge)
tpm <- cpm(dge, normalized.lib.sizes = TRUE)
tpm_export <- cbind(gene = rownames(tpm), as.data.frame(tpm))
write.table(
  tpm_export,
  file = "glahn_mixture_file.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# read in Pan
pan_counts <- read.table(file="/local/projects/idea/rohara/immunophenotype_SCZ/GSE243841_cleaned_counts.tsv", sep="\t", header=T)
rownames(pan_counts) <- pan_counts$Geneid
bulk_matrix <- as.matrix(pan_counts[ , -1])
genes_clean <- sub("\\.\\d+$", "", rownames(pan_counts))

library(org.Hs.eg.db)
library(AnnotationDbi)
gene_map <- mapIds(org.Hs.eg.db,
                  keys = genes_clean,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
idx <- which(!is.na(gene_map))
found_genes <- gene_map[idx]
pan_counts <- pan_counts[idx, ]
ind <- !duplicated(found_genes)
pan_counts <- pan_counts[ind, ]
pan_counts$Geneid <- found_genes[ind]
rownames(pan_counts) <- pan_counts$Geneid

common_genes <- intersect(rownames(final_signature), rownames(bulk_matrix))
final_signature <- final_signature[common_genes, ]  # 1692 genes
bulk_matrix <- bulk_matrix[common_genes, ]

# normalize Pan counts to prep for CIBERSORTx
dge <- DGEList(counts = bulk_matrix)
dge <- calcNormFactors(dge)
tpm <- cpm(dge, normalized.lib.sizes = TRUE)
tpm_export <- cbind(gene = rownames(tpm), as.data.frame(tpm))
write.table(
  tpm_export,
  file = "pan_mixture_file.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

final_export <- cbind(gene = rownames(final_signature), as.data.frame(final_signature))
write.table(
  final_export,
  file = "kelly_signature_matrix_pan.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
 )


# look at CIBERSORT output

# Glahn
prop <- read.table("CIBERSORTx_Glahn_JustKelly.txt", sep="\t", header=T)
colnames(prop)[1] <- "sample"
prop <- prop %>%
	pivot_longer(cols = starts_with("g"),
		     names_to = "celltype",
		     values_to = "prop")

meta <- read_csv("redcap_selection_for_seth.csv")
meta$sample <- paste0(meta$subject_type, meta$id)
meta$group <- gsub("^\\d+\\s*=\\s*", "", meta$subject)
meta$sex <- as.factor(meta$sex)
meta$group <- as.factor(meta$group)

glahn_prop <- merge(prop, meta, by = "sample")
pdf("plots/DAseq_boxplots_glahn.pdf")
for (ct in unique(glahn_prop$celltype)) {
        prop_ct <- glahn_prop %>% filter(celltype == ct)
        p <- ggplot(prop_ct, aes(x = group, y = prop, fill = group)) +
                geom_boxplot() +
                geom_point() +
                stat_compare_means(method = "t.test", label = "p.signif") +
                theme_bw() +
                labs(y = "Proportion", title = paste(ct)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
         }
dev.off()


# load BPRS scores
bprs <- read_excel("GeneNetworksInfluenc_DATA_2021-01-11_1019.xlsx")
bprs$sample <- paste0(bprs$is_this_subject_a_patient, bprs$id)
bprs <- bprs %>% filter(!str_detect(sample, "--"))
glahn <- left_join(meta, bprs, by = "sample")
glahn <- glahn[ , c("sample", "group", "age.x", "sex.x", "bprs_has_this_instrument", "somatic_concerns", "anxiety", "depression", "guilt", "hostility", "suspiciousness", "unusual_thought_content", "grandiosity", "hallucinations", "disorientation", "conceptual_disorganization", "excitement", "motor_retardation", "blunted_affect", "tension", "mannerisms_and_posturing", "uncooperativeness", "emotional_withdrawal", "suicidality", "self_neglect", "bizarre_behaviors", "elevated_mood", "motor_hyperactivity", "distractability", "bprse_complete")]
glahn$neg_avg <- (glahn$blunted_affect + glahn$emotional_withdrawal + glahn$motor_retardation) / 3
glahn$pos_avg <- (glahn$unusual_thought_content + glahn$conceptual_disorganization + glahn$hallucinations + glahn$grandiosity) / 4
glahn <- merge(glahn, prop, by = "sample")

celltypes <- colnames(glahn)[33:69]
scores <- c("neg_avg", "pos_avg")

# center log-ratio transform
glahn[, celltypes] <- glahn[, celltypes] + 1e-6
df_clr <- glahn %>%
  rowwise() %>%
  mutate(
    gm = exp(mean(log(c_across(all_of(celltypes)))))
  ) %>%
  mutate(
    across(all_of(celltypes), ~ log(.x / gm))
  ) %>%
  ungroup() %>%
  select(-gm)

df_clr <- df_clr %>%
  filter(group == "Patient")
df_resid <- df_clr

residualize_vec <- function(y, data) {
  fit <- lm(y ~ age.x + sex.x, data = data, na.action = na.exclude)
  resid(fit)
}

# residualize celltypes
for (ct in celltypes) {
  df_resid[[ct]] <- residualize_vec(df_clr[[ct]], df_clr)
}

# residualize symptoms
for (sc in scores) {
  df_resid[[sc]] <- residualize_vec(df_clr[[sc]], df_clr)
}

cor_df <- expand.grid(celltype = celltypes, symptom = scores) %>%
  rowwise() %>%
  mutate(
    cor = cor(
      df_resid[[celltype]],
      df_resid[[symptom]],
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  ) %>%
  ungroup()

cor_df <- cor_df %>%
  rowwise() %>%
  mutate(
    test = list(cor.test(
      df_resid[[celltype]],
      df_resid[[symptom]],
      method = "spearman"
    )),
    p = test$p.value
  ) %>%
  ungroup() %>%
  mutate(sig = p < 0.05,
	 label = sprintf("%.2f%s", cor, ifelse(sig, "*", "")))


pdf("plots/BPRS_correlations.pdf")
ggplot(cor_df, aes(x = symptom, y = celltype, fill = cor)) +
  geom_tile() +
  geom_text(aes(label = label)) +
  scale_fill_gradient2(midpoint = 0) +
  theme_minimal() +
  labs(
    title = "Cell Type Proportion Correlations with BPRS Subscales (Patients Only)",
    fill = "Spearman r"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# correlations with Kelly
prop_df$cell_type <- paste0("g", prop_df$cell_type)
df <- prop_df %>%
  select(sample, cell_type, prop) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = prop 
  )
df <- df %>% mutate(MPRCID = sub("[_-].*", "", sample))
celltypes <- unique(prop_df$cell_type)
scores <- c("neg_avg", "pos_avg")

sample_metadata <- sample_metadata %>% mutate(MPRCID = sub("[_-].*", "", sample))

library(readxl)
bprs <- read_excel("BPRS-all-entries.xlsx")
bprs$sample <- bprs$MPRCID
bprs <- bprs %>%
  mutate(date = as.Date(Interview_date))
df_latest <- bprs %>%
  group_by(sample) %>%
  slice_max(order_by = date, n = 1, with_ties = FALSE) %>%
  ungroup()
df_latest$neg_avg <- (df_latest$"16_Blunted" + df_latest$"3_Emot_Withdr" + df_latest$"13_MotorRetard") / 3
df_latest$pos_avg <- (df_latest$"15_Unusualthght" + df_latest$"4_ConceptDisorg" + df_latest$"12_Hallucination" + df_latest$"8_Grandiosity") / 4

kelly <- merge(df_latest, df, by = "MPRCID")
kelly <- merge(kelly, sample_metadata, by = "MPRCID")
kelly[, celltypes] <- kelly[, celltypes] + 1e-6
df_clr <- kelly %>%
  rowwise() %>%
  mutate(
    gm = exp(mean(log(c_across(all_of(celltypes)))))
  ) %>%
  mutate(
    across(all_of(celltypes), ~ log(.x / gm))
  ) %>%
  ungroup() %>%
  select(-gm)
aga_levels <- c(124.40, 0.95, 6.40, 129.50, 23.64, 31.87, 32.10, 84.52, 4.82, 22.97)
df_clr <- cbind(df_clr, aga_levels)

df_resid <- df_clr
residualize_vec <- function(y, data) {
  fit <- lm(y ~ Age + Sex, data = data, na.action = na.exclude)
  resid(fit)
}
for (ct in celltypes) {
  df_resid[[ct]] <- residualize_vec(df_clr[[ct]], df_clr)
}
for (sc in scores) {
  df_resid[[sc]] <- residualize_vec(df_clr[[sc]], df_clr)
}
df_resid$aga_levels_resid <- residualize_vec(aga_levels, df_clr)

# cell type proportions vs BPRS subscales
cor_df <- crossing(celltype = celltypes, symptom = scores) %>%
  rowwise() %>%
  mutate(
    cor = cor(
      df_resid[[celltype]],
      df_resid[[symptom]],
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  ) %>%
  ungroup()
cor_df <- cor_df %>%
  rowwise() %>%
  mutate(
    test = list(cor.test(
      df_resid[[celltype]],
      df_resid[[symptom]],
      method = "spearman"
    )),
    p = test$p.value
  ) %>%
  ungroup() %>%
  mutate(sig = p < 0.05,
	 label = sprintf("%.2f%s", cor, ifelse(sig, "*", "")))
pdf("plots/BPRS_correlations_Kelly.pdf")
ggplot(cor_df, aes(x = symptom, y = celltype, fill = cor)) +
  geom_tile() +
  geom_text(aes(label = label)) +
  scale_fill_gradient2(midpoint = 0) +
  theme_minimal() +
  labs(
    title = "Cell Type Proportion Correlations with BPRS Subscales",
    fill = "Spearman r"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


# cell type proportions vs AGA levels correlations
cor_df <- tibble(celltype = celltypes) %>%
  rowwise() %>%
  mutate(
    cor = cor(df_resid[[celltype]], df_resid$aga_levels_resid,
              method = "spearman", use = "pairwise.complete.obs"),
    test = list(cor.test(df_resid[[celltype]], df_resid$aga_levels_resid,
                         method = "spearman")),
    p = test$p.value
  ) %>%
  ungroup() %>%
  mutate(
    log10_p = -log10(p),
    sig = p < 0.05,
    label = sprintf("%.2f%s", cor, ifelse(sig, "*", ""))
  )

cor_df$celltype <- gsub("g", "", cor_df$celltype)
celltype_labels <- c("0" = "SOD^'Hi'*CD14^'+'~Mono", "1" = "'IL-1B'^'Hi'*CD14^'+'~Mono", "10" = "MHC^'Hi'*CD14^'+'~Mono", "12" = "C1Q^'Hi'*CD16^'+'~Mono", "14" = "CD14^'+'~Mono", "24-1" = "CD16^'+'~Mono", "25" = "ISG^'+'*CD14^'+'~Mono", "27" = "Platelets", "21-0" = "DC3", "21-1" = "DC2", "21-2" = "DC1", "26" = "pDC", "29" = "Plasmablasts", "28" = "Ki67^'+'*CD8^'+'*NKT['C1']", "30" = "HSPCs", "32" = "Neutrophils", "24-0" = "Atypical~Memory~B", "15" = "Switched~Memory~B", "19" = "CD5^'+'~B", "11" = "Activated~Naive~B", "18" = "CD56^'bright'~NK", "2" = "CD56^'dim'~NK", "9" = "CD8^'+'*NKT['C1']*A['T']", "17" = "CD161^'+'*CD8^'+'*NKT['C1']*A['T']", "22" = "CD8^'+'*NKT['C1']*A['P']", "8" = "CD8^'+'*T['C1']*M", "13-1" = "CD8^'+'*'MAIT-like'*T['C1']*M", "13-0" = "CD4^'+'*T['H17']*M", "23" = "CD4^'+'*T['FH1']*M", "16" = "CD4^'+'*TM", "6" = "CD4^'+'*T['H2']*M['P']", "31" = "CD4^'-'*CD8^'dim'*T['Reg']*M", "7" = "CD4^'+'*T['FH']*M['P']", "20" = "CD4^'+'*T['Reg']*M", "5" = "CD4^'+'*TN", "3" = "CD26^'Hi'*CD4^'+'*TN", "4" = "CD8^'+'*TN")

pdf("plots/BPRS_correlations_aga_kelly.pdf", width = 10, height = 6)
ggplot(cor_df, aes(x = cor, y = log10_p)) +
  geom_point(aes(color = sig), size = 2) +
  geom_text_repel(
    data = cor_df %>% filter(sig),
    aes(label = celltype_labels[celltype]),
    parse = TRUE,
    size = 5
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0) +
  scale_color_manual(values = c("grey70", "red")) +
  scale_x_continuous(
  limits = c(-1, 1),
  breaks = seq(-1, 1, by = 0.25)
) +
  theme_minimal() +
  labs(
    x = expression(paste("Spearman correlation (", rho, ")")),
    y = "-log10(p-value)",
    color = "Significant",
    title = "Correlations Between Cell Type Proportions and AGA Levels") +
  theme(plot.title = element_text(size = 15))
dev.off()

ct <- "g13-0"
ct_label <- parse(text = celltype_labels[gsub("g", "", ct)])[[1]]
cor_val <- cor_df %>%
  filter(celltype == ct) %>%
  pull(cor)
p_val <- cor_df %>%
  filter(celltype == ct) %>%
  pull(p)
pdf("plots/corr_13-0_aga_scatterplot.pdf", width = 12, height = 8)
ggplot(df_resid, aes(x = aga_levels_resid, y = `g13-0`)) +
	geom_point() +
	geom_smooth(method = "lm", se = FALSE) +
	labs(x = "AGA Levels (adjusted)",
	     y = bquote(.(ct_label) ~ "proportion (center log-ratio transformed & adjusted)"),
	     title = bquote("Correlation between" ~ .(ct_label) ~ "proportion and AGA levels (age + sex adjusted)")) +
	theme_minimal() +
	annotate("text",
		 x = Inf,
		 y = -Inf,
		 label = bquote(rho == .(round(cor_val, 3)) * "," ~ italic(p) == .(signif(p_val, 2))),
		 hjust = 2, 
		 vjust = -6,
		 size = 6) +
	theme(plot.title = element_text(size = 20,
					hjust = 0.5),
	      axis.title = element_text(size = 15))

dev.off()



# Pan
pan <- read.table("CIBERSORTx_Pan_JustKelly.txt", sep="\t", header=T)
colnames(pan)[1] <- "sample"
pan <- pan %>%
        pivot_longer(cols = starts_with("g"),
                     names_to = "celltype",
                     values_to = "prop")

meta <- read_csv("meta_sex_group__Pan.csv")
pan_prop <- merge(pan, meta, by = "sample")
pdf("plots/DAseq_boxplots_Pan.pdf")
for (ct in unique(pan_prop$celltype)) {
        prop_ct <- pan_prop %>% filter(celltype == ct)
        p <- ggplot(prop_ct, aes(x = group, y = prop, fill = group)) +
                geom_boxplot() +
                geom_point() +
                stat_compare_means(method = "t.test", label = "p.signif") +
                theme_bw() +
                labs(y = "Proportion", title = paste(ct)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
         }
dev.off()


# look at the proportions for all studies
glahn_prop$study <- "Glahn"
pan_prop$study <- "Pan"


pdf("plots/CIBERSORTx_prop_kelly_glahn_pan.pdf")
for (ct in unique(prop_df$celltype)) {
        prop_ct <- prop_df %>% filter(celltype == ct)
        p <- ggplot(prop_ct, aes(x = study, y = prop, fill = group)) +
                geom_boxplot(position = position_dodge(0.8)) +
                geom_point() +
                stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif") +
                theme_bw() +
                labs(y = "Proportion", title = paste(ct)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
         }
dev.off()
