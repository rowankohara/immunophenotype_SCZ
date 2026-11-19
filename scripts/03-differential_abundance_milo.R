# Rowan K. O'Hara-Payne
# 5/25/2025
# conda environment: rnaseq
# Differential abundance testing using a k-nearest neighborhood algorithm with the package Milo (PMID: 34594043)
# Some steps are very time consuming and computationally intensive; I would suggest running this script as a job
####################################################################################################################


# Load libraries and data ----
library(miloR)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)

obj <- readRDS("kelly_postqc_clusters_meta.rds")

# convert Seurat object into a single cell experiment object, then into a Milo object
obj_sce <- as.SingleCellExperiment(obj)
obj_milo <- Milo(obj_sce)

# get the UMAP from the single cell experiment object
reducedDim(obj_milo, "UMAP") <- reducedDim(obj_sce, "UMAP")
reducedDim(obj_milo, "PCA") <- Embeddings(obj, "pca")[, 1:30]

# create cell type labels and palette for plotting
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
  "32"   = "Neutrophils",
  "Mixed" = "Mixed"
)

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
  "32"   = "mediumpurple3",
  "Mixed" = "black"
)


# Choosing Parameters ----

# plot variance explained by PCs -- will help choose our d parameter
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

# ran Milo with different parameter combinations to confirm that results were robust to reasonable variation in these parameter
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
    
    # create design
    obj_design <- data.frame(colData(obj_milo_iteration))[, c("sample", "Class")]
    obj_design <- distinct(obj_design)
    rownames(obj_design) <- obj_design$sample
    obj_design <- obj_design[colnames(nhoodCounts(obj_milo_iteration)), , drop=FALSE]
    table(obj_design$Group)

    
    contrastHvSC <- testNhoods(obj_milo_iteration, design=~0+ Class, design.df=obj_design_HvSC, fdr.weighting="graph-overlap", model.contrasts = c("ClassAGA+ - ClassAGA-"))
    table(contrastHvSC$SpatialFDR < 0.3)
    
    obj_milo_iteration <- buildNhoodGraph(obj_milo_iteration)
    
    p2 <- plotNhoodGraphDA(obj_milo_iteration, contrastHvL, alpha=0.3) + ggtitle(label)
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
pdf("plots/milo_k_by_d_grid.pdf", width = 5 * length(d_vals), height = 5 * length(k_vals))
cowplot::plot_grid(
  plotlist = unlist(plot_matrix, recursive = FALSE),
  nrow = length(k_vals),
  ncol = length(d_vals),
  labels = names(plot_info$plot),
  label_size = 8
)
dev.off()


# Differential Abundance Analysis ----
# build neighborhood graph with k = 45 and d = 20
obj_milo_conserv <- buildGraph(obj_milo, k = 45, d = 20)
obj_milo_conserv <- makeNhoods(obj_milo_conserv, prop = 0.2, k = 45, d = 20, refined = TRUE, refinement_scheme="graph")
obj_milo_conserv <- countCells(obj_milo_conserv, meta.data = data.frame(colData(obj_milo_conserv)), samples = "sample")

# plot histogram of how many cells are in each neighborhood
pdf("plots/NHoodSize_milo_kelly_45.pdf")
p1 <- plotNhoodSizeHist(obj_milo_conserv)
print(p1)
dev.off()

# calculate the neighborhood distance
obj_milo_conserv <- calcNhoodDistance(obj_milo_conserv, d=20)
saveRDS(obj_milo_conserv, "miloR_conserv_obj_kelly.RDS")


## AGA+ vs AGA- contrast ----
# using this vignette for contrasts https://bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_contrasts.html

# after you create an initial milo object with your parameters, you can test for differential abundance for multiple contrasts
obj_milo_HvSC <- obj_milo_conserv

# create design
obj_design_HvSC <- data.frame(colData(obj_milo_HvSC))[, c("sample", "Class", "Sex", "Age")]  # class = AGA status
obj_design_HvSC <- distinct(obj_design_HvSC)
rownames(obj_design_HvSC) <- obj_design_HvSC$sample
obj_design_HvSC <- obj_design_HvSC[colnames(nhoodCounts(obj_milo_HvSC)), , drop=FALSE]
table(obj_design_HvSC$Class)  # this should be the number of participants per group, not cells

# create model matrix with AGA status, sex, and age as covariates
model_HvSC <- model.matrix(~ 0 + Class + Sex + Age, data=obj_design_HvSC)

# test the number of neighborhoods between contrasts
contrastHvSC <- testNhoods(obj_milo_HvSC, design=model_HvSC, design.df=obj_design_HvSC, fdr.weighting="graph-overlap", model.contrasts = c("ClassAGA+ - ClassAGA-"))

# build neighborhood graph
obj_milo_HvSC <- buildNhoodGraph(obj_milo_HvSC)

### Fig 2A ----
# plot neighborhood graph colored by the increase or decrease of abundance
pdf("plots/NhoodGraph_milo_kelly.pdf", width = 12, height = 8)
plotNhoodGraphDA(obj_milo_HvSC, contrastHvSC, alpha=0.30) +  # alpha = spatial FDR
  labs(caption = "k = 45, d = 20, SpatialFDR < 0.30")
dev.off()

# annotate the neighborhoods based on our cell type annotations
annot_df <- annotateNhoods(obj_milo_HvSC, contrastHvSC, "seurat_clusters")
write_csv(annot_df, "milo_kelly_results_summary.csv")

# neighborhoods needed to be at least 70% of a single cell type to be labeled as such, otherwise labeled as "Mixed"
annot_df$seurat_clusters <- ifelse(annot_df$seurat_clusters_fraction < 0.7, "Mixed", annot_df$seurat_clusters)

# make into data frame and get the absolute value of the maximum logFC for coordinates
annot_df_base <- as.data.frame(annot_df)
max_logFC <- max(abs(annot_df_base$logFC), na.rm = TRUE)

# order the clusters by their median log fold change
cluster_order <- annot_df_base %>%
  group_by(seurat_clusters) %>%
  summarize(med_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_logFC) %>%
  pull(seurat_clusters)
annot_df_base$seurat_clusters <- factor(annot_df_base$seurat_clusters, levels = cluster_order)

pdf("plots/DA_kelly_frac_neighborhoods_hist_45.pdf")
ggplot(annot_df_base, aes(seurat_clusters_fraction)) +
  geom_histogram(bins=50)
dev.off()

### Fig 2B ----
# plot the beeswarm, ordered by the median log fold change
pdf("plots/DA_kelly_frac_neighborhoods_45.pdf", width = 12, height = 8)
plotDAbeeswarm(annot_df_base, group.by = "seurat_clusters", alpha=0.30) +
  labs(x = "Cell Type",
       y = "logFC") +
  scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  coord_flip(ylim = c(-max_logFC, max_logFC)) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"))
dev.off()

annot_df$seurat_clusters <- factor(
  annot_df$seurat_clusters,
  levels = names(celltype_pal)
)

pdf("plots/DA_results_hist_45_HvSC.pdf")
ggplot(annot_df, aes(PValue)) +
  geom_histogram(bins=50)
dev.off()

### Fig 2C ----
# plot the volcano
pdf("plots/DA_results_milo_kelly_HvSC.pdf", width = 12, height = 8)
ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
  geom_point(aes(color = seurat_clusters)) +
  scale_color_manual(labels = function(x) parse(text = celltype_labels[x]), values = celltype_pal) +
  geom_hline(aes(yintercept = -log10(0.15), 
                 linetype = "spatialFDR < 0.15")) +  # Mark significance threshold (15% FDR)
  geom_hline(aes(yintercept = -log10(0.3), 
                 linetype = "spatialFDR < 0.30")) +
  scale_linetype_manual(name = "threshold", 
                        values = c("spatialFDR < 0.15" = "solid", 
                                   "spatialFDR < 0.30" = "dashed")) +
  theme_classic() +
  labs(color = "cell type") +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 17),
        axis.line = element_line(linewidth = 1.25),
        axis.ticks = element_line(linewidth = 1.25))
dev.off()

# find neighborhood markers
nhood.dge <- findNhoodMarkers(obj_milo_HvSC, contrastHvSC, overlap=1, compute.new=TRUE)

saveRDS(nhood.dge, "milo_kelly_45_markers_HvSC.RDS")


## Celiac vs Healthy Controls (Kornberg et al., 2023; PMID: 37450575) ----
# Repeating Milo workflow on scRNA-seq data (GEO Accession: GSE252762) using celiac vs healthy controls contrast

# read in data, create Milo object
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
plotNhoodSizeHist(obj_milo_conserv)
dev.off()

# create design
obj_design <- data.frame(colData(obj_milo_conserv))[, c("sample", "general_condition", "sample_sex", "sample_age")]  # general_condition = Celiac disease or healthy controls
obj_design <- distinct(obj_design)
rownames(obj_design) <- obj_design$sample
obj_design <- obj_design[colnames(nhoodCounts(obj_milo_conserv)), , drop=FALSE]

# create model matrix
model <- model.matrix(~ 0 + general_condition + sample_sex + sample_age, data=obj_design)

# calculate neighborhood distance and test contrast for differential abundance
obj_milo_conserv <- calcNhoodDistance(obj_milo_conserv, d=20)
contrastCvG <- testNhoods(obj_milo_conserv, design=model, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("general_conditionCELIAC - general_conditionCONTROL"))

obj_milo_conserv <- buildNhoodGraph(obj_milo_conserv)
pdf("plots/NhoodGraph_milo_kornberg_200.pdf")
plotNhoodGraphDA(obj_milo_conserv, contrastCvG, alpha=0.30) + 
  labs(caption = "k = 200, d = 20, SpatialFDR < 0.30")	
dev.off()

# annotate neighborhoods with cell type annotations
annot_df <- annotateNhoods(obj_milo_conserv, contrastCvG, "clusters")
write_csv(annot_df, paste0("milo_kornberg_results_summary_200.csv"))

# label neighborhoods without a majority cell type as "mixed"
annot_df$clusters <- ifelse(annot_df$clusters_fraction < 0.7, "Mixed", annot_df$clusters)

pdf("plots/DA_kornberg_frac_neighborhoods_hist_200_CvG.pdf")
ggplot(annot_df, aes(clusters_fraction)) + 
  geom_histogram(bins=50)
dev.off()

### Supplemental Fig 5 ----
annot_df_base <- as.data.frame(annot_df)
max_logFC <- max(abs(annot_df_base$logFC), na.rm = TRUE)

# order the clusters by their median log fold change
cluster_order <- annot_df_base %>%
  group_by(clusters) %>%
  summarize(med_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_logFC) %>%
  pull(clusters)
annot_df_base$clusters <- factor(annot_df_base$clusters, levels = cluster_order)

pdf("plots/DA_celltype_frac_neighborhoods_kornberg_200.pdf")
plotDAbeeswarm(annot_df_base, group.by = "clusters", alpha=0.30) + 
  labs(x = "Cell Type",
       y = "logFC") +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.6
  ) +
  coord_flip(ylim = c(-max_logFC, max_logFC)) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"))
dev.off()

pdf("plots/DA_results_milo_kornberg_200.pdf")
ggplot(annot_df, aes(PValue)) + geom_histogram(bins=50)
ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
  geom_point(aes(color = clusters)) +
  geom_hline(aes(yintercept = -log10(0.15), 
                 linetype = "spatialFDR < 0.15")) +  # Mark significance threshold (15% FDR)
  geom_hline(aes(yintercept = -log10(0.3), 
                 linetype = "spatialFDR < 0.30")) +
  scale_linetype_manual(name = "threshold",
                        values = c("spatialFDR < 0.15" = "solid", 
                                   "spatialFDR < 0.30" = "dashed")) +
  theme_classic() +
  labs(color = "cell type") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14))
dev.off()

nhood.dge <- findNhoodMarkers(obj_milo_conserv, contrastCvG, overlap=1, compute.new=TRUE)

contrastCvG <- groupNhoods(obj_milo_conserv, contrastCvG, max.lfc.delta = 1)
nhood.dge <- testDiffExp(x=obj_milo_conserv, da.res=contrastCvG, model.contrasts = c("general_conditionCELIAC - general_conditionCONTROL"), design=~ 0 + general_condition + sample_sex + sample_age, meta.data=obj_design)

saveRDS(nhood.dge, "milo_kornberg_200_markers.RDS")  # save object


