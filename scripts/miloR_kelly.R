library(miloR)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(patchwork)

obj <- readRDS("kelly_postqc_clusters_meta.rds")

obj_sce <- as.SingleCellExperiment(obj)
obj_milo <- Milo(obj_sce)

# get the UMAP from the single cell experiment object
reducedDim(obj_milo, "UMAP") <- reducedDim(obj_sce, "UMAP")
reducedDim(obj_milo, "PCA") <- Embeddings(obj, "pca")[, 1:30]

obj_milo_conserv <- buildGraph(obj_milo, k = 45, d = 20)
obj_milo_conserv <- makeNhoods(obj_milo_conserv, prop = 0.2, k = 45, d = 20, refined = TRUE, refinement_scheme="graph")
obj_milo_conserv <- countCells(obj_milo_conserv, meta.data = data.frame(colData(obj_milo_conserv)), samples = "sample")

pdf("plots/NHoodSize_milo_kelly_45.pdf")
p1 <- plotNhoodSizeHist(obj_milo_conserv)
print(p1)
dev.off()

obj_milo_conserv <- calcNhoodDistance(obj_milo_conserv, d=20)
saveRDS(obj_milo_conserv, "miloR_conserv_obj_kelly.RDS")

## SCZHigh vs SCZLow Contrast ----

obj_design <- data.frame(SingleCellExperiment::colData(obj_milo_conserv))[, c("sample", "Group", "Sex", "Age")]
obj_design <- distinct(obj_design)
rownames(obj_design) <- obj_design$sample
obj_design <- obj_design[colnames(nhoodCounts(obj_milo_conserv)), , drop=FALSE]
table(obj_design$Group)

obj_milo_HvL <- obj_milo_conserv

model <- model.matrix(~ 0 + Group + Sex + Age, data=obj_design)
contrastHvL <- testNhoods(obj_milo_HvL, design=model, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("GroupSCZHigh - GroupSCZLow"))

obj_milo_HvL <- buildNhoodGraph(obj_milo_HvL)
pdf("plots/NhoodGraph_milo_kelly_HvL.pdf", width = 12, height = 8)
p2 <- plotNhoodGraphDA(obj_milo_HvL, contrastHvL, alpha=0.30) + 
	scale_fill_gradient2(low = "#2166AC",
                              mid = "white",
                              high = "#B2182B") +
	labs(caption = "k = 45, d = 20, SpatialFDR < 0.30")
p3 <- plotNhoodGraphDA(obj_milo_HvL, contrastHvL, alpha=0.15) +
	scale_fill_gradient2(low = "#2166AC",
                              mid = "white",
                              high = "#B2182B") +	
	labs(caption = "k = 45, d = 20, SpatialFDR < 0.15")
print(p2 + labs(title = "Differential Abundance of Cell Type Neighborhoods Between AGA IgG+ and AGA IgG- Schizophrenia"))
print(p3)
dev.off()

annot_df <- annotateNhoods(obj_milo_HvL, contrastHvL, "seurat_clusters")
write_csv(annot_df, "milo_kelly_results_summary_HvL.csv")

annot_df$seurat_clusters <- ifelse(annot_df$seurat_clusters_fraction < 0.7, "Mixed", annot_df$seurat_clusters)

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


cluster_order <- annot_df %>%
  group_by(seurat_clusters) %>%
  summarize(med_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_logFC) %>%
  pull(seurat_clusters)

annot_df$seurat_clusters <- factor(annot_df$seurat_clusters, levels = cluster_order)

pdf("plots/DA_kelly_frac_neighborhoods_45_HvL.pdf", width = 12, height = 8)
p4 <- ggplot(annot_df, aes(seurat_clusters_fraction)) + 
	geom_histogram(bins=50)
p5 <- plotDAbeeswarm(annot_df, group.by = "seurat_clusters", alpha=0.30) + 
	labs(x = "Cell Type",
	     y = "logFC") +
	scale_color_gradient2(low = "#2166AC",
			      mid = "white",
			      high = "#B2182B") +
	geom_hline(yintercept = 0,
		   linetype = "dashed",
		   color = "black",
		   linewidth = 1) +
	coord_flip(ylim = c(-max_logFC, max_logFC)) + 
	scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) + 
	theme(axis.text = element_text(size = 12),
	      axis.title = element_text(size = 14))
p6 <- plotDAbeeswarm(annot_df, group.by = "seurat_clusters", alpha=0.15) + 
	labs(x = "Cell Type",
	     y = "logFC") + 
	scale_color_gradient2(low = "#2166AC",
                              mid = "white",
                              high = "#B2182B") +
	geom_hline(yintercept = 0,
                   linetype = "dashed",
                   color = "black",
                   linewidth = 1) +
        coord_flip(ylim = c(-max_logFC, max_logFC)) +
	scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) + 
	theme(axis.text = element_text(size = 12),
	      axis.title = element_text(size = 14))
print(p4)
print(p5)
print(p6)
dev.off()

annot_df$seurat_clusters <- factor(
  annot_df$seurat_clusters,
  levels = names(celltype_pal)
)

pdf("plots/DA_results_milo_kelly_HvL.pdf", width = 12, height = 8)
p7 <- ggplot(annot_df, aes(PValue)) + 
	geom_histogram(bins=50)
p8 <- ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
  	geom_point(aes(color = seurat_clusters)) +
  	scale_color_manual(labels = function(x) parse(text = celltype_labels[x]), values = celltype_pal) +
  	geom_hline(aes(yintercept = -log10(0.15), linetype = "spatialFDR < 0.15")) + ## Mark significance threshold (15% FDR)
  	geom_hline(aes(yintercept = -log10(0.3), linetype = "spatialFDR < 0.30")) +
	scale_linetype_manual(name = "threshold", values = c("spatialFDR < 0.15" = "solid", "spatialFDR < 0.30" = "dashed")) +
	theme_classic() +
	labs(color = "cell type") +
	theme(legend.text = element_text(size = 10),
	      legend.title = element_text(size = 14),
	      axis.title = element_text(size = 14))
print(p7)
print(p8)
dev.off()


pdf("plots/DA_results_milo_kelly_combined_HvL.pdf", width = 16, height = 21)
print(
      p2 / (p5 + p8) + 
	      plot_annotation(title = "Differential Abundance of Cell Type Neighborhoods Between AGA IgG+ and AGA IgG- Schizophrenia", 
			      theme = theme(plot.title = element_text(hjust = 0.5,
								      size =20)))
      )
dev.off()


pdf("plots/DA_results_milo_kelly_combined2_HvL.pdf", width = 16, height = 9)
print(
      p5 + p8
      )
dev.off()


nhood.dge <- findNhoodMarkers(obj_milo_HvL, contrastHvL, overlap=1, compute.new=TRUE)

saveRDS(nhood.dge, "milo_kelly_45_markers_HvL.RDS")


## SCZ vs Healthy Controls Contrast ----

obj_milo_CvS <- obj_milo_conserv

obj_design_CvS <- data.frame(colData(obj_milo_CvS))[, c("sample", "Condition", "Sex", "Age")]
obj_design_CvS <- distinct(obj_design_CvS)
rownames(obj_design_CvS) <- obj_design_CvS$sample
obj_design_CvS <- obj_design_CvS[colnames(nhoodCounts(obj_milo_CvS)), , drop=FALSE]
table(obj_design_CvS$Condition)


model_CvS <- model.matrix(~ 0 + Condition + Sex + Age, data=obj_design_CvS)

contrastCvS <- testNhoods(obj_milo_CvS, design=model_CvS, design.df=obj_design_CvS, fdr.weighting="graph-overlap", model.contrasts = c("ConditionSCZ - ConditionHC"))

obj_milo_CvS <- buildNhoodGraph(obj_milo_CvS)
pdf("plots/NhoodGraph_milo_kelly_CvS.pdf", width = 12, height = 8)
p2 <- plotNhoodGraphDA(obj_milo_CvS, contrastCvS, alpha=0.30) +
        labs(caption = "k = 45, d = 20, SpatialFDR < 0.30")
p3 <- plotNhoodGraphDA(obj_milo_CvS, contrastCvS, alpha=0.15) + labs(caption = "k = 45, d = 20, SpatialFDR < 0.15")
print(p2 + labs(title = "Differential Abundance of Cell Type Neighborhoods Between Controls and Schizophrenia"))
print(p3)
dev.off()

annot_df <- annotateNhoods(obj_milo_CvS, contrastCvS, "seurat_clusters")
write_csv(annot_df, "milo_kelly_results_summary_CvS.csv")

annot_df$seurat_clusters <- ifelse(annot_df$seurat_clusters_fraction < 0.7, "Mixed", annot_df$seurat_clusters)

cluster_order <- annot_df %>%
  group_by(seurat_clusters) %>%
  summarize(med_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_logFC) %>%
  pull(seurat_clusters)

annot_df$seurat_clusters <- factor(annot_df$seurat_clusters, levels = cluster_order)

library(ggbeeswarm)
pdf("plots/DA_kelly_frac_neighborhoods_45_CvS.pdf", width = 12, height = 8)
p4 <- ggplot(annot_df, aes(seurat_clusters_fraction)) +
        geom_histogram(bins=50)
p5 <- ggplot(annot_df, aes(x = seurat_clusters, y = logFC)) +
	geom_quasirandom(color = "gray45", width = 0.25) +
	labs(x = "Cell Type",
	     y = "logFC") +
        scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) +
	coord_flip() +
	theme_bw() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
print(p4)
print(p5)
dev.off()

annot_df$seurat_clusters <- factor(
  annot_df$seurat_clusters,
  levels = names(celltype_pal)
)

pdf("plots/DA_results_milo_kelly_CvS.pdf", width = 12, height = 8)
p7 <- ggplot(annot_df, aes(PValue)) +
        geom_histogram(bins=50)
p8 <- ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
        geom_point(aes(color = seurat_clusters)) +
        scale_color_manual(labels = function(x) parse(text = celltype_labels[x]), values = celltype_pal) +
        geom_hline(aes(yintercept = -log10(0.15), linetype = "spatialFDR < 0.15")) + ## Mark significance threshold (15% FDR)
        geom_hline(aes(yintercept = -log10(0.3), linetype = "spatialFDR < 0.30")) +
        scale_linetype_manual(name = "threshold", values = c("spatialFDR < 0.15" = "solid", "spatialFDR < 0.30" = "dashed")) +
	theme_classic() +
        labs(color = "cell type") +
        theme(legend.text = element_text(size = 10),
              legend.title = element_text(size = 14),
              axis.title = element_text(size = 14))
print(p7)
print(p8)
dev.off()


pdf("plots/DA_results_milo_kelly_combined_CvS.pdf", width = 16, height = 21)
print(
      p2 / (p5 + p8) +
              plot_annotation(title = "Differential Abundance of Cell Type Neighborhoods Between Healthy Controls and Schizophrenia",
                              theme = theme(plot.title = element_text(hjust = 0.5,
                                                                      size =20)))
      )
dev.off()


pdf("plots/DA_results_milo_kelly_combined2_CvS.pdf", width = 16, height = 9)
print(
      p5 + p8
      )
dev.off()


#nhood.dge <- findNhoodMarkers(obj_milo_CvS, contrastCvS, overlap=1, compute.new=TRUE)

#saveRDS(nhood.dge, "milo_kelly_45_markers_CvS.RDS")



## SCZHigh vs Control ----

obj_milo_HvC <- obj_milo_conserv

contrastHvC <- testNhoods(obj_milo_HvC, design=model, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("GroupSCZHigh - GroupHC"))

obj_milo_HvC <- buildNhoodGraph(obj_milo_HvC)
pdf("plots/NhoodGraph_milo_kelly_HvC.pdf", width = 12, height = 8)
p2 <- plotNhoodGraphDA(obj_milo_HvC, contrastHvC, alpha=0.30) +
        labs(caption = "k = 45, d = 20, SpatialFDR < 0.30")
p3 <- plotNhoodGraphDA(obj_milo_HvC, contrastHvC, alpha=0.15) + labs(caption = "k = 45, d = 20, SpatialFDR < 0.15")
print(p2 + labs(title = "Differential Abundance of Cell Type Neighborhoods Between AGA IgG+ Schizophrenia and Controls"))
print(p3)
dev.off()

annot_df <- annotateNhoods(obj_milo_HvC, contrastHvC, "seurat_clusters")
write_csv(annot_df, "milo_kelly_results_summary_HvC.csv")

annot_df$seurat_clusters <- ifelse(annot_df$seurat_clusters_fraction < 0.7, "Mixed", annot_df$seurat_clusters)

cluster_order <- annot_df %>%
  group_by(seurat_clusters) %>%
  summarize(med_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_logFC) %>%
  pull(seurat_clusters)

annot_df$seurat_clusters <- factor(annot_df$seurat_clusters, levels = cluster_order)

pdf("plots/DA_kelly_frac_neighborhoods_45_HvC.pdf", width = 12, height = 8)
p4 <- ggplot(annot_df, aes(seurat_clusters_fraction)) +
        geom_histogram(bins=50)
p5 <- plotDAbeeswarm(annot_df, group.by = "seurat_clusters", alpha=0.30) +
        labs(x = "Cell Type",
	     y = "logFC") +
        scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
p6 <- plotDAbeeswarm(annot_df, group.by = "seurat_clusters", alpha=0.15) +
        labs(x = "Cell Type",
	     y = "logFC") +
        scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
print(p4)
print(p5)
print(p6)
dev.off()

annot_df$seurat_clusters <- factor(
  annot_df$seurat_clusters,
  levels = names(celltype_pal)
)


pdf("plots/DA_results_milo_kelly_HvC.pdf", width = 12, height = 8)
p7 <- ggplot(annot_df, aes(PValue)) +
        geom_histogram(bins=50)
p8 <- ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
        geom_point(aes(color = seurat_clusters)) +
        scale_color_manual(labels = function(x) parse(text = celltype_labels[x]), values = celltype_pal) +
        geom_hline(aes(yintercept = -log10(0.15), linetype = "spatialFDR < 0.15")) + ## Mark significance threshold (15% FDR)
        geom_hline(aes(yintercept = -log10(0.3), linetype = "spatialFDR < 0.30")) +
	scale_linetype_manual(name = "threshold", values = c("spatialFDR < 0.15" = "solid", "spatialFDR < 0.30" = "dashed")) +
        theme_classic() +
        labs(color = "cell type") +
        theme(legend.text = element_text(size = 10),
              legend.title = element_text(size = 14),
              axis.title = element_text(size = 14))
print(p7)
print(p8)
dev.off()


pdf("plots/DA_results_milo_kelly_combined_HvC.pdf", width = 16, height = 21)
print(
      p2 / (p5 + p8) +
              plot_annotation(title = "Differential Abundance of Cell Type Neighborhoods Between AGA IgG+ Schizophrenia and Healthy Controls",
                              theme = theme(plot.title = element_text(hjust = 0.5,
                                                                      size =20)))
      )
dev.off()


pdf("plots/DA_results_milo_kelly_combined2_HvC.pdf", width = 16, height = 9)
print(
      p5 + p8
      )
dev.off()


nhood.dge <- findNhoodMarkers(obj_milo_HvC, contrastHvC, overlap=1, compute.new=TRUE)

saveRDS(nhood.dge, "milo_kelly_45_markers_HvC.RDS")


## SCZLow vs Healthy Controls ----

obj_milo_LvC <- obj_milo_conserv

contrastLvC <- testNhoods(obj_milo_LvC, design=model, design.df=obj_design, fdr.weighting="graph-overlap", model.contrasts = c("GroupSCZLow - GroupHC"))

obj_milo_LvC <- buildNhoodGraph(obj_milo_LvC)
pdf("plots/NhoodGraph_milo_kelly_LvC.pdf", width = 12, height = 8)
p2 <- plotNhoodGraphDA(obj_milo_LvC, contrastLvC, alpha=0.30) +
        labs(caption = "k = 45, d = 20, SpatialFDR < 0.30")
p3 <- plotNhoodGraphDA(obj_milo_LvC, contrastLvC, alpha=0.15) + labs(caption = "k = 45, d = 20, SpatialFDR < 0.15")
print(p2 + labs(title = "Differential Abundance of Cell Type Neighborhoods Between AGA IgG- Schizophrenia and Controls"))
print(p3)
dev.off()

annot_df <- annotateNhoods(obj_milo_LvC, contrastLvC, "seurat_clusters")
write_csv(annot_df, "milo_kelly_results_summary_LvC.csv")

annot_df$seurat_clusters <- ifelse(annot_df$seurat_clusters_fraction < 0.7, "Mixed", annot_df$seurat_clusters)

cluster_order <- annot_df %>%
  group_by(seurat_clusters) %>%
  summarize(med_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_logFC) %>%
  pull(seurat_clusters)

annot_df$seurat_clusters <- factor(annot_df$seurat_clusters, levels = cluster_order)



pdf("plots/DA_kelly_frac_neighborhoods_45_LvC.pdf", width = 12, height = 8)
p4 <- ggplot(annot_df, aes(seurat_clusters_fraction)) +
        geom_histogram(bins=50)
p5 <- ggplot(annot_df, aes(x = seurat_clusters, y = logFC)) +
        geom_quasirandom(color = "gray45", width = 0.25) +
        labs(x = "Cell Type",
             y = "logFC") +
        scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) +
        coord_flip() +
        theme_bw() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
print(p4)
print(p5)
dev.off()


annot_df$seurat_clusters <- factor(
  annot_df$seurat_clusters,
  levels = names(celltype_pal)
)

pdf("plots/DA_results_milo_kelly_LvC.pdf", width = 12, height = 8)
p7 <- ggplot(annot_df, aes(PValue)) +
        geom_histogram(bins=50)
p8 <- ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
        geom_point(aes(color = seurat_clusters)) +
        scale_color_manual(labels = function(x) parse(text = celltype_labels[x]), values = celltype_pal) +
        geom_hline(aes(yintercept = -log10(0.15), linetype = "spatialFDR < 0.15")) + ## Mark significance threshold (15% FDR)
        geom_hline(aes(yintercept = -log10(0.3), linetype = "spatialFDR < 0.30")) +
	scale_linetype_manual(name = "threshold", values = c("spatialFDR < 0.15" = "solid", "spatialFDR < 0.30" = "dashed")) +
        theme_classic() +
        labs(color = "cell type") +
        theme(legend.text = element_text(size = 10),
              legend.title = element_text(size = 14),
              axis.title = element_text(size = 14))
print(p7)
print(p8)
dev.off()


pdf("plots/DA_results_milo_kelly_combined_LvC.pdf", width = 16, height = 21)
print(
      p2 / (p5 + p8) +
              plot_annotation(title = "Differential Abundance of Cell Type Neighborhoods Between AGA IgG- Schizophrenia and Healthy Controls",
                              theme = theme(plot.title = element_text(hjust = 0.5,
                                                                      size =20)))
      )
dev.off()


pdf("plots/DA_results_milo_kelly_combined2_LvC.pdf", width = 16, height = 9)
print(
      p5 + p8
      )
dev.off()


#nhood.dge <- findNhoodMarkers(obj_milo_LvC, contrastLvC, overlap=1, compute.new=TRUE)

#saveRDS(nhood.dge, "milo_kelly_45_markers_LvC.RDS")



## SCZ High vs Controls & SCZ Low ----
obj_milo_HvSC <- obj_milo_conserv

obj_design_HvSC <- data.frame(colData(obj_milo_HvSC))[, c("sample", "Class", "Sex", "Age")]
obj_design_HvSC <- distinct(obj_design_HvSC)
rownames(obj_design_HvSC) <- obj_design_HvSC$sample
obj_design_HvSC <- obj_design_HvSC[colnames(nhoodCounts(obj_milo_HvSC)), , drop=FALSE]
table(obj_design_HvSC$Class)


model_HvSC <- model.matrix(~ 0 + Class + Sex + Age, data=obj_design_HvSC)

contrastHvSC <- testNhoods(obj_milo_HvSC, design=model_HvSC, design.df=obj_design_HvSC, fdr.weighting="graph-overlap", model.contrasts = c("ClassPOS - ClassNEG"))

obj_milo_HvSC <- buildNhoodGraph(obj_milo_HvSC)
pdf("plots/NhoodGraph_milo_kelly_HvSC.pdf", width = 12, height = 8)
p2 <- plotNhoodGraphDA(obj_milo_HvSC, contrastHvSC, alpha=0.30) +
        labs(caption = "k = 45, d = 20, SpatialFDR < 0.30")
p3 <- plotNhoodGraphDA(obj_milo_HvSC, contrastHvSC, alpha=0.15) + labs(caption = "k = 45, d = 20, SpatialFDR < 0.15")
print(p2 + labs(title = "Differential Abundance of Cell Type Neighborhoods Between Controls and Schizophrenia"))
print(p3)
dev.off()

annot_df <- annotateNhoods(obj_milo_HvSC, contrastHvSC, "seurat_clusters")
write_csv(annot_df, "milo_kelly_results_summary_HvSC.csv")

annot_df$seurat_clusters <- ifelse(annot_df$seurat_clusters_fraction < 0.7, "Mixed", annot_df$seurat_clusters)

cluster_order <- annot_df %>%
  group_by(seurat_clusters) %>%
  summarize(med_logFC = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(med_logFC) %>%
  pull(seurat_clusters)

annot_df$seurat_clusters <- factor(annot_df$seurat_clusters, levels = cluster_order)


pdf("plots/DA_kelly_frac_neighborhoods_45_HvSC.pdf", width = 12, height = 8)
p4 <- ggplot(annot_df, aes(seurat_clusters_fraction)) +
        geom_histogram(bins=50)
p5 <- plotDAbeeswarm(annot_df, group.by = "seurat_clusters", alpha=0.30) +
        labs(x = "Cell Type",            
             y = "logFC") +    
        scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
p6 <- plotDAbeeswarm(annot_df, group.by = "seurat_clusters", alpha=0.15) +
        labs(x = "Cell Type",
             y = "logFC") +
        scale_x_discrete(labels = function(x) parse(text = celltype_labels[x])) +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
print(p4)
print(p5)
print(p6)
dev.off()

annot_df$seurat_clusters <- factor(
  annot_df$seurat_clusters,
  levels = names(celltype_pal)
)

pdf("plots/DA_results_milo_kelly_HvSC.pdf", width = 12, height = 8)
p7 <- ggplot(annot_df, aes(PValue)) +
        geom_histogram(bins=50)
p8 <- ggplot(annot_df, aes(logFC, -log10(SpatialFDR))) +
        geom_point(aes(color = seurat_clusters)) +
        scale_color_manual(labels = function(x) parse(text = celltype_labels[x]), values = celltype_pal) +
        geom_hline(aes(yintercept = -log10(0.15), linetype = "spatialFDR < 0.15")) + ## Mark significance threshold (15% FDR)
        geom_hline(aes(yintercept = -log10(0.3), linetype = "spatialFDR < 0.30")) +
        scale_linetype_manual(name = "threshold", values = c("spatialFDR < 0.15" = "solid", "spatialFDR < 0.30" = "dashed")) +
        theme_classic() +
        labs(color = "cell type") +
        scale_y_continuous(limits = c(0, NA),
                           expand = expansion(mult = c(0, 0.05))
                             ) +
        theme(plot.title = element_text(hjust = 0.5,
                                        size = 20),
              axis.text = element_text(size = 15),
              axis.title = element_text(size = 17),
              axis.line = element_line(linewidth = 1.25),
              axis.ticks = element_line(linewidth = 1.25))
print(p7)
print(p8)
dev.off()


pdf("plots/DA_results_milo_kelly_combined_HvSC.pdf", width = 16, height = 21)
print(
      p2 / (p5 + p8) +
              plot_annotation(title = "Differential Abundance of Cell Type Neighborhoods Between AGA IgG+ Schizophrenia and AGA IgG- Controls",
                              theme = theme(plot.title = element_text(hjust = 0.5,
                                                                      size =20)))
      )
dev.off()


pdf("plots/DA_results_milo_kelly_combined2_HvSC.pdf", width = 16, height = 9)
print(
      p5 + p8
      )
dev.off()


nhood.dge <- findNhoodMarkers(obj_milo_HvSC, contrastHvSC, overlap=1, compute.new=TRUE)

saveRDS(nhood.dge, "milo_kelly_45_markers_HvSC.RDS")



## testing ----

contrastHvL <- groupNhoods(obj_milo_conserv, contrastHvL, max.lfc.delta = 1)

all(colnames(obj_milo_conserv) == rownames(colData(obj_milo_conserv)))

# Check set equality
setequal(colnames(obj_milo_conserv), rownames(colData(obj_milo_conserv)))

# Find mismatches
setdiff(colnames(obj_milo_conserv), rownames(colData(obj_milo_conserv)))
setdiff(rownames(colData(obj_milo_conserv)), colnames(obj_milo_conserv))

any(duplicated(colnames(obj_milo_conserv)))
any(duplicated(rownames(colData(obj_milo_conserv))))

table(colData(obj_milo_conserv)$Group)

nhood.dge <- testDiffExp(x=obj_milo_conserv, da.res=contrastHvL, design=~ 0 + Group + Sex + Age, meta.data=obj_design)

saveRDS(nhood.dge, "milo_kelly_45_DGE.RDS")  # save object
