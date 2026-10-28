# Rowan K. O'Hara-Payne
# 5/5/2025 
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Kmeans clustering of the differentially expressed genes


# Read in libraries and data ----
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# read in pseudobulk data
elist <- readRDS('dkelly_singlets.pseudobulk.rds')
str(elist)  # list of cell types with 
celltypes <- names(elist)
explist <- vector("list", length(celltypes))
names(explist) <- celltypes
explist <- lapply(elist, function(mat) {
  df <- as.data.frame(mat)
  df
})

# read in differential gene expression data
dge <- read.table('dkelly_singlets.pseudobulk_edgeR.dge.txt', sep='\t')
names(dge) <- dge[1,]
dge <- dge[-1,]

# read in meta data
sample <- colnames(elist[[11]])
MPRCID <- gsub( '\\_(.*)' , '' , sample )
MPRCID <- gsub( '\\-(.*)' , '' , MPRCID )
meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
rownames(meta) <- meta$MPRCID
meta <- meta[ MPRCID , ]


# K-means Clustering: logFC differences between contrasts ----

# store results
all_modules <- list()

# loop over cell types
celltypes <- unique(dge$Celltype)
for (ct in celltypes) {
	
	message("Processing ", ct, "...")

	# filter DGE for this cell type
	dge_ct <- dge %>% filter(Celltype == ct)

	# reshape to wide format: one row per gene, one column per contrast
	logfc_mat <- dge_ct %>%
  	select(Symbol, Contrast, logFC) %>%
  	pivot_wider(names_from = Contrast, values_from = logFC) %>%
  	column_to_rownames("Symbol")

	# calculate variance across contrasts and filter
	gene_vars <- apply(logfc_mat, 1, function(x) var(x, na.rm = TRUE))
	logfc_mat_var <- logfc_mat[order(gene_vars, decreasing = TRUE), ]
	logfc_mat_var <- logfc_mat_var[1:4000, ]  # choose between 3000–5000

	# use kNN imputation
	logfc_mat_imputed <- impute.knn(as.matrix(logfc_mat_var))$data


	# elbow plot to choose K
	set.seed(123)
  	wss <- numeric()
  	k_range <- 25:50

  	for (k in k_range) {
    		km <- kmeans(logfc_mat_imputed, centers = k, nstart = 10)
    		wss[k - min(k_range) + 1] <- km$tot.withinss
  	}

	elbow_file <- paste0("plots/kmeans_elbow_", gsub(" ", "_", ct), ".pdf")
  	pdf(elbow_file)
  	plot(k_range, wss, type = "b", pch = 19, frame = FALSE,
       		xlab = "Number of Clusters K",
       		ylab = "Total Within-Cluster Sum of Squares (WSS)",
       		main = paste("Elbow Plot for", ct))
  	dev.off()
	
	# automatically select optimal K using elbow method (second derivative)
  	diff1 <- diff(wss)
  	diff2 <- diff(diff1)
  	elbow_k_index <- which.max(-diff2) + 1  # +1 to shift from diff index to k_range index
  	k_opt <- k_range[elbow_k_index]

  	message("Selected K for ", ct, ": ", k_opt)
  	
	# run K-means clustering with selected K
  	km <- kmeans(logfc_mat_imputed, centers = k_opt, nstart = 25)

  	modules <- data.frame(
    		gene = rownames(logfc_mat_imputed),
    		module = km$cluster,
    		celltype = ct
  	)
  	all_modules[[ct]] <- modules

  	# compute module centroids
  	module_labels <- km$cluster
  	names(module_labels) <- rownames(logfc_mat_imputed)

  	module_centroids <- as.data.frame(logfc_mat_imputed) %>%
    		mutate(Module = factor(module_labels)) %>%
    		group_by(Module) %>%
    	summarise(across(everything(), mean)) %>%
    	column_to_rownames("Module")

  	# hierarchical clustering of modules
  	dist_mat <- dist(module_centroids)
  	hc_modules <- hclust(dist_mat, method = "ward.D2")
  	ordered_modules <- hc_modules$labels[hc_modules$order]

  	# heatmap of module centroids
  	heatmap_file <- paste0("plots/kmeans_heatmap_", gsub(" ", "_", ct), ".pdf")
  	pdf(heatmap_file)
  	pheatmap(as.matrix(module_centroids),
        	cluster_rows = hc_modules,
        	cluster_cols = FALSE,
        	show_rownames = TRUE,
        	show_colnames = TRUE,
        	main = paste("Hierarchical Clustering of K-means Modules (", ct, ")"),
        	color = colorRampPalette(c("blue", "white", "red"))(100))
  	dev.off()
}

# combine all modules and write to file
modules_df <- do.call(rbind, all_modules)
write.csv(modules_df, "gene_modules_all_celltypes.csv", row.names = FALSE)


# K means modules ----
dge <- read.table("dkelly_Class.pseudobulk_edgeR.dge.txt", sep = "\t", header = TRUE)

clusters_order <- c("3", "5", "4", "16", "6", "7", "31", "20", "13-0", "23", "8", "9", "17","22", "13-1","28", "2", "18", "12","24-1", "25", "0", "1", "10", "14", "27", "21-0", "21-1", "21-2", "26", "11", "15", "19", "24-0", "29", "30", "32")

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



sig.genes <- unique(dge$Symbol[dge$PValue < 0.01])
effects <- clusters_order[clusters_order %in% unique(dge$Celltype)]

nGene <- length(sig.genes)
nEffects <- length(effects)

logfc <- matrix( 0 , ncol = nEffects , nrow = nGene )
colnames(logfc) <- effects
rownames(logfc) <- sig.genes

for( i in 1:nEffects ) {
  dge.i = dge[ dge$Symbol %in% sig.genes & dge$Celltype == effects[i] , ]
  logfc[ dge.i$Symbol , i ] = dge.i$logFC
}

cl <- kmeans( logfc , centers = 10 , iter.max = 10000 , nstart = 10 )
clust <- cl$cluster
centroids <- cl$centers


# get the module gene lists for GO analysis
clusters <- split(
    rownames(logfc),
    clust     
)
enrich_list <- lapply(clusters, function(gene_set) {

    enrichGO(
        gene          = gene_set,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        universe      = dge$Symbol
    )
})
top_terms <- sapply(enrich_list, function(x) {

    if (nrow(x@result) == 0) {
        return("None")
    }

    x@result$Description[1]
})
row_annot <- data.frame(
    pathway = top_terms
)
rownames(row_annot) <- names(top_terms)
pathway_pal <- c("blood vessel morphogenesis" = "#E41A1C",
		 "post-transcriptional regulation of gene expression" = "#377EB8",
		 "cellular lipid catabolic process" = "#4DAF4A",
		 "vitamin transport" = "#984EA3",
		 "purine ribonucleotide catabolic process" = "#FF7F00",
		 "regeneration" = "#FFFF33",
		 "response to bacterium" = "#A65628",
		 "cytoplasmic translation" = "#F781BF",
		 "positive regulation of type I interferon-mediated signaling pathway" = "#999999",
		 "growth hormone receptor signaling pathway via JAK-STAT" = "#FB8072")

palette <- colorRampPalette(
    c("#2166AC", "white", "#B2182B")
)(99)
breaks <- seq(-0.5, 0.5, length.out = 101)

x <- centroids
x[ x > 0.5 ] <- 0.5
x[ x < -0.5 ] <- -0.5

labs <- parse(text = celltype_labels[colnames(x)])


ann_col <- data.frame(
  celltype = colnames(x)
)
rownames(ann_col) <- colnames(x)
ann_colors <- list(
  celltype = celltype_pal,
  pathway = pathway_pal
)

# FIGURE 3C
pdf('plots/dge.sig-genes.kmeans10.centroids.heatmap.pdf', width = 10, height = 6)
pheatmap(
  x,
  color = palette,
  breaks = breaks,
  cellwidth = 20,
  cellheight = 20,
  cluster_cols = TRUE,
  scale = "none",
  annotation_row = row_annot,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  labels_col = labs,
  annotation_legend = FALSE
)
dev.off()

MM = cor( t(logfc) , t(centroids) )


# SUPPLEMENTAL FIGURE
pdf('plots/dge.sig-genes.kmeans10.top10.heatmap.pdf', width = 10, height = 6)
for( i in 1:10 ) {
  pattern = i
  sorted = MM[ order( MM[,pattern] , decreasing =T ) , ]
  top10 = head( rownames( sorted ) , 10 )
  x = logfc[ top10 , ]
  x[ x > 0.5 ] <- 0.5
  x[ x < -0.5 ] <- -0.5
  pheatmap(
  	x,
	main = paste( 'Module' , i ),
  	color = palette,
  	breaks = breaks,
	cellwidth = 20,
	cellheight = 20,
	cluster_rows = FALSE,
  	cluster_cols = TRUE,
  	scale = "none",
  	annotation_col = ann_col,
  	annotation_colors = ann_colors,
  	labels_col = labs,
  	annotation_legend = FALSE)
}
dev.off()


colnames(MM) <- paste( 'kM' , 1:10 , sep = '' )
df <- data.frame(
  gene = sig.genes ,
  clust ,
  MM ,
  logfc )

write.table( df , quote=F , sep='\t' , row.names=F ,
             file = 'kmeans.k10.module_membership.txt' )

saveRDS( cl , file = 'kmeans.k10.rds' )


celltype_map <- c(
  "X0"    = "0",
  "X1"    = "1",
  "X2"    = "2",
  "X3"    = "3",
  "X4"    = "4",
  "X5"    = "5",
  "X6"    = "6",
  "X7"    = "7",
  "X8"    = "8",
  "X9"    = "9",
  "X10"   = "10",
  "X11"   = "11",
  "X12"   = "12",
  "X13.0" = "13-0",
  "X13.1" = "13-1",
  "X14"   = "14",
  "X15"   = "15",
  "X16"   = "16",
  "X17"   = "17",
  "X18"   = "18",
  "X19"   = "19",
  "X20"   = "20",
  "X21.0" = "21-0",
  "X21.1" = "21-1",
  "X21.2" = "21-2",
  "X22"   = "22",
  "X23"   = "23",
  "X24.0" = "24-0",
  "X24.1" = "24-1",
  "X25"   = "25",
  "X26"   = "26",
  "X27"   = "27",
  "X28"   = "28",
  "X29"   = "29",
  "X30"   = "30",
  "X31"   = "31",
  "X32"   = "32"
)
r <- cor(df[, 13:37])

x <- r
x[x > 0.5] <- 0.5
x[x < -0.5] <- -0.5

celltypes <- celltype_map[colnames(x)]

ann_col <- data.frame(celltype = celltypes)
rownames(ann_col) <- colnames(x)

pdf('plots/corplot.pdf')

pheatmap(
  x,
  color = palette,
  breaks = breaks,
  cluster_cols = TRUE,
  scale = "none",
  annotation_row = ann_col,
  annotation_col = ann_col,
  annotation_colors = list(celltype = celltype_pal),
  labels_col = labs,
  labels_row = labs,
  annotation_legend = FALSE
)

dev.off()


