# Rowan K. O'Hara-Payne
# 5/5/2025 
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Kmeans clustering of the differentially expressed genes


# Read in libraries and data ----
library(dplyr)
library(tibble)
library(tidyr)
library(impute)
library(pheatmap)

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
