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
elist <- read.table('dkelly_singlets.pseudobulk.rds')
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

# focus on Tregs
dge_treg <- dge %>% filter(Celltype=="Treg")

# reshape to wide format: one row per gene, one column per contrast
logfc_mat <- dge_treg %>%
  select(Symbol, Contrast, logFC) %>%
  pivot_wider(names_from = Contrast, values_from = logFC) %>%
  column_to_rownames("Symbol")

# calculate variance across contrasts and filter
gene_vars <- apply(logfc_mat, 1, function(x) var(x, na.rm = TRUE))
logfc_mat_var <- logfc_mat[order(gene_vars, decreasing = TRUE), ]
logfc_mat_var <- logfc_mat_var[1:4000, ]  # choose between 3000–5000

# use kNN imputation
logfc_mat_imputed <- impute.knn(as.matrix(logfc_mat_var))$data


## make elbow plot ----
# compute total within-cluster sum of squares (WSS) for each k
set.seed(123)
wss <- numeric()
k_range <- 25:50

for (k in k_range) {
  km <- kmeans(logfc_mat_imputed, centers = k, nstart = 10)
  wss[k - 24] <- km$tot.withinss
}

# plot the elbow plot
pdf("plots/kmeans_elbow_tregs.pdf")
plot(k_range, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters K",
     ylab = "Total Within-Cluster Sum of Squares (WSS)",
     main = "Elbow Plot for K-means Clustering")
dev.off()


## k-means clustering ----
k <- 48
km <- kmeans(logfc_mat_imputed, centers = k, nstart = 25)

# save cluster assignments
modules <- data.frame(
  gene = rownames(logfc_mat_imputed),
  module = km$cluster
)


##  plot heatmap ----

# order rows by cluster
ordered_genes <- rownames(logfc_mat_imputed)[order(km$cluster)]
annotation_row <- data.frame(Module = factor(km$cluster))
rownames(annotation_row) <- rownames(logfc_mat_imputed)

# plot
pdf("plots/kmeans_heatmap_Tregs.pdf")
pheatmap(logfc_mat_imputed[ordered_genes, ],
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_row = annotation_row,
         show_rownames = FALSE,
         main = "K-means Modules Across Contrasts")
dev.off()

write.csv(modules, "gene_modules.csv", row.names = FALSE)


# K-means Clustering: expression differences between subjects ----

## elbow plots ----
for (i in 1:length(explist)){
         wss <- vector()
         max_k <- 10
         set.seed(123)
         for (k in 1:max_k) {
                km <- kmeans(explist[[i]], centers = k, nstart = 10)
                wss[k] <- km$tot.withinss
         }

         filename <- paste0('plots/kmeans_elbow_plot_', celltypes[i], '.pdf')
         pdf(filename)
         plot(1:max_k, wss, type = "b", pch = 19, frame = FALSE,
         xlab = "Number of Clusters K",
         ylab = "Total Within-Cluster Sum of Squares",
         main = "Elbow Plot for K-means")
         dev.off()
}

k <- c(2,2,2,2.2,2,4,2,4,4,5,4,5)

# initialize the list to store clustering results
all_clusters <- list()
for (i in 1:length(k)) {
        deg_expr_scaled <- explist[[i]]
        deg_genes <- dge %>% filter(PValue < 0.01 & Celltype == celltypes[i]) %>% select(Symbol) %>% pull(Symbol)
        deg_expr_scaled <- deg_expr_scaled[rownames(deg_expr_scaled) %in% deg_genes, ]
        k.use <- k[i]
        km <- kmeans(deg_expr_scaled, centers = k.use)
        clusters <- km$cluster

        # Save clustering result
        all_clusters[[celltypes[i]]] <- data.frame(
        gene = names(clusters),
        cluster = clusters)

        filename <- paste0('plots/kmeans_clustering_', celltypes[i], '.pdf')
        pdf(filename)
        pheatmap(deg_expr_scaled[order(clusters), ],
        cluster_cols = TRUE,
        cluster_rows = FALSE,
        main = paste("Cell type:", celltypes[i]),
        annotation_row = data.frame(Cluster = factor(clusters)))
        dev.off()
}
