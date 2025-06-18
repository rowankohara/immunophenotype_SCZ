# Rowan K. O'Hara-Payne
# 6/17/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Decomposition of the Pan et. al data and projection onto our data

# Read in libraries and data ----
library(tidyverse)
library(Seurat)
library(edgeR)
library(CoGAPS)
library(projectR)
library(fgsea)
library(msigdbr)
library(pheatmap)

pan_counts <- read.table(file="GSE243841_cleaned_counts.tsv", sep="\t", header=T)


## data wrangling ----

# check if Pan data is normalized
rownames(pan_counts) <- pan_counts$Geneid  # set rownames to gene ids
pan_counts <- pan_counts[ , -1]  # remove gene id column (matrix can be only one data type)
pan_matrix <- as.matrix(pan_counts)
summary(as.numeric(pan_matrix))  # large range, not normalized

# use edgeR to normalize expression counts
dge <- DGEList(counts = pan_matrix)  # no need to specify groups just for normalization

# filter lowly expressed genes (optional but recommended)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# TMM normalization
dge <- calcNormFactors(dge)

# get logCPM values (recommended for CoGAPS)
logCPM <- cpm(dge, log=TRUE, prior.count=1)


# CoGAPS ----

# log normalization introduced negative expression values but CoGAPS uses NMF
min_val <- min(logCPM)
logCPM_shifted <- logCPM - min_val  # now minimum value is 0

params <- new("CogapsParams")
params@nPatterns <- 5  # trying 5, 10, 15
params@nIterations <- 20000
params@sparseOptimization <- TRUE

result <- CoGAPS(logCPM_shifted, params = params)


## visualization ----

# ow samples = SCZ, YD samples = HC
samples <- colnames(logCPM_shifted)
scz <- grep("^ow", samples)
meta <- data.frame(sample = samples, temp = scz)
meta <- meta %>%
		mutate(group = if_else(temp == T, "SCZ", "HC") %>%
	       		select(-temp)
rownames(meta) <- samples

# extract CoGAPS results
A <- amplitudeMatrix(result)  # genes × patterns
P <- patternMatrix(result)    # patterns × samples

# prepare P matrix for plotting
P_df <- as.data.frame(t(P))  # samples × patterns
P_df$Sample <- rownames(P_df)
P_df <- merge(P_df, meta, by.x = "Sample", by.y = "sample")

# reshape df from wide to long format for ggplot
P_long <- P_df %>%
  pivot_longer(
    cols = starts_with("Pattern"),
    names_to = "Pattern",
    values_to = "Activity"
  )

# plot
pdf("./plots/CoGAPS_bulk_5patterns_boxplot.pdf")
ggplot(P_melt, aes(x = Group, y = Activity, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~ Pattern, scales = "free_y") +
  theme_bw() +
  labs(title = "CoGAPS Pattern Activities by Group", y = "Pattern Activity", x = "")
dev.off()

# get top 10 marker genes per pattern
top_per_pattern <- lapply(unique(markers), function(p) {
  g <- names(which(markers == p))
  g[order(A[g, p], decreasing = TRUE)][1:10]
})
top_genes_all <- unique(unlist(top_per_pattern))

# heatmap of A matrix (gene loadings)
pdf("./plots/CoGAPS_bulk_5patterns_heatmap.pdf")
pheatmap(A[top_genes_all, ], cluster_rows = TRUE, cluster_cols = TRUE,
         scale = "row", main = "Top CoGAPS Marker Genes per Pattern")
dev.off()


