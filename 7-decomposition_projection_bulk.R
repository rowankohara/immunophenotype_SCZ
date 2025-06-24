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

pan_counts <- read.table(file="/local/projects/idea/rohara/immunophenotype_SCZ/GSE243841_cleaned_counts.tsv", sep="\t", header=T)


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
params@nPatterns <- 30  # trying 5, 10, 15
params@nIterations <- 20000
params@sparseOptimization <- TRUE

result <- CoGAPS(logCPM_shifted, params = params)

saveRDS(result, "/local/projects/idea/rohara/immunophenotype_SCZ/bulk_CoGAPS_30results_Pan.rds")

## visualization ----

# ow samples = SCZ, YD samples = HC
samples <- colnames(logCPM_shifted)
scz <- grepl("^ow", samples)
meta <- data.frame(sample = samples, temp = scz)
meta <- meta %>%
		mutate(group = if_else(temp == TRUE, "SCZ", "HC")) %>%
	       		select(-temp)
rownames(meta) <- samples

# extract CoGAPS results
P <- result@sampleFactors  # samples × patterns
A <- result@featureLoadings    # genes × patterns

# prepare P matrix for plotting
P_df <- as.data.frame(P)
P_df$sample <- rownames(P_df)
P_df <- merge(P_df, meta, by = "sample")

# reshape df from wide to long format for ggplot
P_long <- P_df %>%
  pivot_longer(
    cols = starts_with("Pattern"),
    names_to = "Pattern",
    values_to = "Activity"
  )

# plot
pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/CoGAPS_bulk_30patterns_boxplot.pdf")
ggplot(P_long, aes(x = group, y = Activity, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1) +
  facet_wrap(~ Pattern, scales = "free_y") +
  theme_bw() +
  labs(title = "CoGAPS Pattern Activities by Group", y = "Pattern Activity", x = "")
dev.off()

# get top 10 marker genes per pattern
markers <- patternMarkers(result, threshold = "cut")  # returns a named vector of pattern assignment per gene
top_genes <- names(markers[markers == 3])
head(sort(A[top_genes, 3], decreasing = TRUE), 10)
top_per_pattern <- lapply(seq_along(markers$PatternMarkers), function(p) {
  top_genes <- markers$PatternMarkers[[p]]
  # Sort these genes by their weight in pattern p
  top_genes <- top_genes[order(A[top_genes, p], decreasing = TRUE)]
  head(top_genes, 10)
})
names(top_per_pattern) <- paste0("Pattern_", seq_along(top_per_pattern))

# combine top marker genes from all patterns
top_genes_all <- unique(unlist(top_per_pattern))

# heatmap of A matrix (gene loadings)
pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/CoGAPS_bulk_30patterns_heatmap.pdf")
pheatmap(A[top_genes_all, ], cluster_rows = TRUE, cluster_cols = TRUE,
         scale = "row", main = "Top CoGAPS Marker Genes per Pattern")
dev.off()


