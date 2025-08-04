# Rowan K. O'Hara-Payne
# 7/12/2025
# conda activate /autofs/projects-t3/idea/amentlab_software/conda_envs/Seurat5_R4_3_0
# Decomposition of our data to be projected onto the bulk data Pan et. al and Glahn et. al


# Read in libraries and data ----
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

library(CoGAPS)
library(Seurat)


obj <- readRDS(input_file)


# Find variable genes ----

obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
hvg <- VariableFeatures(obj)

expr_mat <- GetAssayData(obj[["RNA"]], layer = "data")  # log-normalized
expr_mat <- expr_mat[hvg, ]                          # subset to HVGs
expr_mat <- as.matrix(expr_mat)

# shift if needed
if (min(expr_mat) < 0) {
  expr_mat <- expr_mat - min(expr_mat) + 1e-6
}


# Run CoGAPS ----

nThrd=16
cells <- colnames(expr_mat)
genes <- rownames(expr_mat)
params <- new("CogapsParams")
params <- setParam(params, "sparseOptimization", FALSE)
params <- setParam(params, "nIterations", 10000)
params <- setParam(params, "nPatterns", 30)
params <- setDistributedParams(params, nSets=1)

xxCoGAPS <- CoGAPS(data=expr_mat, params=params, geneNames=genes, sampleNames=cells, messages=TRUE, transposeData=FALSE, nPatterns=30, nIterations=10000, nThreads=16)

output_file <- paste0("/autofs/projects-t2/idea/rohara/immunophenotype_SCZ/", output_file, "CoGAPS_30Pat.rds")

saveRDS(xxCoGAPS, file = output_file)


# Interpret CoGAPS results ----

library(tidyverse)
library(pheatmap)

# Kelly data wrangling
results_kelly <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kellyCoGAPS_30Pat.rds")
P_kelly <- results_kelly@sampleFactors
A_kelly <- results_kelly@featureLoadings

kelly_counts <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/dkelly_singlets.postqc.azimuth_celltypes.rds")
kelly_meta <- as.data.frame(kelly_counts@meta.data)
kelly_meta$cell_id <- rownames(kelly_meta)
MPRCID <- gsub( '\\_(.*)' , '' , kelly_meta$sample )
meta <- read.csv("/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv")
sample_ids <- unique(kelly_meta$sample)
# match each MPRCID to a sample ID that is a prefix of it
meta$sample <- sapply(meta$MPRCID, function(id) {
  match <- sample_ids[startsWith(sample_ids, id)]
  if (length(match) > 0) match[1] else NA
})
kelly_meta <- left_join(kelly_meta, meta, by = "sample")
rownames(kelly_meta) <- kelly_meta$cell_id
kelly <- merge(P_kelly, kelly_meta, by = "row.names")
rownames(kelly) <- kelly$cell_id
kelly <- kelly[ , -1]

# Kornberg data wrangling
results_kornberg <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kornbergCoGAPS_30Pat.rds")
P_kornberg <- results_kornberg@sampleFactors
A_kornberg <- results_kornberg@featureLoadings

kornberg_counts <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/GSE252762_kornberg_han_celiac_dataset_fig1-2.RDS")
kornberg_counts[['percent.mt']] <- PercentageFeatureSet(kornberg_counts, pattern = '^MT-')
kornberg_counts <- subset(kornberg_counts, nCount_RNA < 10000 & nCount_RNA > 500 & percent.mt < 15)
predictions <- read.delim("/local/projects/idea/rohara/immunophenotype_SCZ/azimuth_pred.tsv", row.names = 1)
kornberg_counts <- AddMetaData(
        object = kornberg_counts,
        metadata = predictions)
kornberg <- merge(P_kornberg, kornberg_counts@meta.data, by = "row.names")
rownames(kornberg) <- kornberg$Row.names
kornberg <- kornberg[ , -1]

# linear modeling
lm_results <- lapply(colnames(kelly)[c(1:30)], function(pat) {
  formula <- as.formula(paste(pat, "~ Group + Sex + Age + predicted.celltype.l2"))
  lm(formula, data = kelly)
})

names(lm_results) <- colnames(kelly)[c(1:30)]

extract_lm_results <- function(lm_list) {
  # lm_list: a named or unnamed list of lm objects
  do.call(rbind, lapply(seq_along(lm_list), function(i) {
    lm_model <- lm_list[[i]]
    coef_summary <- summary(lm_model)$coefficients
    df <- as.data.frame(coef_summary)
    df$term <- rownames(df)
    df$model <- if (!is.null(names(lm_list))) names(lm_list)[i] else paste0("model_", i)
    df <- df[, c("model", "term", "Estimate", "Std. Error", "t value", "Pr(>|t|)")]
    rownames(df) <- NULL
    return(df)
  }))
}

results_df <- extract_lm_results(lm_results)

write.csv(results_df, "/local/projects/idea/rohara/immunophenotype_SCZ/kornberg_proj_model_results.csv", row.names = FALSE)

# plot
P <- kelly[ , c(1:30, 62)]
P_long <- P %>%
  pivot_longer(
    cols = starts_with("Pattern"),
    names_to = "Pattern",
    values_to = "Activity"
  )

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/CoGAPS_kornberg_30patterns_boxplot.pdf")
ggplot(P_long, aes(x = Group, y = Activity, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ Pattern, scales = "free_y") +
  theme_bw() +
  labs(title = "Kornberg CoGAPS Pattern Activities by Group", y = "Pattern Activity", x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

beta_df <- data.frame(
  Pattern = colnames(kelly)[c(1:30)],
  Beta = sapply(lm_results, function(x) coef(summary(x))["general_conditionCONTROL", "Estimate"]),
  P = sapply(lm_results, function(x) coef(summary(x))["general_conditionCONTROL", "Pr(>|t|)"])
)
beta_df$Pattern <- factor(beta_df$Pattern, levels = beta_df$Pattern)

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kelly_pattern_betas.pdf")
ggplot(beta_df, aes(x = Pattern, y = Beta, fill = P < 0.05)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values = c("FALSE" = "#999999", "TRUE" = "#f9421b")) +
  geom_text(aes(label = sprintf("%.2g", P)), vjust = -0.5, size = 3) +
  labs(title = "Kelly Pattern Betas") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

# Extract betas and p-values for each group contrast
contrasts <- c("GroupSCZ Low", "GroupSCZ High")

# Build a long-format data frame of Betas and P-values for each contrast and pattern
beta_list <- lapply(seq_along(lm_results), function(i) {
  model <- lm_results[[i]]
  pattern <- colnames(kelly)[i]
  
  contrast_data <- lapply(contrasts, function(con) {
    if (con %in% rownames(coef(summary(model)))) {
      est <- coef(summary(model))[con, "Estimate"]
      p <- coef(summary(model))[con, "Pr(>|t|)"]
    } else {
      est <- NA
      p <- NA
    }
    data.frame(Pattern = pattern, Contrast = con, Beta = est, P = p)
  })
  
  do.call(rbind, contrast_data)
})

# Combine all patterns into one data frame
beta_df <- do.call(rbind, beta_list)

# Clean up contrast labels for plotting
beta_df$Contrast <- factor(beta_df$Contrast, 
                           levels = contrasts, 
                           labels = c("GroupSCZ Low vs Control", "GroupSCZ High vs Control"))
beta_df$Pattern <- factor(beta_df$Pattern, levels = unique(beta_df$Pattern))

# Plot
pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kelly_pattern_betas.pdf", width = 10, height = 6)
ggplot(beta_df, aes(x = Pattern, y = Beta, fill = P < 0.05)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = ifelse(!is.na(P), sprintf("%.2g", P), "")), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3) +
  facet_wrap(~ Contrast) +
  scale_fill_manual(values = c("FALSE" = "#999999", "TRUE" = "#f9421b")) +
  theme_minimal() +
  labs(title = "Kelly Pattern Betas by Group Contrast",
       y = "Beta",
       fill = "P < 0.05") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# get top markers per pattern
markers <- patternMarkers(results_kornberg, threshold = "cut")  # returns a named vector of pattern assignment per gene
top_per_pattern <- lapply(seq_along(markers$PatternMarkers), function(p) {
  top_genes <- markers$PatternMarkers[[p]]
  # Sort these genes by their weight in pattern p
  top_genes <- top_genes[order(A_kornberg[top_genes, p], decreasing = TRUE)]
  head(top_genes, 10)
})
names(top_per_pattern) <- paste0("Pattern_", seq_along(top_per_pattern))

top_long <- do.call(rbind, lapply(names(top_per_pattern), function(pat) {
  data.frame(
    Pattern = pat,
    Rank = seq_along(top_per_pattern[[pat]]),
    Gene = top_per_pattern[[pat]]
  )
}))
write.csv(top_long, file = "top_marker_genes_long_30Kornberg.csv", row.names = FALSE)


# gsea
library(fgsea)
library(msigdbr)

genesets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gmt_list <- split(genesets$gene_symbol, genesets$gs_name)

# run fgsea
fgsea_results_list <- lapply(1:ncol(A_kornberg), function(pat_id) {
  gene_ranks <- A_kelly[, pat_id]
  names(gene_ranks) <- rownames(A_kornberg)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)

  fgseaRes <- fgsea(
    pathways = gmt_list,
    stats = gene_ranks,
    nperm = 1000
  )

  fgseaRes$pattern <- paste0("Pattern_", pat_id)
  return(fgseaRes)
})

# combine into one data.frame
fgsea_combined <- bind_rows(fgsea_results_list)
fgsea_export <- fgsea_combined %>% select(-leadingEdge)
write.csv(fgsea_export, "fgsea_results_all_30patterns_Kornberg.csv", row.names = FALSE)


## visualization ----

top_gsea <- fgsea_combined %>%
  group_by(pattern) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  ungroup()

plot_gsea_pattern <- function(df, pattern_id, wrap_width=40) {
  df_pattern <- df %>% 
	  	filter(pattern == pattern_id) %>%
    		mutate(pathway = gsub("_", " ", substring(tolower(pathway), 6, nchar(pathway))),
			   pathway_wrapped = str_wrap(pathway, width = wrap_width))

  ggplot(df_pattern, aes(x = reorder(pathway_wrapped, NES), y = NES, fill = NES)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = paste("GSEA for ", pattern_id),
         x = "Pathway",
         y = "Normalized Enrichment Score (NES)") +
    theme_minimal()
}

unique_patterns <- unique(top_gsea$pattern)

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/gsea_by_pattern_Kornberg.pdf", width = 8, height = 6)

for (p in unique_patterns) {
  print(plot_gsea_pattern(top_gsea, p))
}

dev.off()

