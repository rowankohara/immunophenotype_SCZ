# Rowan K. O'Hara-Payne
# 7/21/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Using projectR to project the patterns from CoGAPS onto other datasets


# Load libraries and read in data ----
library(Seurat)
library(CoGAPS)
library(projectR)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(edgeR)
library(tidyverse)


# read in CoGAPS results from single cell datasets

# Kelly data
kelly <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kellyCoGAPS_30Pat.rds")
kelly <- kelly@featureLoadings

kelly_counts <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/dkelly_singlets.postqc.azimuth_celltypes.rds")
kelly_counts$celltype <- kelly_counts$predicted.celltype.l2
kelly_counts$pseudobulk_id <- paste(gsub("_", "-", kelly_counts$sample), kelly_countsc$celltype, sep="_")
kelly_pseudo <- AverageExpression(kelly_counts, group.by = c("sample", "celltype"), slot = "data")
kelly_pseudo_mat <- as.matrix(kelly_pseudo$RNA)

# Kornberg data
kornberg <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kornbergCoGAPS_30Pat.rds")
kornberg <- kornberg@featureLoadings

# read in Kornberg Seurat object and add Azimuth celltype predictions
kornberg_counts <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/GSE252762_kornberg_han_celiac_dataset_fig1-2.RDS")
kornberg_counts[['percent.mt']] <- PercentageFeatureSet(kornberg_counts, pattern = '^MT-')
kornberg_counts <- subset(kornberg_counts, nCount_RNA < 10000 & nCount_RNA > 500 & percent.mt < 15)


predictions <- read.delim("/local/projects/idea/rohara/immunophenotype_SCZ/azimuth_pred.tsv", row.names = 1)
kornberg_counts <- AddMetaData(
        object = kornberg_counts,
        metadata = predictions)

kornberg_counts$celltype <- kornberg_counts$predicted.celltype.l2
kornberg_counts$pseudobulk_id <- paste(gsub("_", "-", kornberg_counts$sample), kornberg_counts$celltype, sep="_")
kornberg_pseudo <- AverageExpression(kornberg_counts, group.by = c("sample", "celltype"), slot = "data")
kornberg_pseudo_mat <- as.matrix(kornberg_pseudo$RNA)


# read in bulk datasets

# Pan data
pan <- read.table(file="/local/projects/idea/rohara/immunophenotype_SCZ/GSE243841_cleaned_counts.tsv", sep="\t", header=T)
rownames(pan) <- pan$Geneid  # set rownames to gene ids
pan <- pan[ , -1]  # remove gene id column (matrix can be only one data type)
pan <- as.matrix(pan)
gene_map <- mapIds(  # map gene names to be the same as the other datasets
  org.Hs.eg.db,
  keys = rownames(pan),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
idx <- which(!is.na(gene_map))
found_genes <- gene_map[idx]
pan <- pan[idx, ]
rownames(pan) <- found_genes
# normalize counts
dge <- DGEList(counts = pan)
dge <- calcNormFactors(dge)
pan <- cpm(dge, log = TRUE)

# Glahn data
glahn <- read.csv("/local/projects/idea/rohara/immunophenotype_SCZ/salmon_merged_gene_counts.csv")
rownames(glahn) <- glahn$gene_id
glahn <- glahn[ , -1]
glahn <- as.matrix(glahn)
# normalize counts
dge <- DGEList(counts = glahn)
dge <- calcNormFactors(dge)
glahn <- cpm(dge, log = TRUE)


# Projection ----

## project single cell datasets onto bulk ----

# project Kelly patterns onto Pan data (schizophrenia immunophenotype)
kelly_pan_proj <- projectR(data=pan, loadings=kelly)
write.csv(kelly_pan_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kelly_pan_proj.csv")

# project Kelly patterns onto Glahn data (psychosis)
kelly_glahn_proj <- projectR(data=glahn, loadings=kelly)
write.csv(kelly_glahn_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kelly_glahn_proj.csv")

# project Kornberg patterns onto Pan data (schizophrenia immunophenotype)
kornberg_pan_proj <- projectR(data=pan, loadings=kornberg)
write.csv(kornberg_pan_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kornberg_pan_proj.csv")

# project Kornberg patterns onto Glahn data (psychosis)
kornberg_glahn_proj <- projectR(data=glahn, loadings=kornberg)
write.csv(kornberg_glahn_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kornberg_glahn_proj.csv")


## project single cell onto single cell ----

# project Kelly patterns onto Kornberg data (gluten intolerance)
kelly_kornberg_proj <- projectR(data=kornberg_pseudo_mat, loadings=kelly)
write.csv(kelly_kornberg_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kelly_kornberg_proj.csv")

# project Kornberg patterns onto Kelly data (schizophrenia and gluten intolerance immunophenotype)
kornberg_kelly_proj <- projectR(data=kelly_pseudo_mat, loadings=kornberg)
write.csv(kornberg_kelly_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kornberg_kelly_proj.csv")


# Visualization ----

# Pan meta data
# ow samples = SCZ, YD samples = HC
samples <- colnames(pan)
scz <- grepl("^ow", samples)
pan_meta <- data.frame(sample = samples, temp = scz)
pan_meta <- pan_meta %>%
		mutate(group = if_else(temp == TRUE, "SCZ", "HC")) %>%
	       		select(-temp)
rownames(pan_meta) <- samples
pan_group <- pan_meta %>% select(-sample)

# determine sex
female_marker <- "XIST"
male_markers <- c("RPS4Y1", "KDM5D", "UTY")
sex_genes <- c(female_marker, male_markers)
expr_sex_genes <- pan[sex_genes, , drop = FALSE]
sex_df <- as.data.frame(t(expr_sex_genes))
sex_df$predicted_sex <- ifelse(
  sex_df$XIST > 1 & rowMeans(sex_df[, male_markers]) < 1, "Female", "Male"
)
pan_meta <- cbind(pan_meta, sex_df$predicted_sex)
colnames(pan_meta) <- c("sample", "group", "sex")

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kelly_pan_proj_heatmap.pdf")
pheatmap(kelly_pan_proj, 
	 scale = "row",
	 main = "Kelly Projection onto Pan",
	 show_colnames = FALSE,
	 annotation_col = pan_group)
dev.off()


pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kornberg_pan_proj_heatmap.pdf")
pheatmap(kornberg_pan_proj, 
         scale = "row",
         main = "Kornberg Projection onto Pan",
         show_colnames = FALSE,
         annotation_col = pan_group)
dev.off()

# Glahn meta data
glahn_meta <- read.csv("/local/projects/idea/rohara/immunophenotype_SCZ/redcap_selection_for_seth.csv")
glahn_meta <- glahn_meta %>%
		mutate(group = if_else(subject_type == 888, "PSY", "HC"),
		       sample = paste0("X", as.character(subject_type), as.character(id)))
rownames(glahn_meta) <- glahn_meta$sample
glahn_group <- glahn_meta %>% select(group)

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kelly_glahn_proj_heatmap.pdf")
pheatmap(kelly_glahn_proj,
         scale = "row",
         main = "Kelly Projection onto Glahn",
         show_colnames = FALSE,
         annotation_col = glahn_group)
dev.off()


pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kornberg_glahn_proj_heatmap.pdf")
pheatmap(kornberg_glahn_proj,
         scale = "row",
         main = "Kornberg Projection onto Glahn",
         show_colnames = FALSE,
         annotation_col = glahn_group)
dev.off()


# Kornberg meta data
kornberg_meta <- as.data.frame(kornberg_counts@meta.data)
kornberg_meta <- kornberg_meta[!duplicated(kornberg_meta$pseudobulk_id), ] %>%
			select(pseudobulk_id, sample, general_condition, condition, sample_age, sample_sex, celltype)
rownames(kornberg_meta) <- kornberg_meta$pseudobulk_id
kornberg_group <- kornberg_meta$general_condition

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kelly_kornberg_proj_heatmap.pdf")
pheatmap(kelly_kornberg_proj,
         scale = "row",
         main = "Kelly Projection onto Kornberg",
	 show_colnames = FALSE)
dev.off()

# Kelly meta data
kelly_meta <- as.data.frame(kelly_counts@meta.data)
MPRCID <- gsub( '\\_(.*)' , '' , kelly_meta$sample )
meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
sample_ids <- unique(kelly_meta$sample)
# match each MPRCID to a sample ID that is a prefix of it
meta$sample <- sapply(meta$MPRCID, function(id) {
  match <- sample_ids[startsWith(sample_ids, id)]
  if (length(match) > 0) match[1] else NA
})
kelly_meta <- left_join(kelly_meta, meta, by = "sample")
kelly_meta$group <- kelly_meta$Group
kelly_meta$sample_sex <- kelly_meta$Sex
kelly_meta$sample_age <- kelly_meta$Age
kelly_meta <- kelly_meta[!duplicated(kelly_meta$pseudobulk_id), ] %>%
                        select(pseudobulk_id, sample, group, sample_age, sample_sex, celltype)
rownames(kelly_meta) <- kelly_meta$pseudobulk_id
kelly_meta <- kelly_meta %>%
		mutate(sample_sex = factor(sample_sex),
		       sample_sex = recode(sample_sex,
					   "Male" = "M",
					   "Female" = "F"))
kelly_group <- kelly_meta$group

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kornberg_kelly_proj_heatmap.pdf")
pheatmap(kornberg_kelly_proj,
         scale = "row",
         main = "Kornberg Projection onto Kelly",
	 show_colnames = FALSE)
dev.off()


# Statistical Modeling ----


kelly_pan_data <- t(kelly_pan_proj)
kelly_pan_data <- cbind(kelly_pan_data, pan_meta)

kornberg_pan_data <- t(kornberg_pan_proj)
kornberg_pan_data <- cbind(kornberg_pan_data, pan_meta)

kelly_glahn_data <- t(kelly_glahn_proj)
kelly_glahn_data <- merge(kelly_glahn_data, glahn_meta, by = "row.names")
rownames(kelly_glahn_data) <- kelly_glahn_data$Row.names
kelly_glahn_data <- kelly_glahn_data[ , -1]

kornberg_glahn_data <- t(kornberg_glahn_proj)
kornberg_glahn_data <- merge(kornberg_glahn_data, glahn_meta, by = "row.names")
rownames(kornberg_glahn_data) <- kornberg_glahn_data$Row.names
kornberg_glahn_data <- kornberg_glahn_data[ , -1]

kelly_kornberg_data <- t(kelly_kornberg_proj)
kelly_kornberg_data <- merge(kelly_kornberg_data, kornberg_meta, by = "row.names")
rownames(kelly_kornberg_data) <- kelly_kornberg_data$Row.names
kelly_kornberg_data <- kelly_kornberg_data[ , -1]

kornberg_kelly_data <- t(kornberg_kelly_proj)
kornberg_kelly_data <- merge(kornberg_kelly_data, kelly_meta, by = "row.names")
rownames(kornberg_kelly_data) <- kornberg_kelly_data$Row.names
kornberg_kelly_data <- kornberg_kelly_data[ , -1]

lm_results <- lapply(colnames(t(kornberg_kelly_proj)), function(pat) {
  formula <- as.formula(paste(pat, "~ group + sample_sex + sample_age + celltype"))
  lm(formula, data = kornberg_kelly_data)
})

names(lm_results) <- colnames(t(kornberg_kelly_proj))

extract_coefs <- function(model, pattern_name) {
  coefs <- summary(model)$coefficients
  
  data.frame(
    Pattern = pattern_name,
    
    groupSCZ_Estimate = coefs["groupSCZ", "Estimate"],
    groupSCZ_StdError = coefs["groupSCZ", "Std. Error"],
    groupSCZ_t        = coefs["groupSCZ", "t value"],
    groupSCZ_p        = coefs["groupSCZ", "Pr(>|t|)"],
    
    sexMale_Estimate  = coefs["sexMale", "Estimate"],
    sexMale_StdError  = coefs["sexMale", "Std. Error"],
    sexMale_t         = coefs["sexMale", "t value"],
    sexMale_p         = coefs["sexMale", "Pr(>|t|)"]

  )
}

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

results_df <- do.call(rbind, Map(extract_coefs, lm_results, names(lm_results)))
write.csv(results_df, "/local/projects/idea/rohara/immunophenotype_SCZ/kornberg_kelly_proj_model_results.csv", row.names = FALSE)


beta_df <- data.frame(
  Pattern = colnames(t(kelly_glahn_proj)),
  Beta = sapply(lm_results, function(x) coef(summary(x))["groupPSY", "Estimate"]),
  P = sapply(lm_results, function(x) coef(summary(x))["groupPSY", "Pr(>|t|)"])
)
beta_df$Pattern <- factor(beta_df$Pattern, levels = beta_df$Pattern)

pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/kelly_glahn_proj_betas.pdf")
ggplot(beta_df, aes(x = Pattern, y = Beta, fill = P < 0.05)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values = c("FALSE" = "#999999", "TRUE" = "#f9421b")) +
  geom_text(aes(label = sprintf("%.2g", P)), vjust = -0.5, size = 3) +
  labs(title = "Kelly Projection onto Glahn") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()



# check patterns
kelly <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kellyCoGAPS_30Pat.rds")
kelly <- kelly@featureLoadings
kelly_kornberg_proj <- read.csv("/local/projects/idea/rohara/immunophenotype_SCZ/kelly_kornberg_proj.csv")






