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

# Kornberg data
kornberg <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kornbergCoGAPS_30Pat.rds")
kornberg <- kornberg@featureLoadings


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
kelly_kornberg_proj <- projectR(data=kornberg, loadings=kelly)
write.csv(kelly_kornberg_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kelly_kornberg_proj.csv")

# project Kornberg patterns onto Kelly data (schizophrenia and gluten intolerance immunophenotype)
kornberg_kelly_proj <- projectR(data=kelly, loadings=kornberg)
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
colnames(pan_meta) <- c("group", "sex")

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
kornberg_seurat <- readRDS("GSE252762_kornberg_han_celiac_dataset_fig1-2.RDS")
kornberg_meta <- as.data.frame(kornberg_seurat@meta.data)
kornberg_meta <- kornberg_meta[!duplicated(kornberg_meta$sample), ] %>%
			select(sample, general_condition, condition, sample_age, sample_sex)
rownames(kornberg_meta) <- kornberg_meta$sample



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


lm_results <- lapply(colnames(t(kelly_glahn_proj)), function(pat) {
  formula <- as.formula(paste(pat, "~ group + sex + age"))
  lm(formula, data = kelly_glahn_data)
})

names(lm_results) <- colnames(t(kelly_glahn_proj))

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

results_df <- do.call(rbind, Map(extract_coefs, lm_results, names(lm_results)))
write.csv(results_df, "/local/projects/idea/rohara/immunophenotype_SCZ/kelly_pan_proj_model_results.csv", row.names = FALSE)


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



