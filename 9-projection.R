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
library(lme4)
library(lmerTest)
library(emmeans)


# read in CoGAPS results from single cell datasets

# Kelly data
kelly <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kellyCoGAPS_30Pat.rds")

kelly_counts <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/dkelly_singlets.postqc.azimuth_celltypes.rds")
kelly_counts$celltype <- kelly_counts$predicted.celltype.l2
kelly_counts$pseudobulk_id <- paste(gsub("_", "-", kelly_counts$sample), kelly_counts$celltype, sep="_")
kelly_pseudo <- AverageExpression(kelly_counts, group.by = c("sample", "celltype"), slot = "data")
kelly_pseudo_mat <- as.matrix(kelly_pseudo$RNA)

# Kornberg data
kornberg <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kornbergCoGAPS_30Pat.rds")

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
kelly_kornberg_proj <- projectR(data=kornberg_pseudo_mat, loadings=kelly) # pseudobulked
write.csv(kelly_kornberg_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kelly_kornberg_proj.csv")

kornberg_rna <- LayerData(object = kornberg_counts, assay = "RNA", layer = "data")
kelly_kornberg_proj <- projectR(data=as.matrix(kornberg_rna), loadings = kelly)
write.csv(kelly_kornberg_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kelly_kornberg_proj_sc.csv", rownames = T)

# project Kelly patterns onto pseudobulked Kelly data for sanity check
kelly_kelly_proj <- projectR(data=kelly_pseudo_mat, loadings=kelly@featureLoadings) # pseudobulked
write.csv(kelly_kelly_proj, file="/local/projects/idea/rohara/immunophenotype_SCZ/kelly_kelly_pseudo_proj.csv")

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

# for pseudobulked
kornberg_meta <- kornberg_meta[!duplicated(kornberg_meta$pseudobulk_id), ] %>%
			select(pseudobulk_id, sample, general_condition, condition, sample_age, sample_sex, celltype)
rownames(kornberg_meta) <- kornberg_meta$pseudobulk_id

# not pseudobulked
kornberg_meta$cell_id <- rownames(kornberg_meta)
kornberg_meta$group <- kornberg_meta$general_condition
kornberg_meta <- kornberg_meta %>% select(cell_id, sample, group, condition, sample_age, sample_sex, celltype)

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

# for pseudobulked
kelly_meta <- kelly_meta[!duplicated(kelly_meta$pseudobulk_id), ] %>%
                        select(pseudobulk_id, sample, group, sample_age, sample_sex, celltype)
rownames(kelly_meta) <- kelly_meta$pseudobulk_id

# not pseudobulked
kelly_meta$cell_id <- kelly_meta$bc_wells
kelly_meta <- kelly_meta %>% select(cell_id, sample, group, sample_age, sample_sex, celltype) 

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
  formula <- as.formula(paste(pat, "~ group + celltype + sample_sex + sample_age"))
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
A <- kelly@featureLoadings
P <- kelly@sampleFactors
P_long <- as.data.frame(P) %>%
		pivot_longer(
    		cols = starts_with("Pattern"),
    		names_to = "Pattern",
    		values_to = "Activity"
  	)
kelly_df <- as.data.frame(P) %>%
  mutate(cell_id = rownames(.)) %>%
  pivot_longer(cols = starts_with("Pattern"), names_to = "Pattern", values_to = "Score") %>%
  left_join(kelly_meta, by = "cell_id") %>%
  mutate(Dataset = "Kelly")
kelly_df$celltype <- factor(kelly_df$celltype, levels = unique(kelly_df$celltype))
kelly_df$sample <- factor(kelly_df$sample, levels = unique(kelly_df$sample))
kelly_df$group <- factor(kelly_df$group, levels = unique(kelly_df$group))

kornberg_df <- kelly_kornberg_data %>%
  pivot_longer(cols = starts_with("Pattern"), names_to = "Pattern", values_to = "Score") %>%
  mutate(Dataset = "Kornberg")
kornberg_df$celltype <- factor(kornberg_df$celltype, levels = unique(kornberg_df$celltype))
kornberg_df$sample <- factor(kornberg_df$sample, levels = unique(kornberg_df$sample))
kornberg_df$group <- factor(kornberg_df$group, levels = unique(kornberg_df$group))
kornberg_df$sample_sex <- factor(kornberg_df$sample_sex, levels = unique(kornberg_df$sample_sex))

combined <- bind_rows(kelly_df, kornberg_df)
combined$Dataset <- factor(combined$Dataset, levels = unique(combined$Dataset))
sorted_levels <- sort(levels(combined$celltype))
combined$celltype <- factor(combined$celltype, levels = sorted_levels)

pdf("pattern_violin_CoGAPS_proj_comparison.pdf", width = 12, height = 6)
for (pattern in unique(combined$Pattern)) {
  print(
    ggplot(subset(combined, Pattern == pattern), aes(x = celltype, y = Score, fill = Dataset)) +
      geom_violin(scale = "width", trim = FALSE, position = position_dodge(0.8)) +
      theme_minimal() +
      labs(title = paste("Comparison for:", pattern), y = "Acivity Score", x = "Cell Type") +
      scale_fill_manual(values = c("Kelly" = "tomato", "Kornberg" = "skyblue")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  )
}
dev.off()


# redo modeling

kelly_model <- kelly_df %>%
	pivot_wider(names_from = "Pattern", values_from = "Score")
kornberg_model <- kornberg_df %>%
	pivot_wider(names_from = "Pattern", values_from = "Score")
patterns <- colnames(kelly_model)[8:37]

lm_results <- lapply(patterns, function(pat) {
  formula <- as.formula(paste(pat, "~ group * celltype + sample_sex + sample_age + ( 1 | sample)"))
  lmer(formula, data = kornberg_model)
})

names(lm_results) <- patterns


extract_lmer_results <- function(models) {
  # --- Fixed effects ---
  fixed_list <- lapply(seq_along(models), function(i) {
    m <- models[[i]]
    s <- summary(m)
    fe <- as.data.frame(s$coefficients)
    fe$term <- rownames(fe)
    fe$pattern <- names(models)[i]
    rownames(fe) <- NULL
    fe
  })
  fixed_df <- do.call(rbind, fixed_list)
  
  # --- Random effects (variance components) ---
  random_list <- lapply(seq_along(models), function(i) {
    m <- models[[i]]
    vc <- as.data.frame(VarCorr(m))
    vc$pattern <- names(models)[i]
    vc
  })
  random_df <- do.call(rbind, random_list)
  
  return(list(fixed = fixed_df, random = random_df))
}

extract_all_results <- function(models) {
  results_list <- lapply(seq_along(models), function(i) {
    m <- models[[i]]
    model_name <- names(models)[i]
    
    ## Fixed effects
    fe <- coef(summary(m))
    fe_df <- data.frame(
      term = rownames(fe),
      estimate = fe[, "Estimate"],
      std_error = fe[, "Std. Error"],
      df = fe[, "df"],
      t_value = fe[, "t value"],
      p_value = fe[, "Pr(>|t|)"],
      effect_type = "fixed",
      pattern = model_name,
      row.names = NULL
    )
    
    ## Random effects
    re <- ranef(m, condVar = TRUE)
    re_list <- lapply(names(re), function(fac) {
      re_df <- as.data.frame(re[[fac]])
      re_df$grouping_factor = fac
      re_df
    })
    re_df <- do.call(rbind, re_list)
    re_df$effect_type <- "random"
    re_df$pattern <- model_name
    
    ## emmeans contrasts
    emm <- emmeans(m, ~ group | celltype)
    contr <- contrast(emm, method = "pairwise", by = "celltype")
    contr_df <- as.data.frame(contr)
    contr_df$effect_type <- "emmeans_contrast"
    contr_df$pattern <- model_name
    
    ## Bind all together (keep same columns)
    combined <- list(
      fe_df,
      re_df,
      contr_df
    )
    do.call(plyr::rbind.fill, combined)  # allows different columns
  })
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  return(results_df)
}


extract_random_var <- function(models) {
  res_list <- lapply(seq_along(models), function(i) {
    m <- models[[i]]
    name <- names(models)[i]
    
    vc <- as.data.frame(VarCorr(m))
    vc$pattern <- name
    
    # Add ICC
    donor_var <- vc$vcov[vc$grp == "sample"]  # replace with your grouping factor name
    resid_var <- vc$vcov[vc$grp == "Residual"]
    vc$icc <- if (!is.na(donor_var) && !is.na(resid_var)) {
      donor_var / (donor_var + resid_var)
    } else {
      NA
    }
    
    vc
  })
  
  do.call(rbind, res_list)
}


results <- extract_lmer_results(lm_results)
all_results <- extract_all_results(lm_results)
random_variances <- extract_random_var(lm_results)

# Save separately
write.csv(results$fixed, "kornberg_lmer_fixed_effects.csv", row.names = FALSE)
write.csv(results$random, "kornberg_lmer_random_effects.csv", row.names = FALSE)
write.csv(all_results, "kornberg_full_model_results.csv", row.names = FALSE)
write.csv(random_variances, "kornbberg_random_effect_variances.csv", row.names = FALSE)


# make violin plots by sample from the Kelly sc projection onto pseudobulked Kelly

# convert to tibble, keeping pattern/sample labels
proj_df <- kelly_kelly_proj %>%
  as.data.frame() %>%
  rownames_to_column("Pattern") %>%
  pivot_longer(
    -Pattern,
    names_to = "pseudobulk_id",
    values_to = "score"
  ) %>%
  mutate(dataset = "pseudobulk")
 kelly_meta$sample <- gsub("_", "-", kelly_meta$sample)
 kelly_meta$pseudobulk_id <- paste(kelly_meta$sample, kelly_meta$celltype, sep="_")

# join with your metadata (assume sample_metadata has Sample + Group columns)
proj_df <- proj_df %>%
  left_join(kelly_meta, by = "pseudobulk_id", multiple = "first")


P$cell_id <- rownames(P)  # make sure P is a data frame
sc_df <- merge(P, kelly_meta, by = "cell_id")  # ran into memory limit error when ran this
sc_df <- sc_df %>%
	pivot_longer(cols = starts_with("Pattern"),
		     names_to = "Pattern",
		     values_to = "score")
sc_df$dataset <- "single cell"

plot_df <- bind_rows(
  sc_df %>% select(sample, Pattern, score, dataset),
  proj_df %>% select(sample, Pattern, score, dataset)
)

pdf("plots/kelly_sc_proj_pattern_violin.pdf", width=10, height=6)
for (pat in colnames(P)[1:30]) {
	plot_df2 <- plot_df %>% filter(Pattern == pat)
	p <- plot_df2 %>%
  		ggplot(aes(x = sample, y = score, fill = dataset)) +
  		geom_violin(position = position_dodge(width = 0.8), trim = FALSE) +
  		theme_classic() +
  		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  		labs(y = "Pattern activity",
		       x = "Samples",
		       title = paste(pat))
	print(p)

}
dev.off()


