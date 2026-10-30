# Rowan K. O'Hara-Payne
# 5/25/2025
# conda environment: rnaseq
# Differential abundance analysis using binomial generalized linear model
####################################################################################################################


# Read in libraries and data ----
library(Seurat)
library(tidyverse)
library(edgeR)
library(pheatmap)
library(emmeans)
library(ggpubr)

# read in Seurat object
obj <- readRDS('kelly_postqc_clusters_meta.rds')

# read in pseudobulk data
elist <- readRDS('dkelly_singlets.pseudobulk.rds')

# read in meta data
sample <- colnames(elist[[11]])
MPRCID <- gsub( '\\_(.*)' , '' , sample )
MPRCID <- gsub( '\\-(.*)' , '' , MPRCID )

meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
rownames(meta) <- meta$MPRCID
meta <- meta[ MPRCID , ]


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


# Differential Abundance Analysis ----

# create an empty data frame to store cell-level metadata
cell_meta <- data.frame(row.names = colnames(obj))

# loop through each metadata column (excluding MPRCID)
meta_filtered <- meta %>% select(MPRCID, Sex, Group, Age)
for (col in colnames(meta_filtered[-1])) {
  vec <- rep(NA, ncol(obj))
  names(vec) <- colnames(obj)
  
  for (i in 1:nrow(meta_filtered)) {
    sample_id <- meta_filtered$MPRCID[i]
    
    # find cells whose sample value starts with the sample_id
    matched_idx <- grep(paste0("^", sample_id), obj@meta.data$sample)
    matched_cells <- colnames(obj)[matched_idx]
    
    # assign sample-level metadata to those cells
    vec[matched_cells] <- meta_filtered[[col]][i]
  }
  
  cell_meta[[col]] <- vec
}

# add the metadata to the Seurat object
obj <- Seurat::AddMetaData(obj, metadata = cell_meta)

# summarize to one row per sample
sample_metadata <- obj@meta.data %>%
  select(sample, Group, Sex, Age) %>%
  distinct()

# cell counts per sample per cell type
counts_table <- table(obj$sample, obj$seurat_clusters)
counts_matrix <- as.matrix(counts_table)  # samples x cell types

# match sample order
sample_metadata <- sample_metadata %>%
  filter(sample %in% rownames(counts_matrix)) %>%
  arrange(match(sample, rownames(counts_matrix)))

# confirm alignment
stopifnot(all(sample_metadata$sample == rownames(counts_matrix)))

# convert categorical variables to factors
sample_metadata$Group <- factor(sample_metadata$Group)
sample_metadata$Group <- forcats::fct_recode(sample_metadata$Group, "SCZLow" = "SCZ Low", "SCZHigh" = "SCZ High")
sample_metadata$Sex <- factor(sample_metadata$Sex)

# convert to long format and compute proportions
prop_df <- counts_table %>%
  as.data.frame() %>%
  rename(sample = Var1, cell_type = Var2, count = Freq) %>%
  group_by(sample) %>%
  mutate(total = sum(count),
         prop = count / total) %>%
  ungroup()

# add sample meta data
prop_df <- left_join(prop_df, sample_metadata, by = "sample")

### Figs 2D-K ----
run_binom_model <- function(df_ct) {

  ct_name <- as.character(unique(df_ct$cell_type)[1])

  group_levels <- c("SCZHigh", "SCZLow", "HC")
  df_ct$Group <- factor(df_ct$Group, levels = group_levels)

  # fit model
  fit <- glm(
    cbind(count, total - count) ~ Group + Sex + Age,
    family = binomial,
    data = df_ct
  )

  # overall test
  anova_res <- car::Anova(fit, type = 3)
  p_group <- anova_res["Group", "Pr(>Chisq)"]

  # Emmeans (response scale)
  emm <- emmeans(fit, ~ Group, type = "response")
  emm_df <- as.data.frame(emm) %>%
    mutate(cell_type = ct_name)


  # split contrast
  contr <- emmeans::contrast(emm, method = "pairwise", adjust = "none")

  pairs_df <- as.data.frame(contr) %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " / ") %>%
    dplyr::mutate(
      group1 = trimws(group1),
      group2 = trimws(group2),
      p_label = sprintf("p=%.3g", p.value),
      cell_type = ct_name
    )
  

  # y position for significance bars
  y_max <- max(c(df_ct$prop, emm_df$prob), na.rm = TRUE)
  y_range <- diff(range(c(df_ct$prop, emm_df$prob), na.rm = TRUE))

  pairs_df <- pairs_df %>%
    dplyr::arrange(group1, group2) %>%
    dplyr::mutate(
      y.position = y_max + 0.08 * y_range * seq_len(n())
    )

  # plot
  p <- ggplot(df_ct, aes(x = Group, y = prop)) +
	  geom_col(data = emm_df,
                aes(x = Group, y = prob, fill = Group), color = "black", lwd = 1.25, width = 0.75) +
	  geom_jitter(width = 0.15, size = 4, alpha = 0.8, color = "black") +
	  geom_errorbar(
    		data = emm_df,
    		aes(x = Group, ymin = asymp.LCL, ymax = asymp.UCL),
    		width = 0.35,
		lwd = 1.25,
    		inherit.aes = FALSE) +
	  scale_x_discrete(
			   labels = c("SCZHigh" = "AGA+ SRD",
				      "SCZLow" = "AGA- SRD",
				      "HC" = "AGA- Control")
			   ) +
	  scale_fill_manual(values = c("#BB5566", "#DDAA33", "#004488")) +
	  labs(title = parse(text = celltype_labels[[ct_name]]),
	       y = "Proportion",
	       fill = "group") +
	  theme_classic() +
	  scale_y_continuous(limits = c(0, NA),
			     expand = expansion(mult = c(0, 0.05))
			     ) +
	  theme(plot.title = element_text(hjust = 0.5,
					  size = 30),
		axis.text = element_text(size = 20),
		axis.title = element_text(size = 17),
		axis.line = element_line(linewidth = 1.25),
		axis.ticks = element_line(linewidth = 1.25))

  # add significance bars
  p <- p + ggpubr::stat_pvalue_manual(
				      pairs_df,
    				      label = "p_label",
				      xmin = "group1",
				      xmax = "group2",
				      y.position = "y.position",
				      tip.length = 0.015,
				      inherit.aes = FALSE,
				      size = 5,
				      bracket.size = 0.7) +
			coord_cartesian(clip = "off")

  return(list(
    summary = data.frame(cell_type = ct_name, p_value = p_group),
    emm = emm_df,
    pairs = pairs_df,
    plot = p
  ))
}

results_list <- prop_df %>%
  group_split(cell_type) %>%
  map(run_binom_model)

summary_df <- map_dfr(results_list, "summary") %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

emm_all <- map_dfr(results_list, "emm")

pairs_all <- map_dfr(results_list, "pairs") %>%
  mutate(FDR = p.adjust(p.value, method = "fdr"))

# print box plot for every cell type
pdf("plots/celltype_proportion_plots.pdf", width = 12, height = 8)
walk(results_list, ~ print(.x$plot))
dev.off()

