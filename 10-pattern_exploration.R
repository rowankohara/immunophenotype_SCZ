# Rowan K. O'Hara-Payne
# 9/24/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Explore the Kelly patterns that were decomposed


# Read in libraries and data ----
library(Seurat)
library(CoGAPS)
library(tidyverse)

kelly_cogaps <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kellyCoGAPS_30Pat.rds")
kelly_obj <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kelly_postqc_celltypes_meta_cogaps.rds")
kelly_obj <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/dkelly_postqc_celltypes_meta.rds")




# get patterns in long format
P <- kelly_cogaps@sampleFactors
P_long <- as.data.frame(P) %>%
		rownames_to_column(var = "cell_id") %>%
		pivot_longer(
    		cols = starts_with("Pattern"),
    		names_to = "Pattern",
    		values_to = "Activity"
  	)


# get cell meta data (including sample meta data like age, sex, and group)
kelly_meta <- as.data.frame(kelly_obj@metadata)
kelly_meta$cell_id <- rownames(kelly_meta)


# join patterns and meta data
kelly_df <- left_join(kelly_meta, P_long, by = "cell_id")


# plot a pattern for the cell types that had significant differences between SCZ High and SCZ Low
pdf("plots/kelly_cogaps_violin_celltype.pdf", width=10, height=6)
df <- kelly_df %>%
	filter(Pattern == "Pattern_30",
	       celltype %in% c("ILC"))
ggplot(df, aes(x = condition, y = Activity, fill = condition)) +
	geom_violin(scale = "width") +
	#facet_wrap(~celltype, scales = "free_y") +
	theme_bw() +
	labs(title = "Pattern 30: ILC",
	     x = "Group",
	     y = "Activity Score",
	     color = "Group")
dev.off()




