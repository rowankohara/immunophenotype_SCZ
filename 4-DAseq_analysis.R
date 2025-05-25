# Rowan K. O'Hara-Payne
# 5/25/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Do differential abundance analysis


# Read in libraries and data ----
library(Seurat)
library(dplyr)
library(speckle)   # contains propeller
library(ggplot2)

# read in Seurat object
obj <- readRDS('dkelly_singlets.postqc.azimuth_celltypes.rds')

# read in pseudobulk data
elist <- readRDS('dkelly_singlets.pseudobulk.rds')

# read in meta data
sample <- colnames(elist[[11]])
MPRCID <- gsub( '\\_(.*)' , '' , sample )
MPRCID <- gsub( '\\-(.*)' , '' , MPRCID )

meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
rownames(meta) <- meta$MPRCID
meta <- meta[ MPRCID , ]


# Differential Abundance Analysis ----

# create a vector to hold the sample-level metadata aligned to each cell
Group <- rep(NA, ncol(obj))
names(Group) <- colnames(obj)  # recommended for clarity

# loop through each sample and assign metadata
for (i in seq_len(nrow(meta))) {
  sample_id <- meta$MPRCID[i]
  matched_cells <- colnames(obj)[grep(paste0("^", sample_id), obj@meta.data$sample)]

  Group[matched_cells] <- meta$Group[i]
}

# assign to Seurat object
obj$Group <- Group

meta_df <- obj@meta.data %>%
        select(predicted.celltype.l2, sample, Group)
counts_table <- table(meta_df$sample, meta_df$predicted.celltype.l2)

# ensure correct order of samples
sample_info <- meta_df %>%
  distinct(sample, Group) %>%
  arrange(match(sample, rownames(counts_table)))

group_vector <- sample_info$Group

# run with logit transformation
prop_res <- speckle::propeller.anova(
  counts = counts_table,
  group = group_vector,
  transform = "logit"  # or "asin"
)

# view results sorted by adjusted p-value
prop_res %>%
  arrange(p.adj) %>%
  filter(p.adj < 0.05)

# compute proportions manually for plotting
plot_df <- as.data.frame(counts_table) %>%
  rename(sample_id = Var1, cell_type = Var2, count = Freq) %>%
  group_by(sample) %>%
  mutate(prop = count / sum(count)) %>%
  left_join(sample_info, by = "sample")

# plot using ggplot2
pdf("plots/celltype_prop_DAseq.pdf")
ggplot(plot_df, aes(x = cell_type, y = prop, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Cell Type Proportion", x = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


