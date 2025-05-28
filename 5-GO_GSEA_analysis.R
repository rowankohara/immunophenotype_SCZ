# Rowan O'Hara-Payne
# 5/28/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Gene Ontology & Gene Set Enrichment Analysis of the K-means clustering


# Read in libraries and data ----
library(dplyr)
library(tibble)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(ggpubr)


# load modules and DGE results
modules_df <- read.csv("gene_modules_all_celltypes.csv")

dge <- read.table('dkelly_singlets.pseudobulk_edgeR.dge.txt', sep='\t')
names(dge) <- dge[1,]
dge <- dge[-1,]


# GO & GSEA Analysis ----

# load hallmark gene sets for fgsea
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# store results
go_results_all <- list()
gsea_results_all <- list()

celltypes <- unique(modules_df$celltype)

# GO and GSEA analysis per module
for (ct in celltypes) {
  dge_ct <- dge %>% filter(Celltype == ct)
  mod_df <- modules_df %>% filter(celltype == ct)

  for (mod in unique(mod_df$module)) {
    message("Analyzing ", ct, " module ", mod, "...")

    genes <- mod_df %>% filter(module == mod) %>% pull(gene)
    eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    # GO enrichment
    go <- enrichGO(gene = eg$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
    go_results_all[[paste(ct, mod, sep = "_")]] <- go

    if (!is.null(go) && nrow(go) > 0) {
      go_file <- paste0("plots/", gsub(" ", "_", ct), "_module_", mod, "_GO.pdf")
      pdf(go_file)
      print(barplot(go, showCategory = 20, title = paste("GO Enrichment -", ct, "Module", mod)))
      dev.off()
    }

    # GSEA
    gene_ranks <- dge_ct %>%
      group_by(Symbol) %>%
      summarise(score = mean(logFC, na.rm = TRUE)) %>%
      deframe()
    gene_ranks <- sort(gene_ranks, decreasing = TRUE)

    gsea <- fgsea(pathways = hallmark_sets,
                  stats = gene_ranks,
                  nperm = 1000)

    gsea_results_all[[paste(ct, mod, sep = "_")]] <- gsea

    if (!is.null(gsea) && nrow(gsea) > 0) {
      top_gsea <- gsea %>% arrange(padj) %>% head(20)
      gsea_file <- paste0("plots/", gsub(" ", "_", ct), "_module_", mod, "_GSEA.pdf")
      pdf(gsea_file)
      print(ggplot(top_gsea, aes(reorder(pathway, NES), NES, fill = NES)) +
              geom_col(show.legend = FALSE) +
              coord_flip() +
              labs(title = paste("GSEA -", ct, "Module", mod), x = "Pathway", y = "NES"))
      dev.off()
    }
  }
}

# save results
saveRDS(go_results_all, "go_results_all.rds")
saveRDS(gsea_results_all, "gsea_results_all.rds")



