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
      ggplot(top_gsea, aes(reorder(pathway, NES), NES, fill = NES)) +
              geom_col(show.legend = FALSE) +
              coord_flip() +
              labs(title = paste("GSEA -", ct, "Module", mod), x = "Pathway", y = "NES")
      dev.off()
    }

    if (nrow(gsea) > 0) {
	# save table of top pathways
  	top_pathways <- gsea %>%
    		arrange(padj) %>%
    		head(10)

  	# write CSV
  	write.csv(top_pathways, file = paste0("gsea_", module_id, ".csv"), row.names = FALSE)

  	# plot top pathways
  	plot_file <- paste0("plots/gsea_dotplot_", module_id, ".pdf")
  	pdf(plot_file, width = 8, height = 6)
  	ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, color = padj)) +
    		geom_point(size = 3) +
    		coord_flip() +
    		scale_color_gradient(low = "red", high = "blue") +
    		labs(title = paste("GSEA -", module_id), x = "Pathway", y = "NES") +
    		theme_bw()
  	dev.off()
    }
  }
}

# save results
saveRDS(go_results_all, "go_results_all.rds")
saveRDS(gsea_results_all, "gsea_results_all.rds")


# Summary Results ----

# summarize all GO results
go_summary <- do.call(rbind, lapply(names(go_results_all), function(key) {
  go <- go_results_all[[key]]
  if (is.null(go) || nrow(go) == 0) return(NULL)

  # Split key into cell type and module
  parts <- strsplit(key, "_")[[1]]
  celltype <- paste(parts[1:(length(parts) - 1)], collapse = "_")
  module <- parts[length(parts)]

  go_df <- as.data.frame(go)
  go_df$Celltype <- celltype
  go_df$Module <- module
  go_df
}))

# reorder columns
go_summary <- go_summary %>%
  dplyr::select(Celltype, Module, ID, Description, GeneRatio, BgRatio,
         pvalue, p.adjust, qvalue, Count, everything())

# save to CSV
write.csv(go_summary, "GO_enrichment_summary.csv", row.names = FALSE)

# loop through cell types
for (ct in unique(go_summary$Celltype)) {
  go_ct <- go_summary %>%
    filter(Celltype == ct) %>%
    group_by(Module) %>%
    slice_min(order_by = p.adjust, n = 3, with_ties = FALSE) %>%
    ungroup()

  # Skip if no rows
  if (nrow(go_ct) == 0) {
    message("Skipping ", ct, " — no GO terms.")
    next
  }

  # shorten GO term descriptions for plot
  go_ct$DescriptionShort <- stringr::str_trunc(go_ct$Description, 60)
  go_ct$Label <- paste0(go_ct$Module, ": ", go_ct$DescriptionShort)

  # order y-axis
  go_ct <- go_ct %>%
    arrange(p.adjust) %>%
    mutate(Label = factor(Label, levels = rev(unique(Label))))

  # define file path
  dotplot_file <- paste0("plots/", gsub(" ", "_", ct), "_go_summary_dotplot.pdf")

  # open PDF and print plot
  pdf(dotplot_file)
  print(ggplot(go_ct, aes(x = -log10(p.adjust), y = Label, size = Count, color = factor(Module))) +
      geom_point(alpha = 0.8) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 7)) +
      labs(x = "-log10 Adjusted P-value",
           y = "GO Term",
           title = paste("Top GO Terms for", ct),
           color = "Module",
           size = "Gene Count"))
  dev.off()
}

# dot plot
pdf("plots/go_summary_dotplot.pdf")
ggplot(top_go, aes(x = -log10(p.adjust), y = Label, size = Count, color = Module)) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7)) +
  labs(x = "-log10 Adjusted P-value", y = "GO Term",
       title = "Top Enriched GO Terms by Cell Type and Module",
       color = "Module", size = "Gene Count")
dev.off()

# Example for one cell type and module
example_go <- go_summary %>% filter(Celltype == "Treg", Module == "8") %>%
  slice_min(p.adjust, n = 10)

# box plot
pdf("plots/go_summary_barplot_Tregs8.pdf")
ggplot(example_go, aes(x = reorder(Description, p.adjust), y = -log10(p.adjust), fill = Count)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 GO BP terms (Tregs Module 8)",
       x = "GO Term", y = "-log10 Adjusted p-value")
dev.off()

# enrichment map
go_enriched <- go_results_all[["Treg_8"]]
go_enriched <- pairwise_termsim(go_enriched)  # compute term similarity
pdf("plots/go_summary_enrichmentmap_Tregs8.pdf")
emapplot(go_enriched, showCategory = 30)
dev.off()

