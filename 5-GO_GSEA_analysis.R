#!/usr/bin/env Rscript
# Rowan O'Hara-Payne
# 5/28/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Gene Ontology & Gene Set Enrichment Analysis of the K-means clustering


# Read in libraries and data ----

# set library path
.libPaths("/autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq/lib/R/library")

library(dplyr)
library(tibble)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(enrichplot)

# load modules and DGE results
modules_df <- read.csv("/local/projects/idea/rohara/immunophenotype_SCZ/gene_modules_all_celltypes.csv")

dge <- read.table("/local/projects/idea/rohara/immunophenotype_SCZ/dkelly_singlets.pseudobulk_edgeR.dge.txt", sep='\t')
names(dge) <- dge[1,]
dge <- dge[-1,]

options(mc.cores = 1) # avoid forking issues


# GO & GSEA Analysis ----

# load hallmark gene sets for fgsea
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# store results
go_results_all <- list()
gsea_results_all <- list()

celltypes <- unique(modules_df$celltype)

# loop through each cell type
for (ct in celltypes) {
  mod_df <- modules_df %>% filter(celltype == ct)
  
  # GO analysis per module
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
      go_file <- paste0("/local/projects/idea/rohara/immunophenotype_SCZ/plots/", gsub(" ", "_", ct), "_module_", mod, "_GO.pdf")
      pdf(go_file)
      print(barplot(go, showCategory = 20, title = paste("GO Enrichment -", ct, "Module", mod)))
      dev.off()
    }
  
  }
  
  # GSEA analysis for each celltype x contrast
  for (contrast in unique(dge$Contrast)) {
	
	  # filter data
    	dge_sub <- dge %>%
      		filter(Celltype == ct, Contrast == contrast)

	if (nrow(dge_sub) < 100) {
      		message("Skipping ", ct, " - ", contrast, " (not enough genes)")
      		next
    	}

	# create ranked list
    	gene_ranks <- dge_sub %>%
		mutate(logFC = as.numeric(logFC)) %>%
      		group_by(Symbol) %>%
      		summarise(score = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
      		filter(!is.na(score)) %>%
		deframe()

    	gene_ranks <- sort(gene_ranks, decreasing = TRUE)
	
	common_genes <- intersect(names(gene_ranks), unique(unlist(hallmark_sets)))
    	gene_ranks <- gene_ranks[common_genes]

    	# skip if all NA or empty
    	if (length(gene_ranks) == 0 || all(is.na(gene_ranks))) {
      		message("Skipping ", ct, " - ", contrast, " (empty or all NA)")
      		next
    	}

	# run fgsea
    	gsea_res <- fgseaMultilevel(pathways = hallmark_sets,
                                stats = gene_ranks)

    	# save
    	id <- paste(gsub(" ", "_", ct), contrast, sep = "_")
    	gsea_results_all[[id]] <- gsea_res

    	# save top results to CSV
	top_gsea <- gsea_res %>% arrange(padj) %>% head(20)
    	write.csv(top_gsea,
              file = paste0("/local/projects/idea/rohara/immunophenotype_SCZ/gsea_", id, ".csv"),
              row.names = FALSE)

	# plot top pathways
    	if (nrow(top_gsea) > 0) {
      		top_gsea <- top_gsea %>%
        	mutate(pathway = str_trunc(pathway, 60),  # shorten names
               		pathway = factor(pathway, levels = rev(pathway)))  # reverse for plotting

      	pdf(paste0("/local/projects/idea/rohara/immunophenotype_SCZ/plots/gsea_dotplot_", id, ".pdf"), width = 8, height = 6)
      	print(
        	ggplot(top_gsea, aes(x = NES, y = pathway, size = -log10(padj), color = padj)) +
          	geom_point() +
          	scale_color_gradient(low = "red", high = "blue") +
          	labs(title = paste("GSEA -", ct, "-", contrast),
               		x = "Normalized Enrichment Score (NES)",
               		y = "Pathway",
               		size = "-log10(padj)",
               		color = "Adjusted P") +
          	theme_bw()
      	)
      	dev.off()
    }
  }
}

# save results
saveRDS(go_results_all, "/local/projects/idea/rohara/immunophenotype_SCZ/go_results_all.rds")
saveRDS(gsea_results_all, "/local/projects/idea/rohara/immunophenotype_SCZ/gsea_results_all.rds")


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
write.csv(go_summary, "/local/projects/idea/rohara/immunophenotype_SCZ/GO_enrichment_summary.csv", row.names = FALSE)

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
  dotplot_file <- paste0("/local/projects/idea/rohara/immunophenotype_SCZ/plots/", gsub(" ", "_", ct), "_go_summary_dotplot.pdf")

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
pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/go_summary_dotplot.pdf")
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
pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/go_summary_barplot_Tregs8.pdf")
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
pdf("/local/projects/idea/rohara/immunophenotype_SCZ/plots/go_summary_enrichmentmap_Tregs8.pdf")
emapplot(go_enriched, showCategory = 30)
dev.off()

# GSEA dot plot
for (name in names(gsea_results_all)) {
  df <- gsea_results_all[[name]] %>% arrange(padj) %>% head(20)
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))

  pdf(paste0("/local/projects/idea/rohara/immunophenotype_SCZ/plots/gsea_dotplot_", name, ".pdf"))
  print(ggplot(df, aes(x = NES, y = pathway, color = padj, size = size)) +
          geom_point() +
          theme_bw() +
          scale_color_gradient(low = "red", high = "blue") +
          labs(title = name, x = "Normalized Enrichment Score", y = "Pathway"))
  dev.off()
}
