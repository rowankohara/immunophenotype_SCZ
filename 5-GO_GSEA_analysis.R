#!/usr/bin/env Rscript
# Rowan O'Hara-Payne
# 5/28/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Gene Ontology & Gene Set Enrichment Analysis of the K-means clustering


# Read in libraries and data ----

# set library path
.libPaths("/autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq/lib/R/library")

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(BiocParallel)

# load DGE results
dge <- read.table("/local/projects/idea/rohara/immunophenotype_SCZ/dkelly_singlets.pseudobulk_edgeR.dge.txt", sep='\t', header = TRUE)


# GO & GSEA Analysis ----

# load gene sets for fgsea
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

immune_sets <- msigdbr(species = "Homo sapiens", category = "C7") %>%
  split(x = .$gene_symbol, f = .$gs_name)

go_sets <- msigdbr(species = "Homo sapiens", category = "C5") %>%
  split(x = .$gene_symbol, f = .$gs_name)


celltypes <- unique(dge$Celltype)

## Gene Ontology Enrichment ----
for (ct in celltypes) {
	dge_sub <- dge %>% filter(Celltype == ct & Contrast == "SCZHi vs. SCZLo")
	
	genes_up <- dge_sub %>%
                  filter(logFC > 0) %>%
		  arrange(PValue) %>%
		  pull(Symbol)

        genes_down <- dge_sub %>%
                  filter(logFC < 0) %>%
		  arrange(PValue) %>%
		  pull(Symbol)


	eg_up <- bitr(genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
	
	eg_down <- bitr(genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


	# GO enrichment
	go_up <- enrichGO(gene = eg_up$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
	 
	go_down <- enrichGO(gene = eg_down$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

	pdf(paste0("plots/GO_barplot_", ct, "HivLo.pdf"))
	print(barplot(go_up, showCategory = 20, title = paste0("GO Enrichment Up: ", ct)))
	print(barplot(go_down, showCategory = 20, title = paste0("GO Enrichment Down: ", ct)))
	dev.off()
}


## Gene Set Enrichment Analysis ----
for (ct in celltypes) {
	dge_sub <- dge %>% filter(Contrast == "SCZHi vs. SCZLo" & Celltype == ct)
	
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

        # run fgsea
        gsea_res <- fgseaMultilevel(pathways = hallmark_sets,
                                    stats = gene_ranks,
                                    minSize = 15,
                                    maxSize = 500,
				    nPermSimple = 10000,
                                    BPPARAM = SerialParam())

        # save
        id <- paste(ct, "HivLo", "hallmark", sep = "_")

        # save top results to CSV
        top_gsea <- gsea_res %>% arrange(padj) %>% head(20) %>% dplyr::select(-leadingEdge)
        write.csv(top_gsea,
              file = paste0("results/gsea_", id, ".csv"),
              row.names = FALSE)

        # reformat pathway names
	clean_names <- function(x) {
		x %>%
			str_remove("^[A-Z]+_") %>%
			str_replace_all("_", " ") %>%
			str_to_lower() %>%
			str_to_title() %>%
			str_replace_all("\\b(Tnfa|Nfkb|Dna|Rna|Myc|Mtorc1|Uv|G2m|Il6|Jak|Stat3)\\b", toupper)
	}

	top_gsea <- top_gsea %>%
		mutate(pathway = clean_names(pathway),
		       pathway = str_trunc(pathway, 60),
		       pathway = factor(pathway, levels = rev(pathway)))

	# plot top pathways
        pdf(paste0("plots/gsea_dotplot_", id, ".pdf"), width = 12, height = 8)
        print(
                ggplot(top_gsea, aes(x = NES, y = pathway, color = -log10(padj), size = size)) +
                geom_point() +
                scale_color_gradient(low = "red", high = "blue") +
                labs(title = paste("GSEA:", ct, "in", "AGA IgG+ vs. AGA IgG-", sep = " "),
                        x = "Normalized Enrichment Score (NES)",
                        y = "Pathway",
                        color = "-log10(padj)",
                        size = "pathway size") +
                theme_bw() +
		theme(axis.text = element_text(size = 10),
		      axis.title = element_text(size = 12))
        )
        dev.off()
}



for (ct in celltypes) {
        dge_sub <- dge %>% filter(Celltype == ct)

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

        # run fgsea
        gsea_res <- fgseaMultilevel(pathways = hallmark_sets,
                                    stats = gene_ranks,
                                    minSize = 15,
                                    maxSize = 500,
                                    nPermSimple = 10000,
                                    BPPARAM = SerialParam())

        # save
        id <- paste(ct, "SCZvHC", "hallmark", sep = "_")

        # save top results to CSV
        top_gsea <- gsea_res %>% arrange(padj) %>% head(20) %>% dplyr::select(-leadingEdge)
        write.csv(top_gsea,
              file = paste0("results/gsea_", id, ".csv"),
              row.names = FALSE)

        # reformat pathway names
        clean_names <- function(x) {
                x %>%
                        str_remove("^[A-Z]+_") %>%
                        str_replace_all("_", " ") %>%
                        str_to_lower() %>%
                        str_to_title() %>%
                        str_replace_all("\\b(Tnfa|Nfkb|Dna|Rna|Myc|Mtorc1|Uv|G2m|Il6|Jak|Stat3|Il2|Stat5)\\b", toupper)
        }

        top_gsea <- top_gsea %>%
                mutate(pathway = clean_names(pathway),
                       pathway = str_trunc(pathway, 60),
                       pathway = factor(pathway, levels = rev(pathway)))

        # plot top pathways
        pdf(paste0("plots/gsea_dotplot_", id, ".pdf"), width = 12, height = 8)
        print(
                ggplot(top_gsea, aes(x = NES, y = pathway, color = -log10(padj), size = size)) +
                geom_point() +
                scale_color_gradient(low = "red", high = "blue") +
                labs(title = paste("GSEA:", ct, "in", "SCZ vs. Controls", sep = " "),
                        x = "Normalized Enrichment Score (NES)",
                        y = "Pathway",
                        color = "-log10(padj)",
                        size = "pathway size") +
                theme_bw() +
                theme(axis.text = element_text(size = 10),
                      axis.title = element_text(size = 12))
        )
        dev.off()
}

