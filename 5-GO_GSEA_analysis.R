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
dge <- read.table("/local/projects/idea/rohara/immunophenotype_SCZ/dkelly_Class.pseudobulk_edgeR.dge.txt", sep='\t', header = TRUE)

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

all_go <- list()

for (ct in celltypes) {
	dge_sub <- dge %>% filter(Celltype == ct)
	
	genes_up <- dge_sub %>%
                  filter(PValue < 0.1, logFC > 0) %>%
		  arrange(PValue)

        genes_down <- dge_sub %>%
                  filter(PValue < 0.1, logFC < 0) %>%
		  arrange(PValue)


	# GO enrichment
	go_up <- enrichGO(gene = genes_up$Symbol,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   universe = dge$Symbol)
	
	#if (nrow(go_up) > 1) {
		#go_up <- as.data.frame(go_up)
		#go_up$Celltype <- ct
		#all_go[[ct]] <- go_up
	#}

	
	go_down <- enrichGO(gene = genes_down$Symbol,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
         	   universe = dge$Symbol)

	if (nrow(go_down) > 1) {
		go_down <- as.data.frame(go_down)
		go_down$Celltype <- unique(dge_sub$Celltype)
		all_go[[ct]] <- go_down
	}

	
}

# combine all cell types
go_df <- bind_rows(all_go)

# reorder pathways globally
go_df <- go_df %>%
        mutate(
                Description = str_trunc(Description, 50),
                Description = factor(
                        Description,
                        levels = rev(unique(Description))
                )
        )

go_df$Celltype <- factor(
    go_df$Celltype,
    levels = c("23", "17", "13-1", "18", "16", "13-0", "2", "6", "7", "3", "4", "15", "20", "11", "5", "22", "9", "8", "19", "10", "21-0", "0", "1", "12", "14")
)

go_df2 <- go_df %>%
		group_by(Celltype) %>%
                arrange(p.adjust) %>%
                slice_head(n = 5) %>%
                ungroup()
                

# facet plot
p <- ggplot(
        go_df2,
        aes(
                x = Celltype,
                y = Description,
                color = -log10(p.adjust),
                size = Count
        )
) +
        geom_point() +
        scale_x_discrete(
                labels = function(x) parse(text = celltype_labels[x])
        ) +
        scale_size_continuous(range = c(2, 10)) +
        theme_bw() +
        labs(
                title = "",
                x = NULL,
                y = "Gene Ontology Biological Process",
                color = "-log10(p.adjust)",
                size = "Count"
        ) +
        theme(axis.text.x = element_text(
                        angle = -90,
                        vjust = 0,
                        hjust = 0,
                        size = 15,
                        color = "black"
                ),

                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 15,
                                           color = "black")
        )

ggsave(
        "plots/go_all_celltypes.pdf",
        p,
        width = 18,
        height = 18
)

write.csv(go_df, "go_results_all_hallmark.csv")


## Gene Set Enrichment Analysis ----

clusters_order <- c("3", "5", "4", "16", "6", "7", "31", "20", "13-0", "23", "8", "9", "17","22", "13-1","28", "2", "18", "12","24-1", "25", "0", "1", "10", "14", "27", "21-0", "21-1", "21-2", "26", "11", "15", "19", "24-0", "29", "30", "32")

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

all_gsea <- list()

clean_names <- function(x) {
        x %>%
                str_remove("^[A-Z]+_") %>%
                str_replace_all("_", " ") %>%
                str_to_lower() %>%
                str_to_title() %>%
                str_replace_all(
                        "\\b(Tnfa|Nfkb|Dna|Rna|Myc|Mtorc1|Uv|G2m|Il6|Jak|Stat3|Il2|Stat5|Tgf)\\b",
                        toupper
                )
}


for (ct in celltypes) {

        dge_sub <- dge %>%
                filter(Celltype == ct)

        # ranked gene list
        gene_ranks <- dge_sub %>% 
		mutate(score = -log10(PValue) * sign(logFC)) %>%
		dplyr::select(Symbol, score) %>%
		deframe()

        gene_ranks <- sort(gene_ranks, decreasing = TRUE)

        common_genes <- intersect(
                names(gene_ranks),
                unique(unlist(hallmark_sets))
        )

        gene_ranks <- gene_ranks[common_genes]

        # fgsea
        gsea_res <- fgseaMultilevel(
                pathways = hallmark_sets,
                stats = gene_ranks,
                minSize = 15,
                maxSize = 500,
                nPermSimple = 10000,
                BPPARAM = SerialParam()
        )

        # keep top pathways
        top_gsea <- gsea_res %>%
                arrange(padj) %>%
                slice_head(n = 20) %>%
                mutate(
                        Celltype = ct,
                        pathway = clean_names(pathway)
                )

        all_gsea[[ct]] <- top_gsea
}

# combine all cell types
gsea_df <- bind_rows(all_gsea)

# reorder pathways globally
gsea_df <- gsea_df %>%
        mutate(
                pathway = str_trunc(pathway, 50),
                pathway = factor(
                        pathway,
                        levels = rev(unique(pathway))
                )
        )


gsea_df$Celltype <- factor(
    gsea_df$Celltype,
    levels = c("23", "17", "13-1", "18", "16", "13-0", "2", "6", "7", "3", "4", "15", "20", "11", "5", "22", "9", "8", "19", "10", "21-0", "0", "1", "12", "14")
)


# facet plot
p <- ggplot(
        gsea_df,
        aes(
                x = Celltype,
                y = pathway,
                color = NES,
                size = -log10(padj)
        )
) +
        geom_point() +
        scale_color_gradient2(low = "#2166AC",
			      mid = "white",
			      high = "#B2182B",
			      midpoint = 0
) +
	scale_x_discrete(
                labels = function(x) parse(text = celltype_labels[x])
        ) +
	scale_size_continuous(range = c(2, 10),
			      breaks = c(0.5, 2, 4, 8, 16)) +
        theme_bw() +
        labs(
                title = "",
                x = NULL,
                y = "Pathway",
                color = "NES",
                size = "-log10(padj)"
        ) +
        theme(axis.text.x = element_text(
                        angle = -90,
                        vjust = 0,
                        hjust = 0,
                        size = 15,
			color = "black"
                ),

                axis.ticks.x = element_blank(),
                axis.text.y = element_text(size = 15,
					   color = "black")
        )

ggsave(
        "plots/gsea_facetwrap_all_celltypes.pdf",
        p,
        width = 18,
        height = 14
)

gsea_df %>%
	dplyr::select(-leadingEdge) %>%
	write.csv("gsea_results_all_hallmark.csv")
