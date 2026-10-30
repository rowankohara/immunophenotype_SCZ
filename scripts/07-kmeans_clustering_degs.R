# Rowan K. O'Hara-Payne
# 5/5/2025 
# conda enviroment: rnaseq
# Kmeans clustering of the differentially expressed genes
####################################################################################################################


# Read in libraries and data ----
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# read in differential gene expression data
dge <- read.table("dkelly_Class.pseudobulk_edgeR.dge.txt", sep="\t", header = TRUE)

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

celltype_pal <- c(
  "3"    = "coral3",
  "5"    = "orangered",
  "4"    = "chocolate1",
  "16"   = "goldenrod",
  "6"    = "darkgoldenrod",
  "7"    = "sandybrown",
  "31"   = "darkorange3",
  "20"   = "indianred1",
  "13-0" = "red2",
  "23"   = "orange",
  "8"    = "rosybrown1",
  "9"    = "tomato",
  "17"   = "tan1",
  "22"   = "salmon",
  "13-1" = "palevioletred2",
  "28"   = "deeppink3",
  "2"    = "gold",
  "18"   = "yellow",
  "12"   = "dodgerblue2",
  "24-1" = "cyan",
  "25"   = "aquamarine3",
  "0"    = "darkturquoise",
  "1"    = "steelblue2",
  "10"   = "cadetblue2",
  "14"   = "darkcyan",
  "27"   = "orchid3",
  "21-0" = "dodgerblue3",
  "21-1" = "deepskyblue4",
  "21-2" = "darkblue",
  "26"   = "royalblue3",
  "11"   = "springgreen2",
  "15"   = "limegreen",
  "19"   = "olivedrab",
  "24-0" = "seagreen",
  "29"   = "mediumseagreen",
  "30"   = "navy",
  "32"   = "mediumpurple3"
)

clusters_order <- c("3", "5", "4", "16", "6", "7", "31", "20", "13-0", "23", "8", "9", "17","22", "13-1","28", "2", "18", "12","24-1", "25", "0", "1", "10", "14", "27", "21-0", "21-1", "21-2", "26", "11", "15", "19", "24-0", "29", "30", "32")


# K-means modules ----
sig.genes <- unique(dge$Symbol[dge$PValue < 0.01])
effects <- clusters_order[clusters_order %in% unique(dge$Celltype)]

nGene <- length(sig.genes)
nEffects <- length(effects)

logfc <- matrix(0, ncol = nEffects, nrow = nGene)
colnames(logfc) <- effects
rownames(logfc) <- sig.genes

for( i in 1:nEffects ) {
  dge.i = dge[ dge$Symbol %in% sig.genes & dge$Celltype == effects[i] , ]
  logfc[ dge.i$Symbol , i ] = dge.i$logFC
}

cl <- kmeans(logfc, centers = 10, iter.max = 10000, nstart = 10)
clust <- cl$cluster
centroids <- cl$centers

# get the module gene lists for GO analysis
clusters <- split(
    rownames(logfc),
    clust     
)

# find enrichment of gene ontology pathways for each module
enrich_list <- lapply(clusters, function(gene_set) {

    enrichGO(
        gene          = gene_set,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        universe      = dge$Symbol
    )
})

# grab pathway names
top_terms <- sapply(enrich_list, function(x) {

    if (nrow(x@result) == 0) {
        return("None")
    }

    x@result$Description[1]
})

# create row annotations
row_annot <- data.frame(
    pathway = top_terms
)
rownames(row_annot) <- names(top_terms)

# color palettes
pathway_pal <- c("blood vessel morphogenesis" = "#E41A1C",
		 "post-transcriptional regulation of gene expression" = "#377EB8",
		 "cellular lipid catabolic process" = "#4DAF4A",
		 "vitamin transport" = "#984EA3",
		 "purine ribonucleotide catabolic process" = "#FF7F00",
		 "regeneration" = "#FFFF33",
		 "response to bacterium" = "#A65628",
		 "cytoplasmic translation" = "#F781BF",
		 "positive regulation of type I interferon-mediated signaling pathway" = "#999999",
		 "growth hormone receptor signaling pathway via JAK-STAT" = "#FB8072")
palette <- colorRampPalette(
    c("#2166AC", "white", "#B2182B")
)(99)
breaks <- seq(-0.5, 0.5, length.out = 101)
ann_colors <- list(
  celltype = celltype_pal,
  pathway = pathway_pal
)

x <- centroids
x[ x > 0.5 ] <- 0.5
x[ x < -0.5 ] <- -0.5

labs <- parse(text = celltype_labels[colnames(x)])

# create column annotations
ann_col <- data.frame(
  celltype = colnames(x)
)
rownames(ann_col) <- colnames(x)


### Fig 3C ----
pdf("plots/dge.sig-genes.kmeans10.centroids.heatmap.pdf", width = 10, height = 6)
pheatmap(
  x,
  color = palette,
  breaks = breaks,
  cellwidth = 20,
  cellheight = 20,
  cluster_cols = TRUE,
  scale = "none",
  annotation_row = row_annot,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  labels_col = labs,
  annotation_legend = FALSE
)
dev.off()

MM = cor( t(logfc) , t(centroids) )


### Supplemental Figure 6 ----
pdf("plots/dge.sig-genes.kmeans10.top10.heatmap.pdf", width = 10, height = 6)
for( i in 1:10 ) {
  pattern = i
  sorted = MM[ order( MM[,pattern] , decreasing =T ) , ]
  top10 = head( rownames( sorted ) , 10 )
  x = logfc[ top10 , ]
  x[ x > 0.5 ] <- 0.5
  x[ x < -0.5 ] <- -0.5
  pheatmap(
  	x,
	main = paste( 'Module' , i ),
  	color = palette,
  	breaks = breaks,
	cellwidth = 20,
	cellheight = 20,
	cluster_rows = FALSE,
  	cluster_cols = TRUE,
  	scale = "none",
  	annotation_col = ann_col,
  	annotation_colors = ann_colors,
  	labels_col = labs,
  	annotation_legend = FALSE)
}
dev.off()


colnames(MM) <- paste( 'kM' , 1:10 , sep = '' )
df <- data.frame(
  gene = sig.genes ,
  clust ,
  MM ,
  logfc )

### Supplemental Data Table: Kmeans_Clustering ----
write.table( df , quote=F , sep='\t' , row.names=F ,
             file = 'kmeans.k10.module_membership.txt' )

saveRDS( cl , file = 'kmeans.k10.rds' )
