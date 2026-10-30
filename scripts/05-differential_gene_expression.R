# Rowan O'Hara-Payne
# 5/5/2025
# conda environment: Seurat5
# Based on /local/projects/idea/sament2/dkelly_tcell/dge/edgeR.R 
# Pseudobulk expression data, then do differential gene expression analysis using edgeR
####################################################################################################################


# Read in libraries and data ----
library(Seurat)
library(edgeR)
library(tidyverse)
library(ggrepel)

obj <- readRDS("kelly_postqc_clusters_meta.rds")


# Make pseudobulk ----
DefaultAssay(obj) = 'RNA'

freq = table( Idents(obj) )
types = names(freq)[ freq >= 500 ]

elist = list()

sample = unique( obj$sample )

for( i in 1:length(types) ) {
  cat( types[i] , '\n' )
  expr = matrix( 0 , ncol = length(sample) , nrow = nrow(obj) )
  sub = subset( obj , idents = types[i] )
  for( j in 1:length(sample) ) {
    idx = which( sub$sample == sample[j] )
    if(length(idx)==0) next
    sub.j = subset( sub , cells = idx )
    expr[,j] = rowSums( GetAssayData( sub.j , layer = 'counts' ) )
    rownames(expr) = rownames(obj)
    colnames(expr) = sample
  }
  elist[[i]] = expr
}

names(elist) = types

saveRDS( elist , file = 'dkelly_singlets.pseudobulk.rds' )


# Differential gene expression analysis ----
sample = colnames(elist[[11]])
MPRCID = gsub( '\\_(.*)' , '' , sample )
MPRCID = gsub( '\\-(.*)' , '' , MPRCID )

# read in meta data
meta = read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
rownames(meta) = meta$MPRCID
meta = meta[ MPRCID , ]


## plots ----

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

dge$logP = -log10(dge$PValue)


dge <- dge %>%
       mutate(Sig = if_else(PValue < 0.01 & logFC > 1, "UP",
			    if_else(PValue < 0.01 & logFC < -1, "DOWN", "NS")))

dge$Celltype <- factor(dge$Celltype, levels = clusters_order)

meta <- meta %>% mutate(Class = if_else(Group == "SCZ High", "AGA+", "AGA-"))
meta$Class <- factor(meta$Class, levels = c("AGA-", "AGA+"))

design = model.matrix( ~ Class + Sex + Age , data = meta )

dgelist = list()

# dge using AGA status contrast
types = names(elist)
n = length(types)
for( i in 1:n ) {
  cat( types[i] , '\n' )
  counts = elist[[i]]
  y = DGEList( counts , group = meta$Class )
  keep = filterByExpr( y , min.count = 3 , min.total.count = 15 )
  y = y[ keep , ]
  y = calcNormFactors(y)
  y = estimateDisp( y , design )

  fit = glmQLFit( y )
  res = glmQLFTest( fit , coef = "ClassAGA+" )
  res = as.data.frame(topTags( res , n = Inf ))
  res$Symbol = rownames(res)
  res$Contrast = 'AGA+ v. AGA-'
  res = res[ order( res$PValue ) , ]
  res$Celltype = types[i]
  dgelist[[i]] = res
}

dge = do.call( rbind , dgelist )
dge = dge[ order( dge$PValue ) , ]

### Supplemental Data Table: DGE ----
write.table( dge , quote=F , sep='\t' , row.names=F ,
             file = 'dkelly_Class.pseudobulk_edgeR.dge.txt' )

dge$logP = -log10(dge$PValue)


dge <- dge %>%
       mutate(Sig = if_else(PValue < 0.01 & logFC > 1, "UP",
                            if_else(PValue < 0.01 & logFC < -1, "DOWN", "NS")))

dge$Celltype <- factor(dge$Celltype, levels = clusters_order)


### Fig 3A ----
pdf("plots/AGA_v_AGA_volcano_dge.pdf", width = 12, height = 8)
ggplot(dge, aes(x = logFC, y = logP)) +
  geom_point(color = "lightgray",
             size = 1.5) +
  geom_point(data = dge %>% filter(Sig != "NS"),
             aes(color = Celltype),
             size = 1.5) +
  geom_text_repel(data = dge %>% filter(FDR < 0.05, abs(logFC) > 1),
                  aes(label = Symbol),
                  size = 4,
                  color = "black",
                  max.overlaps = 30) +
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed", 
             linewidth = 0.75) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             linewidth = 0.75) +
  geom_vline(xintercept = -1,
             linetype = "dashed", 
             linewidth = 0.75) +
  scale_color_manual(values = celltype_pal,
                     labels = function(x) parse(text = celltype_labels[x]),
                     name = "cell types") +
  labs(title = parse(text = "Differential~Gene~Expression~Between~AGA^'+'~and~AGA^'-'"),
       x = "logFC",
       y = parse(text = "-log[10]*(PValue)")) +
  theme_classic() +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05))) +
  theme(plot.title = element_text(size = 17,
                                  hjust = 0.5),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 10),
		    axis.line = element_line(linewidth = 1.25, 
		                             color = "black"),
		    axis.ticks = element_line(linewidth = 1.25))
dev.off()


# per cell type
pdf("plots/celltypes_AGA_v_AGA_volcano_dge.pdf", width = 12, height = 8)
for (ct in unique(dge$Celltype)) {
        p <- ggplot(dge %>% filter(Celltype == ct), aes(x = logFC, y = -log10(PValue), color = Sig)) +
            geom_point(size = 1.5) +
            geom_text_repel(data = dge %>% filter(Celltype == ct & PValue < 0.01 & abs(logFC) > 1),
                            aes(label = Symbol),
                            size = 4,
                            color = "black",
                            max.overlaps = 15) +
            geom_hline(yintercept = -log10(0.01), 
                       linetype = "dashed", 
                       linewidth = 1) +
            geom_vline(xintercept = 1, 
                       linetype = "dashed", 
                       linewidth = 1) +
            geom_vline(xintercept = -1, 
                       linetype = "dashed", 
                       linewidth = 1) +
            scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "gray70"),
                               guide = "none") +
            labs(title = parse(text = celltype_labels[ct]),
                 x = "logFC",
                 y = parse(text = "-log[10]*(PValue)")) +
            theme_classic() +
	    scale_y_continuous(limits = c(0, NA),
                               expand = expansion(mult = c(0, 0.05))) +
            theme(plot.title = element_text(size = 15,
                                       hjust = 0.5),
                  axis.title = element_text(size = 17),
                  axis.text = element_text(size = 15),
		  axis.line = element_line(linewidth = 1.25, color = "black"),
                  axis.ticks = element_line(linewidth = 1.25))
        print(p)
}
dev.off()


