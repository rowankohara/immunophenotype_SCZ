# Rowan O'Hara-Payne
# 5/5/2025
# conda activate /local/projects-t3/idea/amentlab_software/conda_envs/Seurat5
# Based on /local/projects/idea/sament2/dkelly_tcell/dge/edgeR.R 
# Make pseudobulk, do differential gene expression analysis


# Read in libraries and data ----
library(Seurat)
library(edgeR)
library(dplyr)
library(tibble)
library(tidyr)

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


meta = read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')
rownames(meta) = meta$MPRCID
meta = meta[ MPRCID , ]

# Treg

design = model.matrix( ~ 0 + Group + Sex + Age , data = meta )

dgelist = list()

types = names(elist)
n = length(types)
for( i in 1:n ) {
  cat( types[i] , '\n' )
  counts = elist[[i]]
  y = DGEList( counts , group = meta$Group )
  keep = filterByExpr( y , min.count = 3 , min.total.count = 15 )
  y = y[ keep , ]
  y = calcNormFactors(y)
  y = estimateDisp( y , design )

  fit = glmQLFit( y )
  HvL = glmQLFTest( fit , contrast = c(0,1,-1,0,0) )
  HvC = glmQLFTest( fit , contrast = c(-1,1,0,0,0) )
  LvC = glmQLFTest( fit , contrast = c(-1,0,1,0,0) )
  HvL = as.data.frame(topTags( HvL , n = Inf ))
  HvL$Symbol = rownames(HvL)
  HvL$Contrast = 'SCZHi vs. SCZLo'
  HvC = as.data.frame(topTags( HvC , n = Inf ))
  HvC$Symbol = rownames(HvC)
  HvC$Contrast = 'SCZHi vs. HC'
  LvC = as.data.frame(topTags( LvC , n = Inf ))
  LvC$Symbol = rownames(LvC)
  LvC$Contrast = 'SCZLo vs. HC'
  res = rbind( HvL , HvC , LvC )
  res = res[ order( res$PValue ) , ]
  res$Celltype = types[i]
  dgelist[[i]] = res
}

dge = do.call( rbind , dgelist )
dge = dge[ order( dge$PValue ) , ]

write.table( dge , quote=F , sep='\t' , row.names=F ,
             file = 'dkelly_singlets.pseudobulk_edgeR.dge.txt' )


# Make plots ----

# barplot
sig = dge[ dge$PValue < 0.01 , ]
sig$dir = 'up'
sig$dir[ sig$logFC < 0 ] = 'down'

sig$cond = paste( sig$Celltype , sig$Contrast )

nsig = table( sig$cond , sig$dir )
nsig[,1] = -1 * nsig[,1]

types = sort(  unique( sig$Celltype ) )

pdf('nsig.barplot.pdf',width=8,height=5)
par( mar = c(8,4,1,1) , bty = 'l' )
barplot( height = nsig[,1] , ylim = c(-300,200) , space = c(1,0,0) ,
         xaxt = 'n' , col = c('maroon','green','blue') )
barplot( add = T , height = nsig[,2] , space = c(1,0,0) ,
         col = c('maroon','green','blue') , xaxt = 'n' )
axis( side = 1 , at = 1:13 * 4 - 1.5 , las = 2 , label = types )
abline( h = c(0,-300) )
legend( 'bottomright' , bty = 'n' , fill = c('maroon','green','blue') ,
        legend = c('SZhigh vs. HC','SZlow vs. HC','SZhigh vs. SZlow') )
mtext( side = 2 , line = 2.5 , '# of DEGs (P < 0.01)' )
dev.off()


# boxplots
meta$xpos = 1
meta$xpos[ meta$Group == 'SCZ High' ] = 2
meta$xpos[ meta$Group == 'SCZ Low' ] = 3
meta$color = 'maroon'
meta$color[ meta$Group == 'SCZ Low' ] = 'green'
meta$color[ meta$Group == 'HC' ] = 'blue'

e = elist[['Treg']]
e = cpm( e , log = T )

pdf('dkelly_pilot.treg.top_degs.boxplot.pdf',width=10,height=2.5)
# boxplot( y ~ grp , outline = F )
par( bty = 'l' , mfrow = c(1,5) , mar = c(3,2.5,2.5,0) , oma = c(1,1,0,0) )
for( gene.i in c('TNFAIP3','NR4A2','ELL2','FMN1','RAP2A') ) {
y = as.vector(t(e[gene.i,]))
plot( y = y ,
      x = jitter( meta$xpos ) ,
      col = meta$color ,
      pch = 19 , cex = 1.5 ,
      xlim = c(0.5,3.5) ,
      xaxt = 'n' ,
      xlab = '' , ylab = '' )
if( gene.i == 'TNFAIP3' ) mtext( outer = T , side = 1 ,
                              'SCZ AGA group' )
if( gene.i == 'TNFAIP3' ) mtext( outer = T , side = 2 , 'logCPM' )
boxplot( add = T , y ~ meta$xpos , outline = F , col = NA ,
         notch = F , xaxt = 'n' , border = c('blue','maroon','green') )
axis( side = 1 , at = 1:3 ,
      labels = c('HC','SZ Hi','SZ Lo') )
mtext( side = 3 , font = 4 , adj = 0 , cex = 1 , line = 1 , gene.i )
}
dev.off()


# volcano plot
custom_order <- c("0", "1", "10", "12", "14", "24-1", "25", "27", "21-0", "21-1", "21-2", "26", "29", "28", "30", "32", "24-0", "15", "19", "11", "18", "2", "9", "17", "22", "8", "13-1", "13-0", "23", "16", "6", "31", "7", "20", "5", "3", "4")
col_pal <- c("darkturquoise", "darkslategrey", "cadetblue2", "blue3", "darkcyan", "cyan", "aquamarine3", "cornflowerblue", "dodgerblue3", "deepskyblue4", "darkblue", "royalblue3", "lightblue", "slateblue1", "navy", "powderblue", "seagreen", "limegreen", "olivedrab", "springgreen4", "yellow", "violetred1", "tomato", "tan1", "salmon", "rosybrown1", "red3", "palevioletred2", "orange", "maroon", "magenta3", "coral3", "gold", "firebrick", "deeppink", "darkorange3", "brown2")
names(col_pal) <- custom_order
types = levels(obj)


dge$logP = -log10(dge$PValue)

useTypes = unique( dge$Celltype )
pal = pal[ useTypes ]

pdf('dge.volcano.alltypes.pdf')
par( bty = 'l' )
plot( x = dge$logFC ,
      y = dge$logP ,
      pch = 19 ,
      col = 'lightgrey' )
for( t in useTypes ) {
      dge.t = dge[ dge$Celltype == t &
                   dge$PValue < 0.01 , ]
      points( x = dge.t$logFC ,
              y = dge.t$logP ,
              pch = 19 ,
              col = pal[t] )
}
abline( h = 2 , v = c(-.5,.5) , lty = 2 )
legend( 'bottomright' , pch = 19 , col = pal , legend = names(pal) )
dev.off()

x = dge[ dge$Celltype == 'Treg' & dge$Contrast == 'SCZHi vs. HC' , ]
pdf('dge.volcano.treg.HvC.pdf')
par( bty = 'l' )
plot( x = x$logFC ,
      y = x$logP ,
      pch = 19 ,
      col = 'lightgrey' )
up = x[ x$PValue < 0.01 & x$logFC > 0.5 , ]
points( x = up$logFC ,
        y = up$logP ,
        col = 'maroon' ,
        pch = 19 )
sel = head( up , 5 )
text( x = sel$logFC + 0.1 ,
      y = sel$logP , adj = 0 ,
      sel$Symbol )
dn = x[ x$PValue < 0.01 & x$logFC < -0.5 , ]
points( x = dn$logFC ,
        y = dn$logP ,
        col = 'blue' ,
        pch = 19 )
sel = head( dn , 5 )
text( x = sel$logFC - 0.1 ,
      y = sel$logP , adj = 1 ,
      sel$Symbol )
abline( h = 2 , v = c(-.5,.5) , lty = 2 )
dev.off()


# Redo with l3 celltype ----

# read in Seurat object, get RNA counts, set Idents to L3 cell type
obj <- readRDS("kelly_postqc_celltypes_meta.rds")
DefaultAssay(obj) = 'RNA'
counts <- GetAssayData(obj, layer = "counts")
Idents(obj) <- obj$predicted.celltype.l3


# look at the summary counts for each cell type
meta_df <- obj@meta.data %>% 
	mutate(celltype = as.character(predicted.celltype.l3), 
	       sample = as.character(sample), 
	       group = as.character(Group),
	       sex = as.factor(Sex),
	       age = as.numeric(Age))
cell_counts <- meta_df %>% 
	group_by(celltype, sample, group) %>% 
	summarize(cells = n(), .groups = "drop")
summary_per_ct <- cell_counts %>% 
	group_by(celltype) %>% 
	summarize(total_cells = sum(cells), n_samples_total = n_distinct(sample), median_cells_per_sample = median(cells), p25_cells_per_sample = quantile(cells, .25), p75_cells_per_sample = quantile(cells, .75)) %>% 
	arrange(desc(total_cells))
print(summary_per_ct)


# filter cell types to make sure they have at least 100 total cells, 4 samples represented, and 2 samples per group
min_total_cells <- 100
min_samples_total <- 6
min_samples_per_group <- 2
ct_group_summary <- cell_counts %>%
  group_by(celltype, group) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = group, values_from = n_samples, values_fill = 0)
ct_ok <- summary_per_ct %>%
  left_join(ct_group_summary %>% 
	    group_by(celltype) %>% 
	    summarize(across(everything(), ~sum(.x))), by = "celltype")
write.csv(ct_ok, file="cell_counts_summary_per_sample.csv")
keep_cts <- summary_per_ct %>%
  filter(total_cells >= min_total_cells,
         n_samples_total >= min_samples_total) %>%
  pull(celltype)


# filter out samples for a cell type if they have fewer than 10 cells
min_cells_per_sample <- 10
samples_to_use <- cell_counts %>%
  filter(cells >= min_cells_per_sample) %>%
  group_by(celltype) %>%
  summarize(samples = list(unique(sample)), .groups = "drop")


samples_to_use_list <- setNames(samples_to_use$samples, samples_to_use$celltype)
samples_to_use_list <- samples_to_use_list[keep_cts]


make_pseudobulk_for_celltype <- function(ct, samples_keep) {
	# filter for cell type
	cells_ct <- rownames(meta_df)[meta_df$celltype == ct]
	if (length(cells_ct) == 0) return(NULL)

	# subset to cells from samples that have at least 10 cells for that cell type
	cells_keep <- cells_ct[meta_df[cells_ct, "sample"] %in% samples_keep]
	if (length(cells_keep) == 0) return(NULL)

	# subset counts to those cells
	mat <- counts[, cells_keep, drop = FALSE]  # genes x cells

	# get sample vector aligned to columns
	sample_vec <- meta_df[cells_keep, "sample"]

	# sum counts per sample, convert to simple matrix
	if (is(mat, "dgCMatrix")) mat_dense <- as.matrix(mat) else mat_dense <- mat
	pseudobulk <- t(rowsum(t(mat_dense), group = sample_vec))
	
	# order columns
	pseudobulk <- pseudobulk[, sort(colnames(pseudobulk)), drop = FALSE]
	return(pseudobulk)
}


celltypes_available <- names(samples_to_use_list)
pseudobulk_list <- list()

# filter for cell types that pass our thresholds
for (ct in celltypes_available) {
	pb <- make_pseudobulk_for_celltype(ct, samples_to_use_list[[ct]])
	if (!is.null(pb)) pseudobulk_list[[ct]] <- pb
}

sample_table <- meta_df %>%
        as.data.frame() %>%
        distinct(sample, group, age, sex) %>%
        arrange(group)
rownames(sample_table) <- sample_table$sample

# go through pseudobulk and get the sample names and groups
ct_summary <- lapply(names(pseudobulk_list), function(ct) {
			samp_names <- colnames(pseudobulk_list[[ct]])
			gr <- unique(sample_table[samp_names, "group"])
			per_group <- sample_table[samp_names, ] %>%
				as.data.frame() %>%
				group_by(group) %>%
				summarize(n = n_distinct(sample), .groups = "drop")
			tibble(celltype = ct,
				total_samples = length(samp_names),
				per_group = list(per_group))
		}) %>% bind_rows()


keep_celltypes <- c()
# filter for cell types that have two samples for each group (HC, SCZ Low, SCZ High) to make comparison easy
for (i in seq_len(nrow(ct_summary))) {
	ct <- ct_summary$celltype[i]
	total_samps <- ct_summary$total_samples[i]
	per_group_df <- ct_summary$per_group[[i]]

	n_groups_meeting <- sum(per_group_df$n >= min_samples_per_group)

	if(total_samps >= min_samples_total && n_groups_meeting >= 2) {
		keep_celltypes <- c(keep_celltypes, ct)
	}
}

pseudobulk_list <- pseudobulk_list[keep_celltypes]


results_list <- list()

for (ct in names(pseudobulk_list)) {
  message("Running edgeR for: ", ct)
  pb <- pseudobulk_list[[ct]]   # genes x samples for this ct
  sample_names <- colnames(pb)

  # use sample metadata for these samples
  samp_meta_ct <- sample_table %>% filter(sample %in% sample_names) %>%
    mutate(group = factor(group)) %>%
    arrange(sample)

  # ensure order matches pb columns
  samp_meta_ct <- samp_meta_ct[match(sample_names, samp_meta_ct$sample), ]

  # build DGEList
  y <- DGEList(counts = pb, samples = samp_meta_ct)

  # filter lowly expressed genes (keep genes expressed in at least some samples)
  keep <- filterByExpr(y, group = y$samples$group, min.count = 3, min.total.count = 15)
  y <- y[keep, ]
  y <- calcNormFactors(y)

  # design
  design <- model.matrix(~ 0 + group + sex + age, data = y$samples)
  colnames(design) <- c("HC", "SCZ_High", "SCZ_Low", "sexMale", "age")

  # fit
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  # contrasts
  groups_present <- colnames(design)[1:3]
  contrs <- list()
  if (all(c("SCZ_High","SCZ_Low") %in% groups_present)) {
    contrs$SCZHigh_vs_SCZLow <- makeContrasts(SCZ_High - SCZ_Low, levels = design)
  }
  if (all(c("SCZ_High","HC") %in% groups_present)) {
    contrs$SCZHigh_vs_HC <- makeContrasts(SCZ_High - HC, levels = design)
  }
  if (all(c("SCZ_Low","HC") %in% groups_present)) {
    contrs$SCZLow_vs_HC <- makeContrasts(SCZ_Low - HC, levels = design)
  }

  # run tests for each contrast and collect tables
  res_ct <- list()
  for (cn in names(contrs)) {
    qlf <- glmQLFTest(fit, contrast = contrs[[cn]])
    tab <- topTags(qlf, n = Inf)$table
    tab <- tab %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "gene") %>%
      mutate(contrast = cn, celltype = ct)
    res_ct[[cn]] <- tab
  }

  results_list[[ct]] <- bind_rows(res_ct)
}

# combine all results
dge_results <- bind_rows(results_list)
dge_results <- dge_results[order(dge_results$PValue), ]
write.csv(dge_results, file="kelly_l3_celltypes_edgeR_dge.csv")


