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

obj <- readRDS('dkelly_singlets.postqc.azimuth_celltypes.rds')



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
    expr[,j] = rowSums( GetAssayData( sub.j , slot = 'counts' ) )
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
pal = DiscretePalette( 24 , 'polychrome' )
types = levels(obj)
names(pal) = types

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

