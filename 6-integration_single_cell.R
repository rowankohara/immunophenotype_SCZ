# Rowan K O'Hara-Payne
# 6/3/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# conda activate /autofs/projects-t3/idea/amentlab_software/conda_envs/Seurat5_R4_3_0
# Do a decomposition and projection analysis using two outside data sets


# Read in libraries and data ----
library(Seurat)
library(future)
library(SeuratDisk)
library(glmGamPoi)
library(harmony)
library(tidyverse)

# plot function
subplots <- function(sobj,colName,filename,raster=FALSE,legend=TRUE,reduction){
	classes <- unique(sobj@meta.data[colName])[,1]
	pdf(filename,width=12,height=8)
	for(i in classes){
		sobj@meta.data[i] <- "Other"
		sobj@meta.data[which(sobj@meta.data[colName]==i),i] <- i
		if(legend==FALSE){
		print(DimPlot(sobj, reduction = reduction, group.by = i,order=c(i,'other'),raster=raster) + ggplot2::labs(title=i) + NoLegend())
    		} else {
 		print(DimPlot(sobj, reduction = reduction, group.by = i,order=c(i,'other'),raster=raster) + ggplot2::labs(title=i))
    		}
	}
	dev.off()
}

# Kornberg et al. GEO data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252762
kornberg_counts <- readRDS("GSE252762_kornberg_han_celiac_dataset_fig1-2.RDS")
kornberg_counts$study <- 'Kornberg'

# Our data
kelly_counts <- readRDS("dkelly_singlets.postqc.azimuth_celltypes.rds")
kelly_counts$study <- 'Kelly'
kelly_counts$celltype <- kelly_counts$predicted.celltype.l2

meta <- read.csv('/local/projects/idea/sament2/dkelly_tcell/dge/sample_metadata.csv')


# Data Cleaning & Wrangling ----

# Kelly data
# create an empty data frame to store cell-level metadata
cell_meta <- data.frame(row.names = colnames(kelly_counts))

# loop through each metadata column (excluding MPRCID)
meta_filtered <- meta %>% select(MPRCID, Sex, Group, Age)
for (col in colnames(meta_filtered[-1])) {
  vec <- rep(NA, ncol(kelly_counts))
  names(vec) <- colnames(kelly_counts)
  
  for (i in 1:nrow(meta_filtered)) {
    sample_id <- meta_filtered$MPRCID[i]
    
    # find cells whose sample value starts with the sample_id
    matched_idx <- grep(paste0("^", sample_id), kelly_counts@meta.data$sample)
    matched_cells <- colnames(kelly_counts)[matched_idx]
    
    # assign sample-level metadata to those cells
    vec[matched_cells] <- meta_filtered[[col]][i]
  }
  
  cell_meta[[col]] <- vec
}

# add the metadata to the Seurat object
kelly_counts <- Seurat::AddMetaData(kelly_counts, metadata = cell_meta)
kelly_counts$sample_sex <- kelly_counts$Sex
kelly_counts$condition <- kelly_counts$Group
kelly_counts$sample_age <- kelly_counts$Age
kelly_counts@meta.data$sample_sex[kelly_counts@meta.data$sample_sex == "Female"] <- "F"
kelly_counts@meta.data$sample_sex[kelly_counts@meta.data$sample_sex == "Male"] <- "M"
kelly_counts$clusters <- kelly_counts$seurat_clusters

# summarize to one row per sample
sample_metadata <- kelly_counts@meta.data %>%
  select(sample, Group, Sex, Age) %>%
  distinct()


# Harmony Integration ----
#combDat <- merge(kelly_counts, kornberg_counts, add.cell.ids = c("kelly", "kornberg"))
#combDat <- NormalizeData(combDat)
#combDat <- FindVariableFeatures(combDat)
#combDat <- ScaleData(combDat)
#combDat <- RunPCA(combDat)
#combDat <- RunHarmony(combDat, group.by.vars = "study")
#combDat <- RunUMAP(combDat, reduction = "harmony", dims = 1:30)
#combDat <- FindNeighbors(combDat, reduction = "harmony", dims = 1:30)
#combDat <- FindClusters(combDat, resolution = 0.5)

IntName <- 'Kelly_Kornberg'
PrintColCat <- c('sample','celltype','condition')
PrintColCont <- c('nCount_RNA','nFeature_RNA','percent.mt')

#pdf(file=paste("plots/",IntName,"_study.pdf",sep=''))
#DimPlot(combDat, group.by = "study")         # view study of origin
#dev.off()

#pdf(file=paste("plots/",IntName,"_clusters.pdf",sep=''))
#DimPlot(combDat, label = T, group.by = "seurat_clusters")  # view clusters
#dev.off()

#pdf(file=paste("plots/",IntName,"_celltype.pdf",sep=''), width=12,height=8)
#DimPlot(combDat, label = T, group.by = "celltype")  # view celltype
#dev.off()

#pdf(file=paste("plots/",IntName,"_nCount_RNA.pdf",sep=''))
#VlnPlot(obj, group.by = 'sample', split.by = 'doublet_finder',
	#features = 'nCount_RNA',
	#combine = FALSE, pt.size = 0)  # nCount_RNA
#dev.off()


#DefaultAssay(combDat) <- "RNA"


#saveRDS(combDat, file = paste(IntName,'_integrated_harmony_V5.rds', sep=''), compress=TRUE)



# Integration ----


options(future.globals.maxSize = 240000 * 1024^2)

features <- SelectIntegrationFeatures(object.list = combDat, nfeatures = 2000)
combDat <- PrepSCTIntegration(object.list = combDat, anchor.features = features)

All.anchors <- FindIntegrationAnchors(object.list = combDat, normalization.method = "SCT",anchor.features = features,reduction='rpca', k.anchor=15) #,reference=refInd #,reference=ctrlInd
All.integrated <- IntegrateData(anchorset = All.anchors, normalization.method = "SCT", k.weight=20)

#saveRDS(All.integrated,file=paste('./SeuratObj/',IntName,'_integrated_V5.rds',sep=''),compress=TRUE)


All.integrated <- ScaleData(All.integrated, features = rownames(All.integrated))
All.integrated <- RunPCA(All.integrated, verbose = FALSE)
All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:30,reduction.name='umap.integrated'
All.integrated <- FindNeighbors(All.integrated, dims = 1:30, verbose = FALSE)
All.integrated <- FindClusters(All.integrated, verbose = FALSE, resolution=0.5)

DefaultAssay(All.integrated) <- "RNA"

All.integrated[['RNA']] = JoinLayers(All.integrated[['RNA']])
All.integrated <- NormalizeData(All.integrated, assay = "RNA") ## adds data
All.integrated <- ScaleData(All.integrated, features = rownames(All.integrated), assay = "RNA")

saveRDS(All.integrated, 'Kelly_Kornberg_integrated_V5.rds', compress=TRUE)



