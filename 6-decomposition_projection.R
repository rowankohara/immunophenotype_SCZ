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

library(tidyverse)
#library(Seurat)
#library(CoGAPS)
#library(projectR)
#library(fgsea)
#library(msigdbr)

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

# Pan et al. GEO data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243841
raw_counts <- read_tsv("GSE243841_raw_counts_all.txt", comment = "#")

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

# Pan data
pan_counts <- raw_counts %>%
  select(-Chr, -Start, -End, -Strand, -Length)

# clean sample column names (remove "./BAM/" and "_sorted.bam")
clean_names <- colnames(pan_counts)[-1] %>%
  str_replace_all("\\./BAM/", "") %>%
  str_replace_all("_sorted.bam", "")
colnames(pan_counts)[-1] <- clean_names

# remove rows where Geneid has merged info (i.e., rows with ";" in Geneid)
pan_counts <- pan_counts %>%
  filter(!str_detect(Geneid, ";"))

# save cleaned matrix (optional)
write_tsv(pan_counts, "GSE243841_cleaned_counts.tsv")

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

# summarize to one row per sample
sample_metadata <- kelly_counts@meta.data %>%
  select(sample, Group, Sex, Age) %>%
  distinct()


# Integration ----
kellyCD <- SplitObject(object = kelly_counts , split.by = "sample")
kornbergCD <- SplitObject(object = kornberg_counts , split.by = "sample")
combDat <- c(kellyCD,kornbergCD)

IntName <- 'Kelly_Kornberg'
PrintColCat <- c('sample','celltype','condition')
PrintColCont <- c('nCount_RNA','nFeature_RNA','percent.mt')


for(i in 1:length(combDat)){
	if(i==1){
    		cellCounts = ncol(combDat[[i]])
	} else {
    		cellCounts = c(cellCounts,ncol(combDat[[i]]))
	}	
}


options(future.globals.maxSize = 240000 * 1024^2)

if(length(which(cellCounts<=100))>0){
    combDat <- combDat[-which(cellCounts<=100)]
} 

for (i in 1:length(combDat)) {
  combDat[[i]] <- SCTransform(combDat[[i]], verbose = FALSE)
}
for (i in 1:length(combDat)) {
  DefaultAssay(combDat[[i]]) <- "SCT"
}

features <- SelectIntegrationFeatures(object.list = combDat, nfeatures = 3000)
combDat <- PrepSCTIntegration(object.list = combDat, anchor.features = features)

All.anchors <- FindIntegrationAnchors(object.list = combDat, normalization.method = "SCT",anchor.features = features,reduction='cca') #,reference=refInd #,reference=ctrlInd

All.integrated <- IntegrateData(anchorset = All.anchors, normalization.method = "SCT") 


All.integrated <- ScaleData(All.integrated, features = rownames(All.integrated))
All.integrated <- RunPCA(All.integrated, verbose = FALSE)
All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:30,reduction.name='umap.integrated')
All.integrated <- FindNeighbors(All.integrated, dims = 1:30, verbose = FALSE)
All.integrated <- FindClusters(All.integrated, verbose = FALSE, resolution=0.5)

DefaultAssay(All.integrated) <- "RNA"

All.integrated[['RNA']] = JoinLayers(All.integrated[['RNA']])
All.integrated <- NormalizeData(All.integrated, assay = "RNA") ## adds data
All.integrated <- ScaleData(All.integrated, features = rownames(All.integrated), assay = "RNA")

saveRDS(All.integrated,file=paste(IntName,'_integrated_V5.rds',sep=''),compress=TRUE)


# plots
subplots(All.integrated,colName='celltype',filename=paste0('plots/', IntName, 'integrated_celltype.pdf'),raster=TRUE,reduction='umap.integrated')


subplots(All.integrated,colName='study',filename=paste0('plots/', IntName, 'integrated_study.pdf'),raster=TRUE,reduction='umap.integrated')

pdf(file=paste("plots/",IntName,"_Check_Clustering_V5.pdf",sep=''),width=12,height=8)
for(k in PrintColCat){
	print(DimPlot(All.integrated, reduction = "umap.integrated", group.by = k, label = TRUE, repel = TRUE,order=TRUE) + NoLegend() )
}
for(m in PrintColCont){
	print(FeaturePlot(All.integrated, features = m, reduction = "umap.integrated", pt.size = 0.2))
}
dev.off()




