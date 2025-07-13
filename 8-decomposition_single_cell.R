# Rowan K. O'Hara-Payne
# 7/12/2025
# conda activate /autofs/projects-t3/idea/amentlab_software/conda_envs/Seurat5_R4_3_0
# Decomposition of our data to be projected onto the bulk data Pan et. al and Glahn et. al


# Read in libraries and data ----
library(CoGAPS)
library(Seurat)

combDat <- readRDS("/autofs/projects-t2/idea/rohara/immunophenotype_SCZ/Kelly_Kornberg_integrated_harmony_V5.rds")

DefaultAssay(combDat) <- "RNA"


# Find variable genes ----

combDat <- FindVariableFeatures(combDat, selection.method = "vst", nfeatures = 3000)
hvg <- VariableFeatures(combDat)

expr_mat <- GetAssayData(combDat[["RNA"]], layer = "data")  # log-normalized
expr_mat <- expr_mat[hvg, ]                          # subset to HVGs
expr_mat <- as.matrix(expr_mat)

# shift if needed
if (min(expr_mat) < 0) {
  expr_mat <- expr_mat - min(expr_mat) + 1e-6
}


# Run CoGAPS ----

nThrd=16
cells <- colnames(expr_mat)
genes <- rownames(expr_mat)
params <- new("CogapsParams")
params <- setParam(params, "sparseOptimization", FALSE)
params <- setParam(params, "nIterations", 10000)
params <- setParam(params, "nPatterns", 30)
params <- setDistributedParams(params, nSets=1)

xxCoGAPS <- CoGAPS(data=expr_mat, params=params, geneNames=genes, sampleNames=cells, messages=TRUE, transposeData=FALSE, nPatterns=30, nIterations=10000, nThreads=16)

save(xxCoGAPS, '/autofs/projects-t2/idea/rohara/immunophenotype_SCZ/Kelly_Kornberg_CoGAPS_30Pat.rda')
