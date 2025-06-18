# Rowan K. O'Hara-Payne
# 6/17/2025
# Decomposition of the Pan et. al data and projection onto our data

# Read in libraries and data ----
library(tidyverse)
library(Seurat)
library(edgeR)
library(CoGAPS)
library(projectR)
library(fgsea)
library(msigdbr)

pan_counts <- read_tsv("GSE243841_cleaned_counts.tsv")

# check if Pan data is normalized
pan_counts <- as.numeric(as.matrix(pan_counts))
summary(pan_counts)  # large range, not normalized

# use edgeR to normalize expression counts
dge <- DGEList(counts = as.numeric(pan_counts))

# filter lowly expressed genes (optional but recommended)
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.sizes=FALSE]

# TMM normalization
dge <- calcNormFactors(dge)

# get logCPM values (recommended for CoGAPS)
logCPM <- cpm(dge, log=TRUE, prior.count=1)


