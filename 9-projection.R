# Rowan K. O'Hara-Payne
# 7/21/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Using projectR to project the patterns from CoGAPS onto other datasets


# Load libraries and read in data----
library(Seurat)
library(CoGAPS)
library(projectR)
library(pheatmap)

# read in CoGAPS results from single cell datasets
kelly <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kellyCoGAPS_30Pat.rds")
kornberg <- readRDS("/local/projects/idea/rohara/immunophenotype_SCZ/kornbergCoGAPS_30Pat.rds")

# read in bulk datasets
pan <- read.table(file="/local/projects/idea/rohara/immunophenotype_SCZ/GSE243841_cleaned_counts.tsv", sep="\t", header=T)
glahn <- read.csv("salmon_merged_gene_counts.csv")


# Projection ----

## project single cell datasets onto bulk ----

# project Kelly patterns onto Pan data (schizophrenia immunophenotype)
kelly_pan_proj <- projectR(data=pan, loadings=kelly)

# project Kelly patterns onto Glahn data (psychosis)
kelly_glahn_proj <- projectR(data=glahn, loadings=kelly)

# project Kornberg patterns onto Pan data (schizophrenia immunophenotype)
kornberg_pan_proj <- projectR(data=pan, loadings=kornberg)

# project Kornberg patterns onto Glahn data (psychosis)
kornberg_glahn_proj <- projectR(data=glahn, loadings=kornberg)


## project single cell onto single cell ----

# project Kelly patterns onto Kornberg data (gluten intolerance)
kelly_kornberg_proj <- projectR(data=kornberg, loadings=kelly)

# project Kornberg patterns onto Kelly data (schizophrenia and gluten intolerance immunophenotype)
kornberg_kelly_proj <- projectR(data=kelly, loadings=kornberg)


# Visualization ----



