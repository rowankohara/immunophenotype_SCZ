# Rowan K O'Hara-Payne
# 6/3/2025
# conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
# Do a decomposition and projection analysis using two outside data sets


# Read in libraries and data ----
library(tidyverse)
library(Seurat)
library(CoGAPS)
library(projectR)
library(fgsea)
library(msigdbr)

# Pan et al. GEO data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243841
raw_counts <- read_tsv("GSE243841_raw_counts_all.txt", comment = "#")

# Kornberg et al. GEO data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252762
kornberg_counts <- readRDS("GSE252762_kornberg_han_celiac_dataset_fig1-2.RDS")

# Our data
kelly_counts <- readRDS("dkelly_singlets.pseudobulk.rds")


# Data Cleaning & Wrangling ----
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


