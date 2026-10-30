# Molecular adaptations of activated T-cells in an inflammation-associated schizophrenia subgroup
**Deepak Salem, Rowan K. O’Hara-Payne, Sarah M. Clark, Marcia Cortes-Gutierrez, Nevil J. Singh, Daniel J. Roche, Deanna L. Kelly, Seth A. Ament**

Thirty-five percent of people with schizophrenia-related disorders (SRD) form a high-inflammation subgroup defined by elevated anti-gliadin antibodies (AGA+) and inflammatory proteins and associated with an increased severity of negative symptoms. 
However, the immune mechanisms mediating these effects remain poorly defined. 
Here, we characterized transcriptional signatures of peripheral immune cells in AGA+ SRD (n=7) compared to AGA-negative (AGA-) SRD (n=3) and healthy controls (HCs; n=5), using single-cell RNA-sequencing (scRNA-seq) of peripheral blood mononuclear cells (PBMCs; n=61,001). 
AGA+ SRD was associated with increased abundance of T-helper-17 cells (Th17), T-follicular helper-1 (Tfh1), CD5+ B cells, plasmacytoid dendritic cells (pDCs), and several CD8+ T cell subsets, including memory and Natural Killer-T-like activated subsets. 
In contrast, AGA- SRD showed higher abundances of several monocyte subsets than either AGA+ SRD or HCs. 
Differential gene expression and pathway analyses identified altered metabolic and inflammatory signaling, including JAK/STAT, IL-6, WNT/β-catenin, and type I interferon-associated pathways across distinct immune cell populations. 
Cell-cell communication analysis further identified altered signaling between Th17 and myeloid populations, alongside differential myeloid signaling states and predicted LGALS9-mediated interactions between C1Qhi CD16+ monocytes and T cell populations in AGA+ SRD. 
Collectively, these results define a unique peripheral immune signature associated with AGA status in SRD, characterized by adaptive immune cell enrichment and altered immune cell communication in AGA+ SRD and a more prominent monocyte signature in AGA- SRD. 
These findings contribute to our current understanding of negative symptom-predominant, AGA+ SRD and highlight potential targets for therapeutic intervention for this subgroup.

Keywords: *schizophrenia*, *single cell*, *T cell*, *antigliadin antibody*, *RNA sequencing*, *inflammation*.

[Pre-print](https://www.biorxiv.org/content/10.64898/2026.05.22.727194v1)

## Scripts
1. [01-qc_and_preprocessing.R](scripts/01-qc_and_preprocessing.R)
2. [03-differential_abundance_milo.R](scripts/03-differential_abundance_milo.R)
3. [04-differential_abundance_model.R](scripts/04-differential_abundance_model.R)
4. [05-differential_gene_expression.R](scripts/05-differential_gene_expression.R)
5. [06-gene_set_enrichment_analysis.R](scripts/06-gene_set_enrichment_analysis.R)
6. [07-kmeans_clustering_degs.R](scripts/07-kmeans_clustering_degs.R)
7. [08-differential_signaling_domino.R](scripts/08-differential_signaling_domino.R)
8. [08.1-pySCENIC.sh](scripts/08.1-pySCENIC.sh)

## Plots
1. [Figure 2A](plots/fig2A-NhoodGraph_milo_kelly_45.pdf)
2. [Figure 2B](plots/fig2b-beeswarm_milo_kelly_45.pdf)
3. [Figure 2C](plots/fig2c-volcano_milo_kelly_45.pdf)
4. [Figure 2D-K](plots/fig2d-k-celltype_proportion_plots.pdf)
5. [Figure 3A](plots/fig3a-volcano_dge.pdf)
6. [Figure 3B](plots/fig3b-gsea_hallmark_all_celltypes.pdf)
7. [Figure 3C](plots/fig3c-kmeans10_centroids_heatmap_dge.pdf)

## Supplemental Data
Also available an Excel spreadsheet with each data set as a separate tab

1. [Differentially Expressed Genes](data/dge_pseudobulk_edgeR.txt)
   Differential gene expression analysis results with p-value < 0.01
2. [Gene Set Enrichment Analysis of DEGs](data/gsea_results_all_hallmark.csv)
   Gene Set Enrichment Analysis results with p-value < 0.05
3. [K-Means Clustering of DEGs](data/kmeans.k10.module_membership.txt)
   K means clustering analysis results for each gene, including the each cluster's and cell type's expression of the gene
4. [Differential Signaling Linkages](data/sig_diff_linkages_all.csv)
   Differential Signaling analysis results from DominoSignal with p < 0.05
5. [Incoming Ligand Senders](data/incoming_ligand_senders.csv)
    Top three cell types expressing each ligand in a receptor-ligand pair
