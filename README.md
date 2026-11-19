# Molecular adaptations of peripheral immune cells in an inflammation-associated schizophrenia subgroup
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
   - Graph-based representation of neighborhoods identified by Milo. Each node is a circle and depicts a neighborhood that is colored by their log fold change (if significant at FDR 30%***) comparing AGA+ SRD vs. AGA- participants (AGA- SRD plus AGA- HC). The size of each node depicts the neighborhood size (how many cells fall into that neighborhood) and edges depict the number of cells shared between adjacent neighborhoods. Non-differential abundance neighborhoods are colored white (FDR 30%). The layout of nodes corresponds to the position of the neighborhood index single cell in its respective UMAP embedding.
2. [Figure 2B](plots/fig2b-beeswarm_milo_kelly_45.pdf)
   - Beeswarm plot showing the distribution of log fold change in abundance between conditions in neighborhoods from different cell type clusters, with the most positive log fold change at the top of the plot and the most negative log fold change at the bottom of the plot. Positive log fold change indicates higher abundance in AGA+ SRD compared to AGA- participants, negative log fold change indicates lower abundance in AGA+ SRD compared to AGA- participants. Differential abundance neighborhoods are colored at FDR 30%.
3. [Figure 2C](plots/fig2c-volcano_milo_kelly_45.pdf)
   - Volcano plot showing each neighborhood’s respective log fold change and -log10(FDR) between AGA+ SRD compared to AGA- participants, with the color of each dot corresponding to a different immune cell subset. Solid line indicates FDR 15% and dashed line indicates FDR 30%.
4. [Figure 2D-K](plots/fig2d-k-celltype_proportion_plots.pdf)
   - Bar plots depicting pairwise comparisons using a binomial generalized linear model with a logit link, controlling for age and sex, of the cell type proportion (%) of each immune cell subset in AGA+ SRD (red), AGA- SRD (yellow), and HC (blue)
5. [Figure 3A](plots/fig3a-volcano_dge.pdf)
   - Volcano plot depicting the differential gene expression log fold change between AGA+ SRD and AGA- participants (AGA- SRD plus AGA- HC). Differentially expressed genes (DEGs) that are colored all have a log fold change magnitude of at least 1 and were nominally significant (p-value < 0.01). DEGs that passed multiple testing corrections are labeled.
6. [Figure 3B](plots/fig3b-gsea_hallmark_all_celltypes.pdf)
   - Dot plot of normalized enrichment score (NES) from gene set enrichment analysis (GSEA) using DEGs in AGA+ SRD compared to AGA- participants. The size of each dot corresponds to the -log10(FDR-adjusted p-value) of each Molecular Signature Database (MSigDB) Hallmark pathway for select immune cell types. The color of each dot corresponds to the NES of the pathway for that immune cell type, with red indicating positive enrichment in AGA+ SRD, and dark blue indicating negative enrichment in AGA+ SRD.
7. [Figure 3C](plots/fig3c-kmeans10_centroids_heatmap_dge.pdf)
   - Heatmap of K-means clustering of DEGs into immune transcriptional programs across immune cell types. The rows are annotated with the most overrepresented Gene Ontology biological process pathway for each module. The columns are annotated with immune cell types. Red  colors represent modules upregulated in that cell type in AGA+ SRD compared to AGA- participants. 

## Supplemental Data
Also available an Excel spreadsheet with each data set as a separate tab

1. [Differentially Expressed Genes](data/dge_pseudobulk_edgeR.txt)
   - Differential gene expression analysis results with p-value < 0.01
2. [Gene Set Enrichment Analysis of DEGs](data/gsea_results_all_hallmark.csv)
   - Gene Set Enrichment Analysis results with p-value < 0.05
3. [K-Means Clustering of DEGs](data/kmeans.k10.module_membership.txt)
   - K means clustering analysis results for each gene, including the each cluster's and cell type's expression of the gene
4. [Differential Signaling Linkages](data/sig_diff_linkages_all.csv)
   - Differential Signaling analysis results from DominoSignal with p < 0.05
5. [Incoming Ligand Senders](data/incoming_ligand_senders.csv)
   - Top three cell types expressing each ligand in a receptor-ligand pair
