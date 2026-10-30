#!/bin/bash
#SBATCH --account=sament-lab
#SBATCH --job-name=domino
#SBATCH --output=domino_%j.log
#SBATCH --error=domino_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/domino

module load apptainer/1.3.6

CURDIR=$(pwd)

apptainer exec \
    --bind $CURDIR:$CURDIR \
    scenic/aertslab-pyscenic-0.12.1.sif \
    pyscenic grn \
    "$CURDIR/scenic/pbmc13k_counts.tsv" \
    "$CURDIR/scenic/allTFs_hg38.txt" \
    -o "$CURDIR/scenic/pbmc_adj_13k.tsv" \
    --num_workers 6 \
    --seed 123


apptainer exec \
    --bind $CURDIR:$CURDIR \
    scenic/aertslab-pyscenic-0.12.1.sif \
    pyscenic ctx \
    "$CURDIR/scenic/pbmc_adj_13k.tsv" \
    "$CURDIR/scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    "$CURDIR/scenic/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
    --annotations_fname "$CURDIR/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
    --expression_mtx_fname "$CURDIR/scenic/pbmc13k_counts.tsv" \
    --mode "dask_multiprocessing" \
    --output "$CURDIR/scenic/regulons_pbmc_13k.csv" \
    --num_workers 1


apptainer exec \
    --bind $CURDIR:$CURDIR \
    scenic/aertslab-pyscenic-0.12.1.sif \
    pyscenic aucell \
    "$CURDIR/scenic/pbmc60k_counts.tsv" \
    "$CURDIR/scenic/regulons_pbmc_13k.csv" \
    -o "$CURDIR/scenic/auc_pbmc_60k.csv"
