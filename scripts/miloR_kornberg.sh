#!/bin/bash
#SBATCH --account=sament-lab
#SBATCH --job-name=miloR_korn
#SBATCH --output=logs/kornberg_%j.log
#SBATCH --error=output/kornberg_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --mail-type=END,FAIL

source ~/.bashrc
source /usr/local/packages/miniconda3/etc/profile.d/conda.sh
conda activate /autofs/projects-t3/idea/rowan.ohara/software/conda_pkg/rnaseq
Rscript /local/projects/idea/rohara/immunophenotype_SCZ/miloR_kornberg.R
