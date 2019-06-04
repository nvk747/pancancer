
#!/bin/bash

# Pipeline to reproduce pancancer pathways XYZ machine learning classifier
#
# Usage: bash XYZ_analysis.sh
#
# Output: Will train a pan cancer model to detect XYZ aberration. Will also
#         train a unique classifier within each specific cancer type

#--genes 'data/abc_genes.csv' \
        

# Set Constants
#abc_diseases='BLCA,BRCA,CESC,COAD,ESCA,GBM,HNSC,KICH,LGG,LIHC,LUAD,LUSC,'\
#'PAAD,PRAD,READ,SARC,SKCM,STAD,UCEC'
alphas='0.1,0.13,0.15,0.18,0.2,0.3,0.4,0.6,0.7'
l1_mixing='0.1,0.125,0.15,0.2,0.25,0.3,0.35'
abc_dir='classifiers/PI3K_TOTAL_MOD'

###############
# Step 1. Pan Cancer ABC classification
# abc_genes : Comma separated string of HUGO gene symbols file
# abc_diseases : Comma separated string of TCGA disease acronyms
###############
python scripts/pancancer_classifier.py \
        --genes 'data/abc_genes.csv' \
        --diseases 'data/abc_diseases.csv' \
        --drop \
        --copy_number \
        --remove_hyper \
        --alt_folder $abc_dir \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --keep_intermediate \
        --shuffled \
        --x_as_raw \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' \
        --alt_genes 'data/alt_abc_genes.csv' \
        --alt_diseases 'data/alt_abc_diseases.csv'     

###############
# Step 2. Within Disease type ABC classification
###############
python scripts/within_tissue_analysis.py \
        --genes 'data/abc_genes.csv' \
        --diseases 'data/abc_diseases.csv' \
        --remove_hyper \
        --alt_folder $abc_dir'/within_disease' \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --x_matrix 'data/pancan_rnaseq_freeze.tsv' \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' 

###############
# Step 3. Get scores for all samples and visualize distribution of scores
###############
python scripts/apply_weights.py \
        --x_matrix 'data/pancan_rnaseq_freeze.tsv' \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' \
        --classifier $abc_dir \
        --copy_number

python scripts/visualize_decisions.py \
        --scores $abc_dir 

python scripts/map_mutation_class.py \
        --scores $abc_dir \
        --path_genes 'data/pi3k_genes.csv'\
        --copy_number \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_raw_mut 'data/raw/mc3.v0.2.8.PUBLIC.maf'

python scripts/abc_alternative_genes_pathwaymapper.py \
        --alt_folder $abc_dir

python scripts/copy_burden_merge.py \
        --classifier_folder $abc_dir \
        --filename_burden 'data/seg_based_scores.tsv'

python scripts/abc_pathway_count_heatmaps.py \
        --genes 'data/abc_genes.csv' \
        --path_genes 'data/pi3k_genes.csv' \
        --alt_folder $abc_dir \
        --x_matrix 'data/pancan_rnaseq_freeze.tsv' \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv'

Rscript --vanilla scripts/viz/abc_summary_figures.R \
        --alt_folder $abc_dir