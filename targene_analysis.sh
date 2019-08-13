
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
#alphas='0.1,0.13,0.15,0.18,0.2,0.25,0.3'
#l1_mixing='0.15,0.155,0.16,0.2,0.25,0.3,0.4'
gene_dir='classifiers/PI3K_GAIN_TEST'
ex_vlog='/mnt/isilon/data/w_gmi/blanked2lab/vijay/HGDC_GSE69240/vlog_trans.csv' 
sign='/mnt/isilon/data/w_gmi/blanked2lab/vijay/HGDC_GSE69240/sign.csv'
ccle_rnaseq='/data/vijay/git/pancancer/data/ccle_rnaseq_genes_rpkm_20180929_mod.gct'
ccle_mut='/data/vijay/git/pancancer/data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct'
ccle_maf='/data/vijay/git/pancancer/data/CCLE_DepMap_18Q1_maf_20180207.txt'
phar_data='/data/vijay/git/pancancer/data/CCLE_NP24.2009_Drug_data_2015.02.24.csv'

###############
# Step 1. Pan Cancer alternate gene classification
###############
 
python scripts/pancancer_classifier.py \
        --genes 'NFE2L2' \
        --diseases 'BLCA,HNSC,LUSC,UCEC' \
        --drop \
        --copy_number \
        --remove_hyper \
        --alt_folder 'classifiers/NFE2L2_NEW' \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --keep_intermediate \
        --shuffled \
        --x_as_raw \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' \ 

###############
# Step 2. Pan Cancer target gene classification
############### 

python scripts/pancancer_classifier.py \
        --genes 'ERBB2,PIK3CA,KRAS,AKT1' \
        --diseases 'BLCA,BRCA,CESC,COAD,ESCA,LUAD,LUSC,OV,PRAD,READ,STAD,UCEC,UCS' \
        --drop \
        --copy_number \
        --remove_hyper \
        --alt_folder $gene_dir \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --keep_intermediate \
        --shuffled \
        --x_as_raw \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' \
        --alt_genes 'NFE2L2' \
        --alt_diseases 'BLCA,HNSC,LUSC,UCEC'   

###############
# Step 3. Within Disease type gene/alt_gene classification
###############

python scripts/within_tissue_analysis.py \
        --genes 'NFE2L2' \
        --diseases 'BLCA,HNSC,LUSC,UCEC' \
        --remove_hyper \
        --alt_folder 'classifiers/NFE2L2_NEW/within_disease' \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --x_as_raw \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv'

Rscript scripts/compare_within_models.R --pancan_summary 'classifiers/NFE2L2_NEW' \
        --within_dir 'classifiers/NFE2L2_NEW' \

python scripts/within_tissue_analysis.py \
        --genes 'ERBB2,PIK3CA,KRAS,AKT1' \
        --diseases 'BLCA,BRCA,CESC,COAD,ESCA,LUAD,LUSC,OV,PRAD,READ,STAD,UCEC,UCS' \
        --remove_hyper \
        --alt_folder $gene_dir'/within_disease' \
        --alphas $alphas \
        --l1_ratios $l1_mixing \
        --x_as_raw \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' 

Rscript scripts/compare_within_models.R --pancan_summary $gene_dir \
        --within_dir $gene_dir'/within_disease' \
        --alt_gene 'classifiers/NFE2L2_NEW'

###############
# Step 4. Get scores for all samples and visualize distribution of scores
###############
python scripts/apply_weights.py \
        --x_matrix 'data/pancan_rnaseq_freeze.tsv' \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' \
        --classifier $gene_dir \
        --copy_number

python scripts/visualize_decisions.py \
        --scores $gene_dir 

python scripts/apply_weights.py \
        --x_matrix 'data/pancan_rnaseq_freeze.tsv' \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv' \
        --classifier 'classifiers/NFE2L2_NEW' \
        --copy_number

python scripts/visualize_decisions.py \
        --scores 'classifiers/NFE2L2_NEW'

python scripts/map_mutation_class.py \
        --scores $gene_dir \
        --path_genes 'data/pi3k_genes.csv'\
        --copy_number \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_raw_mut 'data/raw/mc3.v0.2.8.PUBLIC.maf'

python scripts/targene_alternative_genes_pathwaymapper.py \
        --genes 'ERBB2,PIK3CA,KRAS,AKT1' \
        --scores $gene_dir \
        --path_genes 'data/pi3k_genes.csv'\
        --copy_number \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_sample 'data/sample_freeze.tsv' 

python scripts/copy_burden_merge.py \
        --classifier_folder $gene_dir \
        --filename_burden 'data/seg_based_scores.tsv'

python scripts/targene_pathway_count_heatmaps.py \
        --genes 'ERBB2,PIK3CA,KRAS,AKT1' \
        --path_genes 'data/pi3k_genes.csv' \
        --scores $gene_dir \
        --x_matrix 'data/pancan_rnaseq_freeze.tsv' \
        --filename_mut 'data/pancan_mutation_freeze.tsv' \
        --filename_mut_burden 'data/mutation_burden_freeze.tsv' \
        --filename_copy_loss 'data/copy_number_loss_status.tsv' \
        --filename_copy_gain 'data/copy_number_gain_status.tsv' \
        --filename_sample 'data/sample_freeze.tsv' \
        --filename_cancer_gene_classification 'data/vogelstein_cancergenes.tsv'

Rscript --vanilla scripts/viz/targene_summary_figures.R \
        --classifier_folder $gene_dir

python scripts/viz/targene_cell_line_predictions.py \
        --classifier $gene_dir \
        --targenes 'ERBB2_MUT,PIK3CA_MUT,KRAS_MUT,AKT1_MUT'  \
        --path_genes 'data/pi3k_genes.csv' \
        --ccle_rnaseq $ccle_rnaseq \
        --ccle_maf $ccle_maf \
        --ccle_mut $ccle_mut \
        --phar_data $phar_data 

# classifier prediction on external samples with HGDC_GSE69240 data 
python scripts/viz/external_sample_pred_targene_classifier.py \
        --classifier $gene_dir \
        --ex_vlog $ex_vlog \
        --sign $sign 

Rscript --vanilla scripts/viz/targene_ccle_pharmacology.R \
        --classifier $gene_dir \