#!/usr/bin/env python
# coding: utf-8

# # Cell Line Analysis
# 
# We sought to validate the Ras classifier trained on TCGA pan-cancer data by generating predictions on cell line data. A good classifier should generalize to predicting Ras status in other samples. We apply the classifier on two datasets:
# 
# 1. [Cancer Cell Line Encyclopedia (CCLE)](https://software.broadinstitute.org/software/cprg/?q=node/11) Gene Expression data.
#   * 1020 cell lines with matching gene expression and mutation calls
#   * Pharmacologic profiling of 24 drugs over 504 cell lines
# 2. GDSC data: These data were accessed via publicly available resources with help from links in the [UCSD-CCAL Onco-GPS github repository](https://github.com/UCSD-CCAL/onco-gps-paper-analysis)
#  data from : A landscape of pharmacogenomic interactions in cancer', Iorio F et al. Cell. 2016
#  (https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip) Gene Expression data.
#  390/1000 cell lines with matching gene expression and mutation calls from CCLE data
#  Pharmacologic profiling of ~500 drugs over 1000 cell lines from two different studies GDSC1 and GDSC2 datasets
#  GDSC1_data: ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC1_fitted_dose_response_15Oct19.xlsx
#  GDSC1_data: ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/GDSC2_fitted_dose_response_15Oct19.xlsx
#  we replaced all the GDSC celllines with CCLE cell line names and for convenient processing. 
import os
import numpy as np
import pandas as pd
from decimal import Decimal
from scipy.stats import ttest_ind
from statsmodels.stats.proportion import proportions_chisquare
from sklearn.preprocessing import StandardScaler
from Bio.SeqUtils import IUPACData
import matplotlib.pyplot as plt
import seaborn as sns
import plotnine as gg
import argparse

# Store protein change dictionary
aa = IUPACData.protein_letters_1to3_extended
#get_ipython().run_line_magic('matplotlib', 'inline')

parser = argparse.ArgumentParser()

parser.add_argument('-t', '--targenes', default= 'KRAS_MUT,NRAS_MUT,HRAS_MUT',
                    help='string of the genes to extract or gene list file')
parser.add_argument('-p', '--path_genes',
                    help='pathway gene list file')
parser.add_argument('-c', '--classifier', default= None,
                    help='location of classifier_summary file')
parser.add_argument('-r', '--ccle_rnaseq',default= None,
                    help='path for ccle_rnaseq data file')
parser.add_argument('-m', '--ccle_mut',
                    help='path for ccle mutational data file')
parser.add_argument('-a', '--ccle_maf',
                    help='path for ccle variant data file')
parser.add_argument('-n', '--gdsc_rnaseq',default= None,
                    help='path for gdsc_rnaseq data file')
parser.add_argument('-u', '--gdsc_mut',
                    help='path for gdsc/ccle common mutational data file')
parser.add_argument('-e', '--gdsc1_phar',
                    help='path for GDSC1 pharmacological data file')
parser.add_argument('-f', '--gdsc2_phar',
                    help='path for GDSC2 pharmacological data file')
args = parser.parse_args()

# Load PI3K_gain Classifier Coefficients
# classifier_file = os.path.join('..', 'classifiers', 'ERBB2_PIK3CA_KRAS_AKT1', 'classifier_summary.txt')
# with open(classifier_file) as class_fh:
#    for line in class_fh:
#        line = line.strip().split('\t')
#        if line[0] == 'Coefficients:':
#            all_coef_df = pd.read_table(os.path.join('..', line[1]), index_col=0)

# Only non-zero coefficients contribute to model performance

classifier = args.classifier
classifier_file = os.path.join( classifier , "classifier_summary.txt")
all_coef_df = pd.read_table(os.path.join( classifier , "classifier_coefficients.tsv"), index_col=0)
coef_df = all_coef_df[all_coef_df['abs'] > 0]
coef_df.head(10)

# ## Part 1: CCLE
# 
# Note - This data was also retrieved from the Onco-GPS paper analysis repository

#ccle_file_name = os.path.join('..', '..', 'onco-gps-paper-analysis', 'data',
#                              'rpkm__gene_x_ccle_cellline.gct')

ccle_file_name = args.ccle_rnaseq or os.path.join('..','data','ccle_rnaseq_genes_rpkm_20180929_mod.gct')
ccle_df = pd.read_csv(ccle_file_name, skiprows=2, index_col=0)
#ccle_df = ccle_df.drop_duplicates(subset='Description',keep = 'first')
ccle_df = ccle_df[~ccle_df.index.duplicated()]

# Subset to common genes in the classifier and CCLE data
common_genes = list(set(coef_df['feature']) & set(ccle_df.index))
common_ccle_coef = coef_df[coef_df['feature'].isin(common_genes)]

ccle_df = ccle_df.loc[common_ccle_coef['feature'], ccle_df.columns[1:]]

scaled_fit = StandardScaler().fit(ccle_df.T)
ccle_df = pd.DataFrame(scaled_fit.transform(ccle_df.T),
                            index=ccle_df.columns,
                            columns=ccle_df.index)

ccle_df = ccle_df.T

# Get the weights ready for applying the classifier
apply_weights = pd.DataFrame(common_ccle_coef['weight'])
apply_weights.index = common_ccle_coef.feature

# Apply a logit transform [y = 1/(1+e^(-wX))] to output probabilities
result_ccle = apply_weights.T.dot(ccle_df)
result_ccle = 1 / (1 + np.exp(-1 * result_ccle))

# Distribution of predictions of the Ras Classifier applied to CCLE data
result_ccle.T.hist();
r = os.path.join(classifier,'figures','ccle_histogram.png')
plt.savefig(r)
plt.close()

# Load CCLE Mutation Data
#ccle_mut_file_name = os.path.join('..', '..', 'onco-gps-paper-analysis', 'data', 
#                                  'mutation__gene_x_ccle_cellline.gct')
ccle_mut_file_name = args.ccle_mut or os.path.join('..','data','CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct')
ccle_all_mut_df = pd.read_table(ccle_mut_file_name, skiprows=2, index_col=0)

# Load CCLE Variant Data
#ccle_maf_file = 'https://data.broadinstitute.org/ccle/CCLE_DepMap_18Q1_maf_20180207.txt'
ccle_maf_file = args.ccle_maf or os.path.join('..','data','CCLE_DepMap_18Q1_maf_20180207.txt')
ccle_maf_df = pd.read_table(ccle_maf_file, index_col=15)

# Identify all cell lines with mutations in targene genes, also subset braf mutant samples

targenes = args.targenes.split(',')
targene_status = ccle_all_mut_df.loc[targenes, :].T.apply(max, axis=1)

# BRAF mutations do not contribute to targene status in this case
ccle_mut_df = (
    ccle_all_mut_df.loc[targenes, :].T
    .assign(targene_status=targene_status).drop(['Description'])
    )

# Join classifier scores with mutation status
ccle_full_df = ccle_mut_df.join(result_ccle.T).dropna()
ccle_full_df = ccle_full_df.assign(sample_name = ccle_full_df.index)
ccle_full_df = ccle_full_df.sort_values(by='weight', ascending=False)
ccle_full_df.index.name = 'cell_line'

# Write CCLE Scores to file
results_folder = os.path.join(classifier, 'results')
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
ccle_scores_file = os.path.join(classifier, 'results', 'ccle_targene_classifier_scores.tsv')
ccle_full_df.to_csv(ccle_scores_file, sep='\t')

# ### Perform a t-test on classifier weights across groups
# targene mutant vs. targene wildtype 
targene_mutant = ccle_full_df[ccle_full_df['targene_status'] == 1]
targene_wt = ccle_full_df[ccle_full_df['targene_status'] == 0]

# Output t-test results
t_results_ccle_targene = ttest_ind(a = targene_mutant['weight'],
                          b = targene_wt['weight'], equal_var = False)
print('targene Status:')
print(t_results_ccle_targene)

# Use Seaborn for the 2nd plot
import seaborn as sns
import matplotlib.pyplot as pl

sns.set_style("whitegrid")
sns.set_context("paper", rc={"font.size":11, "axes.titlesize":11, "axes.labelsize":16,
                             'xtick.labelsize':11, 'ytick.labelsize':11,  'figure.facecolor': 'white'})


cell_line_folder = os.path.join(classifier, 'figures', 'cell_line')
if not os.path.exists(cell_line_folder):
    os.makedirs(cell_line_folder)

# Plot Results for targene alone
x1, x2 = 0, 1
y1, y2,h = 1.05, 1.0, 0.03

plt.rcParams['figure.figsize']=(3.5, 4)
ax2 = sns.boxplot(x="targene_status", y="weight", data=ccle_full_df,
                 palette = {0: "lightgreen",1: 'yellow'},
                 fliersize=0)
ay2 = sns.stripplot(x='targene_status', y='weight', data=ccle_full_df, 
                   dodge=False,
                   palette = {0: "blue", 1: 'red'},
                   jitter=0.12, size=2, alpha=0.65)

ax2.axes.set_ylim(0, 1.2)
ax2.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1,''])
ax2.set_xticklabels(['Targene_WT', 'Targene_Mut'])
ax2.set_ylabel('Targene Classifier Score')
ax2.set_xlabel('CCLE Data')
ax2.legend
pl.axhline(0.5, color='black', linestyle='dashed', linewidth=1)
ay2.axes.set_ylim(0, 1.2)

# Add targene T-Test Results
pl.plot([x1, x1, x2, x2], [y1, y1+h, y1+h, y1], lw=1.2, c='black')
pl.text(.6, y1+h, "{:.2E}".format(Decimal(t_results_ccle_targene.pvalue)),
         ha='center', va='bottom', color="black")

pl.tight_layout()
ccle_fig_file = os.path.join(classifier, 'figures', 'cell_line', 'ccle_targene_WT_MUT_predictions.pdf')
pl.savefig(ccle_fig_file)
plt.close()

# ### What percentage of correct classifications in CCLE data?

# Assign a label to what the predictions are given classifier scores
ccle_full_df = ccle_full_df.assign(predictions = 'wild-type')
ccle_full_df.loc[ccle_full_df['weight'] > 0.5, 'predictions'] = 'mutant'

# Stratify cell lines based on predictions and ground truth status
positive_targene_predictions_ccle = ccle_full_df[ccle_full_df['weight'] > 0.5]
negative_targene_predictions_ccle = ccle_full_df[ccle_full_df['weight'] <= 0.5]

positive_targene_lines_ccle = ccle_full_df[ccle_full_df['targene_status'] == 1]
negative_targene_lines_ccle = ccle_full_df[ccle_full_df['targene_status'] == 0]

# Of wild-type targene cell lines, how many are predicted correctly?
# True Negative Rate, Specificity
negative_targene_lines_ccle['predictions'].value_counts()

# Of mutated targene cell lines, how many are predicted correctly?
# True Positive Rate (TPR), Recall, Sensitivity
positive_targene_lines_ccle['predictions'].value_counts()

# Of the wild-type predictions, how many are actually wild-type?
# Negative Predictive Value (NPV)
neg_ccle_results = negative_targene_predictions_ccle['targene_status'].value_counts()
true_neg = neg_ccle_results[0]
predicted_condition_neg = neg_ccle_results.sum()

print('{} out of {} Targene wild-type predictions '
      'are true ({:.1f}%)'.format(true_neg, predicted_condition_neg,
                                  true_neg * 100 / predicted_condition_neg))

# Of the mutated predictions, how many are actually mutated?
# Positive Predictive Value (PPV) -or- precision
pos_ccle_results = positive_targene_predictions_ccle['targene_status'].value_counts()
false_pos, true_pos = pos_ccle_results
predicted_condition_pos = pos_ccle_results.sum()

print('{} out of {} Targene mutation predictions '
      'are true ({:.1f}%)'.format(true_pos, predicted_condition_pos,
                                  true_pos * 100 / predicted_condition_pos))

total_correct = true_pos + true_neg
print('{} of {} Total cell lines '
      'predicted correctly ({:.1f}%)'.format(total_correct, ccle_full_df.shape[0],
                                             total_correct * 100 / ccle_full_df.shape[0]))

# ### Add CCLE Variant Scores (nucleotide and amino acid) to Supplementary Data Files

# Load TCGA PanCanAtlas Core targene Pathway genes
path_genes_file = args.path_genes or os.path.join('..', 'data', 'pi3k_genes.csv')
path_core_df = pd.read_table(path_genes_file)

# Subset MAF file to targene pathway variants and merge with CCLE classifier scores
path_pathway_genes = path_core_df['genes'].tolist()
all_common_lines = set(ccle_maf_df.index).intersection(set(ccle_full_df.index))

# Subset to common cell lines
subset_maf = ccle_maf_df.loc[list(all_common_lines), :]
subset_maf = (
    subset_maf.query('Hugo_Symbol in @path_pathway_genes')
    .loc[:, ['Hugo_Symbol', 'Protein_Change', 'cDNA_Change']]
    .merge(ccle_full_df, left_index=True, right_index=True)
)

subset_maf.head(3)

# Get the mean classifier scores for CCLE nucleotide variants
mean_nuc_data = (
    pd.DataFrame(subset_maf
                 .groupby(['cDNA_Change', 'Hugo_Symbol'])['weight']
                 .mean())
)
mean_nuc_data.columns = ['ccle_mean_weight']
mean_nuc_data = mean_nuc_data.reset_index()

# Get the sd classifier scores for CCLE variants
sd_nuc_data = (
    pd.DataFrame(subset_maf
                 .groupby(['cDNA_Change', 'Hugo_Symbol'])['weight']
                 .std())
)
sd_nuc_data.columns = ['ccle_sd_weight']
sd_nuc_data = sd_nuc_data.reset_index()

# Counts of CCLE variants altering amino acids
count_nuc_data = (
    pd.DataFrame(subset_maf
                 .groupby(['cDNA_Change', 'Hugo_Symbol'])['weight']
                 .count())
)
count_nuc_data.columns = ['ccle_count']
count_nuc_data = count_nuc_data.reset_index()

# Merge protein data
nuc_merge_on = ['Hugo_Symbol', 'cDNA_Change']
nuc_change_df = (
    mean_nuc_data.merge(sd_nuc_data,
                        left_on=nuc_merge_on, right_on=nuc_merge_on)
    .merge(count_nuc_data, left_on=nuc_merge_on, right_on=nuc_merge_on)
)

nuc_change_df.sort_values('ccle_count').tail(5)

#data_s4_file = os.path.join('..', 'classifiers', 'ERBB2_PIK3CA_KRAS_AKT1', 'tables',
#                            'nucleotide_mutation_scores.tsv')
data_s4_file = os.path.join(classifier, 'tables','nucleotide_mutation_scores.tsv')
data_s4_df = pd.read_table(data_s4_file)

# Merge the CCLE nucleotide scores
data_s4_df = data_s4_df.merge(nuc_change_df, left_on = ['Hugo_Symbol', 'HGVSc'],
                                    right_on = ['Hugo_Symbol', 'cDNA_Change'],
                              how='outer')
updated_data_s4_df = data_s4_df.sort_values(by='count', ascending=False)

#updated_data_s4_file = os.path.join('..', 'classifiers', 'ERBB2_PIK3CA_KRAS_AKT1', 'tables', 'updated_Data_S4.csv')
updated_data_s4_file = os.path.join(classifier,'tables','updated_Data_S4.csv')
updated_data_s4_df.to_csv(updated_data_s4_file, sep=',', index=False)

# Get the mean classifier scores for CCLE variants
mean_protein_data = (
    pd.DataFrame(subset_maf
                 .groupby(['Protein_Change', 'Hugo_Symbol'])['weight']
                 .mean())
)
mean_protein_data.columns = ['ccle_mean_weight']
mean_protein_data = mean_protein_data.reset_index()

# Get the sd classifier scores for CCLE variants
sd_protein_data = (
    pd.DataFrame(subset_maf
                 .groupby(['Protein_Change', 'Hugo_Symbol'])['weight']
                 .std())
)
sd_protein_data.columns = ['ccle_sd_weight']
sd_protein_data = sd_protein_data.reset_index()

# Counts of CCLE variants altering amino acids
count_protein_data = (
    pd.DataFrame(subset_maf
                 .groupby(['Protein_Change', 'Hugo_Symbol'])['weight']
                 .count())
)
count_protein_data.columns = ['ccle_count']
count_protein_data = count_protein_data.reset_index()

# Merge protein data
merge_on = ['Hugo_Symbol', 'Protein_Change']
protein_change_df = (
    mean_protein_data.merge(sd_protein_data,
                            left_on=merge_on, right_on=merge_on)
    .merge(count_protein_data, left_on=merge_on, right_on=merge_on)
)
protein_change_df.sort_values('ccle_count').tail(5)

# Convert amino acid to 3 letters
protein_convert = [''.join([aa[x] if x in aa.keys() else x for x in y]) 
                   for y in protein_change_df['Protein_Change']]
protein_change_df = protein_change_df.assign(conversion = protein_convert)

#data_s5_file = os.path.join('..', 'classifiers', 'ERBB2_PIK3CA_KRAS_AKT1', 'tables',
#                            'amino_acid_mutation_scores.tsv')

data_s5_file = os.path.join(classifier, 'tables',
                            'amino_acid_mutation_scores.tsv')
data_s5_df = pd.read_table(data_s5_file)

# Merge the CCLE protein scores
data_s5_df = data_s5_df.merge(protein_change_df, left_on = ['Hugo_Symbol', 'HGVSp'],
                                    right_on = ['Hugo_Symbol', 'conversion'],
                              how='outer')

# Sort by the total number of mutations observed
updated_data_s5_df = (
    data_s5_df.drop(['Protein_Change'], axis=1).sort_values(by='count', ascending=False)
)

updated_data_s5_file = os.path.join(classifier,'tables', 'updated_Data_S5.csv')
updated_data_s5_df.to_csv(updated_data_s5_file, sep=',', index=False)

# GDSC cellline expression and pharmacological evaluation:

# GDSC 1000 cell line expression is downloaded from 
# https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/
# /Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip and RNA expression for 17737 genes from 382 cellines
# among CCLE and GDSC was used. All GDSC cellines names are replaced by CCLE celllines

gdsc_file_name = args.gdsc_rnaseq or os.path.join('..','data','GDSC_EXP_CCLE_converted_name.csv')
gdsc_df = pd.read_csv(gdsc_file_name,sep=',', index_col=0)

# Subset to common genes in the classifier and gdsc data
common_genes = list(set(coef_df['feature']) & set(gdsc_df.index))
common_gdsc_coef = coef_df[coef_df['feature'].isin(common_genes)]
gdsc_df = gdsc_df.loc[common_gdsc_coef['feature'], gdsc_df.columns[1:]]
scaled_fit = StandardScaler().fit(gdsc_df.T)
gdsc_df = pd.DataFrame(scaled_fit.transform(gdsc_df.T),
                            index=gdsc_df.columns,
                            columns=gdsc_df.index)

# Get the weights ready for applying the classifier
apply_weights = pd.DataFrame(common_gdsc_coef['weight'])
apply_weights.index = common_gdsc_coef.feature

# Apply a logit transform [y = 1/(1+e^(-wX))] to output probabilities
result_gdsc = apply_weights.T.dot(gdsc_df.T)
result_gdsc = 1 / (1 + np.exp(-1 * result_gdsc))

# Distribution of predictions of the targene Classifier applied to GDSC cell-line data
result_gdsc.T.hist();
r = os.path.join(classifier,'figures','gdsc_scores_histogram.png')
plt.savefig(r)
plt.close()

gdsc_mut_file_name = args.gdsc_mut or os.path.join('..','data', 'GDSC_CCLE_common_mut_cnv_binary.csv')
gdsc_all_mut_df = pd.read_csv(gdsc_mut_file_name, index_col=0)

# Identify all cell lines with mutations in Ras genes, also subset BRAF mutant samples
targene_status = gdsc_all_mut_df.loc[targenes, :].T.apply(max, axis=1)
gdsc_mut_df = (
    gdsc_all_mut_df.loc[targenes, :].T
    .assign(targene_status=targene_status)
    )
# Join classifier scores with mutation status
gdsc_full_df = gdsc_mut_df.join(result_gdsc.T).dropna()
gdsc_full_df = gdsc_full_df.assign(sample_name = gdsc_full_df.index)
gdsc_full_df = gdsc_full_df.sort_values(by='weight', ascending=False)
gdsc_full_df.index.name = 'cell_line'

# Write gdsc Scores to file
gdsc_scores_file = os.path.join(classifier, 'results', 'gdsc_targene_classifier_scores.tsv')
gdsc_full_df.to_csv(gdsc_scores_file, sep='\t')

# ### Perform a t-test on classifier weights across groups
# targene mutant vs. targene wildtype 
targene_mutant = gdsc_full_df[gdsc_full_df['targene_status'] == 1]
targene_wt = gdsc_full_df[gdsc_full_df['targene_status'] == 0]

# Output t-test results
t_results_gdsc_targene = ttest_ind(a = targene_mutant['weight'],
                          b = targene_wt['weight'], equal_var = False)
print('targene Status:')
print(t_results_gdsc_targene)
# Plotting GDSC targene predictions
import seaborn as sns
import matplotlib.pyplot as pl

sns.set(style="whitegrid")
sns.set_context("paper", rc={"font.size":11, "axes.titlesize":11, "axes.labelsize":16,
                             'xtick.labelsize':11, 'ytick.labelsize':11})
# Plot Results for targene alone
x1, x2 = 0, 1
y1, y2,h = 1.05, 1.0, 0.03

plt.rcParams['figure.figsize']=(3.5, 4)
ax2 = sns.boxplot(x="targene_status", y="weight", data=gdsc_full_df,
                 palette = {0: "lightgreen",1: 'yellow'},
                 fliersize=0)
ay2 = sns.stripplot(x='targene_status', y='weight', data=gdsc_full_df, 
                   dodge=False,
                   palette = {0: "blue", 1: 'red'},
                   jitter=0.12, size=2, alpha=0.65)

ax2.axes.set_ylim(0, 1.2)
ax2.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1,''])
ax2.set_xticklabels(['Targene_WT', 'Targene_Mut'])
ax2.set_ylabel('Targene Classifier Score')
ax2.set_xlabel('GDSC Data')
ax2.legend
pl.axhline(0.5, color='black', linestyle='dashed', linewidth=1)
ay2.axes.set_ylim(0, 1.2)

# Add targene T-Test Results
pl.plot([x1, x1, x2, x2], [y1, y1+h, y1+h, y1], lw=1.2, c='black')
pl.text(.6, y1+h, "{:.2E}".format(Decimal(t_results_gdsc_targene.pvalue)),
         ha='center', va='bottom', color="black")
pl.tight_layout()
gdsc_fig_file = os.path.join(classifier, 'figures', 'cell_line', 'gdsc_targene_WT_MUT_predictions.pdf')
pl.savefig(gdsc_fig_file)
plt.close()

# What percentage of correct classifications in GDSC data?
# Assign a label to what the predictions are given classifier scores
gdsc_full_df = gdsc_full_df.assign(predictions = 'wild-type')
gdsc_full_df.loc[gdsc_full_df['weight'] > 0.5, 'predictions'] = 'mutant'

# Stratify cell lines based on predictions and ground truth status
positive_targene_predictions_gdsc = gdsc_full_df[gdsc_full_df['weight'] > 0.5]
negative_targene_predictions_gdsc = gdsc_full_df[gdsc_full_df['weight'] <= 0.5]

positive_targene_lines_gdsc = gdsc_full_df[gdsc_full_df['targene_status'] == 1]
negative_targene_lines_gdsc = gdsc_full_df[gdsc_full_df['targene_status'] == 0]

# Of wild-type targene cell lines, how many are predicted correctly?
# True Negative Rate, Specificity
negative_targene_lines_gdsc['predictions'].value_counts()

# Of mutated targene cell lines, how many are predicted correctly?
# True Positive Rate (TPR), Recall, Sensitivity
positive_targene_lines_gdsc['predictions'].value_counts()

# Of the wild-type predictions, how many are actually wild-type?
# Negative Predictive Value (NPV)
neg_gdsc_results = negative_targene_predictions_gdsc['targene_status'].value_counts()
true_neg = neg_gdsc_results[0]
predicted_condition_neg = neg_gdsc_results.sum()

print('{} out of {} targene wild-type predictions '
      'are true ({:.1f}%)'.format(true_neg, predicted_condition_neg,
                                  true_neg * 100 / predicted_condition_neg))

# Of the mutated predictions, how many are actually mutated?
# Positive Predictive Value (PPV) -or- precision
pos_gdsc_results = positive_targene_predictions_gdsc['targene_status'].value_counts()
false_pos, true_pos = pos_gdsc_results
predicted_condition_pos = pos_gdsc_results.sum()

print('{} out of {} targene mutation predictions '
      'are true ({:.1f}%)'.format(true_pos, predicted_condition_pos,
                                  true_pos * 100 / predicted_condition_pos))

total_correct = true_pos + true_neg
print('{} of {} Total cell lines '
      'predicted correctly ({:.1f}%)'.format(total_correct, gdsc_full_df.shape[0],
                                             total_correct * 100 / gdsc_full_df.shape[0]))

# GDSC classifier scores for GDSC1_cell_line_phramacological_evaluation
# Load in gdsc1 pharmacological results
pharm1_file = args.gdsc1_phar or os.path.join('..', 'data', 'gdsc1_ccle_pharm_fitted_dose_data.txt')
pharm_df = pd.read_table(pharm1_file, index_col=0)
pharm_df = pharm_df.assign(tissue = [' '.join(x[1:]) for x in pharm_df.index.str.split('_')])

pharm_full_df = pharm_df.merge(gdsc_full_df, left_index=True, right_index=True)

common_celllines_pharm = set(gdsc_full_df.index).intersection(set(pharm_df.index))
print('There are {} cell lines in common'.format(len(common_celllines_pharm)))

pharm_full_df['Compound'].value_counts()
pharm_full_df['tissue'].value_counts()

# What is the cell line tissue representation?
compound_heatmap = pd.pivot_table(pharm_full_df[['tissue', 'Compound']],
                                  columns='tissue', index='Compound',
                                  aggfunc=len)

compound_heatmap = pd.DataFrame(compound_heatmap.unstack()).reset_index()
compound_heatmap.columns = ['tissue', 'Compound', 'count']
compound_heatmap = compound_heatmap.sort_values(by=['tissue', 'Compound'])

# save compound_heatmap as table
gdsc1_comp_heatmap = os.path.join(classifier,'tables', 'gdsc1_compound_heatmap.csv')
compound_heatmap.to_csv(gdsc1_comp_heatmap, sep=',', index=False)

pharm_file = os.path.join(classifier, 'tables', 'gdsc1_targene_pharmacology_predictions.tsv')
pharm_full_df.to_csv(pharm_file, sep='\t')

# GDSC classifier scores for GDSC2_cell_line_phramacological_evaluation
# Load in gdsc2 pharmacological results
pharm2_file = args.gdsc2_phar or os.path.join('..', 'data', 'gdsc2_ccle_pharm_fitted_dose_data.txt')
pharm_df = pd.read_table(pharm1_file, index_col=0)
pharm_df = pharm_df.assign(tissue = [' '.join(x[1:]) for x in pharm_df.index.str.split('_')])

pharm_full_df = pharm_df.merge(gdsc_full_df, left_index=True, right_index=True)

common_celllines_pharm = set(gdsc_full_df.index).intersection(set(pharm_df.index))
print('There are {} cell lines in common'.format(len(common_celllines_pharm)))

pharm_full_df['Compound'].value_counts()
pharm_full_df['tissue'].value_counts()

# What is the cell line tissue representation?
compound_heatmap = pd.pivot_table(pharm_full_df[['tissue', 'Compound']],
                                  columns='tissue', index='Compound',
                                  aggfunc=len)

compound_heatmap = pd.DataFrame(compound_heatmap.unstack()).reset_index()
compound_heatmap.columns = ['tissue', 'Compound', 'count']
compound_heatmap = compound_heatmap.sort_values(by=['tissue', 'Compound'])

# save compound_heatmap as table
gdsc2_comp_heatmap = os.path.join(classifier,'tables', 'gdsc2_compound_heatmap.csv')
compound_heatmap.to_csv(gdsc2_comp_heatmap, sep=',', index=False)

pharm_file = os.path.join(classifier, 'tables', 'gdsc2_targene_pharmacology_predictions.tsv')
pharm_full_df.to_csv(pharm_file, sep='\t')

# CCLE classifier scores for GDSC1_cell_line_phramacological_evaluation
# Load in pharmacological results
pharm1_file = args.gdsc1_phar or os.path.join('..', 'data', 'gdsc1_ccle_pharm_fitted_dose_data.txt')
pharm_df = pd.read_table(pharm1_file, index_col=0)
pharm_df = pharm_df.assign(tissue = [' '.join(x[1:]) for x in pharm_df.index.str.split('_')])

pharm_full_df = pharm_df.merge(ccle_full_df, left_index=True, right_index=True)

common_celllines_pharm = set(ccle_full_df.index).intersection(set(pharm_df.index))
print('There are {} cell lines in common'.format(len(common_celllines_pharm)))

pharm_full_df['Compound'].value_counts()
pharm_full_df['tissue'].value_counts()

# What is the cell line tissue representation?
compound_heatmap = pd.pivot_table(pharm_full_df[['tissue', 'Compound']],
                                  columns='tissue', index='Compound',
                                  aggfunc=len)

compound_heatmap = pd.DataFrame(compound_heatmap.unstack()).reset_index()
compound_heatmap.columns = ['tissue', 'Compound', 'count']
compound_heatmap = compound_heatmap.sort_values(by=['tissue', 'Compound'])

pharm_file = os.path.join(classifier, 'tables', 'gdsc1_ccle_targene_pharmacology_predictions.tsv')
pharm_full_df.to_csv(pharm_file, sep='\t')

# CCLE classifier scores for GDSC2_cell_line_phramacological_evaluation
# Load in pharmacological results
pharm2_file = args.gdsc2_phar or os.path.join('..', 'data', 'gdsc2_ccle_pharm_fitted_dose_data.txt')
pharm_df = pd.read_table(pharm1_file, index_col=0)
pharm_df = pharm_df.assign(tissue = [' '.join(x[1:]) for x in pharm_df.index.str.split('_')])

pharm_full_df = pharm_df.merge(ccle_full_df, left_index=True, right_index=True)

common_celllines_pharm = set(ccle_full_df.index).intersection(set(pharm_df.index))
print('There are {} cell lines in common'.format(len(common_celllines_pharm)))

pharm_full_df['Compound'].value_counts()
pharm_full_df['tissue'].value_counts()

# What is the cell line tissue representation?
compound_heatmap = pd.pivot_table(pharm_full_df[['tissue', 'Compound']],
                                  columns='tissue', index='Compound',
                                  aggfunc=len)

compound_heatmap = pd.DataFrame(compound_heatmap.unstack()).reset_index()
compound_heatmap.columns = ['tissue', 'Compound', 'count']
compound_heatmap = compound_heatmap.sort_values(by=['tissue', 'Compound'])

pharm_file = os.path.join(classifier, 'tables', 'gdsc2_ccle_targene_pharmacology_predictions.tsv')
pharm_full_df.to_csv(pharm_file, sep='\t')
