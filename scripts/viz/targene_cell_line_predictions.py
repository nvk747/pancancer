#!/usr/bin/env python
# coding: utf-8

# # Cell Line Analysis
# 
# We sought to validate the Ras classifier trained on TCGA pan-cancer data by generating predictions on cell line data. A good classifier should generalize to predicting Ras status in other samples. We apply the classifier on two datasets:
# 
# 1. [GSE94937](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94937) from [Kim et al. 2017](http://doi.org/10.1016/j.cels.2017.08.002)
#   * Illumina NextSeq 5000 of Human Small Airway Epithelial Cells expressing KRAS G12V and wild-type KRAS
# 2. [Cancer Cell Line Encyclopedia (CCLE)](https://software.broadinstitute.org/software/cprg/?q=node/11) Gene Expression data.
#   * 737 cell lines with matching gene expression and mutation calls
#   * Pharmacologic profiling of 24 drugs over 504 cell lines
#   
# These data were accessed via publicly available resources with help from links in the [UCSD-CCAL Onco-GPS github repository](https://github.com/UCSD-CCAL/onco-gps-paper-analysis)

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
parser.add_argument('-d', '--ccle_rnaseq',default= None,
                    help='path for ccle_rnaseq data file')
parser.add_argument('-m', '--ccle_mut',
                    help='path for ccle mutational data file')
parser.add_argument('-a', '--ccle_maf',
                    help='path for ccle variant data file')
parser.add_argument('-r', '--phar_data',
                    help='path for ccle pharmacological data file')
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

# ## Part 2: CCLE
# 
# Note - This data was also retrieved from the Onco-GPS paper analysis repository

#ccle_file_name = os.path.join('..', '..', 'onco-gps-paper-analysis', 'data',
#                              'rpkm__gene_x_ccle_cellline.gct')

ccle_file_name = args.ccle_rnaseq or os.path.join('..','data','ccle_rnaseq_genes_rpkm_20180929_mod.gct')
ccle_df = pd.read_csv(ccle_file_name, skiprows=2, index_col=0)
ccle_df = ccle_df.drop_duplicates(subset='Description',keep = 'first')

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

# BRAF mutations do not contribute to Ras status in this case
ccle_mut_df = (
    ccle_all_mut_df.loc[targenes + ['BRAF_MUT'], :].T
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
ccle_scores_file = os.path.join(classifier, 'results', 'ccle_targene_BRAF_classifier_scores.tsv')
ccle_full_df.to_csv(ccle_scores_file, sep='\t')

# Use Seaborn for the 2nd plot
sns.set_style("whitegrid")
sns.set_context("paper", rc={"font.size":11, "axes.titlesize":11, "axes.labelsize":16,
                             'xtick.labelsize':11, 'ytick.labelsize':11,  'figure.facecolor': 'white'})

# ### Perform a t-test on classifier weights across groups
# pi3k mutant vs. pi3k wildtype
targene_mutant = ccle_full_df[ccle_full_df['targene_status'] == 1]
targene_wt = ccle_full_df[ccle_full_df['targene_status'] == 0]

# Also interested in BRAF status within Ras wildtype samples
BRAF_mutant = targene_wt[targene_wt['BRAF_MUT'] == 1]
BRAF_wt = targene_wt[targene_wt['BRAF_MUT'] == 0]

# Also interested in PTEN status within PI3k gain_mutant samples
BRAF_mutant_2 = targene_mutant[targene_mutant['BRAF_MUT'] == 1]
BRAF_wt_2 = targene_mutant[targene_mutant['BRAF_MUT'] == 0]

# Output t-test results
t_results_targene = ttest_ind(a = targene_mutant['weight'],
                          b = targene_wt['weight'], equal_var = False)
print('targene Status:')
print(t_results_targene)

t_results_BRAF = ttest_ind(a = BRAF_mutant['weight'],
                           b = BRAF_wt['weight'], equal_var = False)
print('\nBRAF Status in TARGENE Wild-Type Samples:')
print(t_results_BRAF)

t_results_BRAF_2 = ttest_ind(a = BRAF_mutant_2['weight'],
                           b = BRAF_wt_2['weight'], equal_var = False)
print('\nBRAF Status in TARGENE mutant-Type Samples:')
print(t_results_BRAF_2)

cell_line_folder = os.path.join(classifier, 'figures', 'cell_line')
if not os.path.exists(cell_line_folder):
    os.makedirs(cell_line_folder)

# Plot Results for targene and BRAF mutation status
x1, x2 = 0, 1
x3, x4 = -0.2, 0.2
x5, x6 = 0.8, 1.2
y1, y2,y3, h = 1.17, 1.0,1.02, 0.03

plt.rcParams['figure.figsize']=(3.5, 4)
ax1 = sns.boxplot(x="targene_status", y="weight", data=ccle_full_df,
                 hue='BRAF_MUT', palette = {0: "whitesmoke", 1: 'gainsboro'},
                 fliersize=0)
ax1 = sns.stripplot(x='targene_status', y='weight', hue='BRAF_MUT',
                   data=ccle_full_df, 
                   dodge=True, edgecolor='gray',
                   palette = {1: "seagreen", 0: 'goldenrod'},
                   jitter=0.25, size=2, alpha=0.65)
handles, labels = ax1.get_legend_handles_labels()
l = plt.legend(handles[2:4], ['Wild-Type', 'Mutant'], bbox_to_anchor=(.63, 0.2), loc=2, borderaxespad=0.)
l.set_title("BRAF")
ax1.axes.set_ylim(0, 1.3)
ax1.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1, ''])
ax1.set_xticklabels(['Targene Wild-Type', 'Targene Mutant'])
ax1.set_ylabel('Targene Classifier Score')
ax1.set_xlabel('CCLE Data')
ax1.legend
plt.axhline(0.5, color='black', linestyle='dashed', linewidth=1)

# Add targene T-Test Results
plt.plot([x1, x1, x2, x2], [y1, y1+h, y1+h, y1], lw=1.2, c='black')
plt.text(.5, y1+h, "{:.2E}".format(Decimal(t_results_targene.pvalue)),
         ha='center', va='bottom', color="black")

# Add BRAF t-test results
plt.plot([x3, x3, x4, x4], [y2, y2+h, y2+h, y2], lw=1.2, c='black')
plt.text(0, y2+h, "{:.2E}".format(Decimal(t_results_BRAF.pvalue)),
         ha='center', va='bottom', color="black")

# Add BRAF t-test results
plt.plot([x5, x5, x6, x6], [y3, y3+h, y3+h, y3], lw=1.2, c='black')
plt.text(1, y3+h, "{:.2E}".format(Decimal(t_results_BRAF_2.pvalue)),
         ha='center', va='bottom', color="black")

plt.tight_layout()
#ccle_fig_file = os.path.join('..', 'figures', 'cell_line', 'ccle_pi3k_GAIN_wild_mut_BRAF_predictions.pdf')
ccle_fig_file1 = os.path.join(classifier, 'figures', 'cell_line', 'ccle_targene_WT_MUT_BRAF_predictions.pdf')
plt.savefig(ccle_fig_file1)
plt.close()

import seaborn as sn
import matplotlib.pyplot as pl
# Plot Results for targene alone
x1, x2 = 0, 1
y1, y2,h = 1.05, 1.0, 0.03

plt.rcParams['figure.figsize']=(3.5, 4)
ax2 = sn.boxplot(x="targene_status", y="weight", data=ccle_full_df,
                 palette = {0: "lightgreen",1: 'yellow'},
                 fliersize=0)
ay2 = sn.stripplot(x='targene_status', y='weight', data=ccle_full_df, 
                   dodge=False,
                   palette = {0: "blue", 1: 'red'},
                   jitter=0.12, size=2, alpha=0.65)

ax2.axes.set_ylim(0, 1.2)
ax2.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1,''])
ax2.set_xticklabels(['Targene Wild-Type', 'Targene Mutant'])
ax2.set_ylabel('Targene Classifier Score')
ax2.set_xlabel('CCLE Data')
ax2.legend
pl.axhline(0.5, color='black', linestyle='dashed', linewidth=1)
ay2.axes.set_ylim(0, 1.2)

# Add PI3K T-Test Results
pl.plot([x1, x1, x2, x2], [y1, y1+h, y1+h, y1], lw=1.2, c='black')
pl.text(.6, y1+h, "{:.2E}".format(Decimal(t_results_targene.pvalue)),
         ha='center', va='bottom', color="black")

pl.tight_layout()
#ccle_fig_file = os.path.join('..', 'figures', 'cell_line', 'ccle_pi3k_GAIN_wild_mut_predictions.pdf')
ccle_fig_file2 = os.path.join(classifier, 'figures', 'cell_line', 'ccle_targene_WT_MUT_predictions.pdf')
pl.savefig(ccle_fig_file2)
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

# Of wild-type Ras cell lines, how many are predicted correctly?
# True Negative Rate, Specificity
negative_targene_lines_ccle['predictions'].value_counts()

# Of mutated Ras cell lines, how many are predicted correctly?
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

# Of the False positives, how many are BRAF mutant?
wt_targene_braf_ccle = positive_targene_predictions_ccle[positive_targene_predictions_ccle['targene_status'] == 0]
braf_neg, braf_pos = wt_targene_braf_ccle['BRAF_MUT'].value_counts()

print('{} of {} total false positives '
      'have BRAF mutations ({:.1f}%)'.format(braf_pos, false_pos,
                                             braf_pos * 100 / false_pos))

# If include BRAF mutations, how many correct
correct_braf = wt_targene_braf_ccle['BRAF_MUT'].value_counts()[1]
true_pos_with_braf = true_pos + correct_braf
print('Including BRAF mutants, {} of {} TARGENE mutation predictions '
      'have targene pathway mutations ({:.1f}%)'.format(true_pos_with_braf,
                                                    predicted_condition_pos,
                                                    true_pos_with_braf * 100 / predicted_condition_pos))

print('Of the false positives, there are {} BRAF mutated cell lines '
      'and {} BRAF wild-type cell lines'.format(braf_pos, braf_neg))

total_braf_wildtype, total_braf_mut = ccle_full_df['BRAF_MUT'].value_counts()
print('In all of CCLE, there are {} BRAF mutated cell lines '
      'and {} BRAF wild-type cell lines'.format(total_braf_mut, total_braf_wildtype))


# ### Add CCLE Variant Scores (nucleotide and amino acid) to Supplementary Data Files

# Load TCGA PanCanAtlas Core Ras Pathway genes
path_genes_file = args.path_genes or os.path.join('..', 'data', 'pi3k_genes.csv')
path_core_df = pd.read_table(path_genes_file)

# Subset MAF file to Ras pathway variants and merge with CCLE classifier scores
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

#updated_data_s5_file = os.path.join('..', 'classifiers', 'ERBB2_PIK3CA_KRAS_AKT1','tables', 'updated_Data_S5.csv')
updated_data_s5_file = os.path.join(classifier,'tables', 'updated_Data_S5.csv')
updated_data_s5_df.to_csv(updated_data_s5_file, sep=',', index=False)

# ## CCLE - Pharmacologic Efficacy
# 
# Here, we process drug efficacy data on the CCLE dataset. Data obtained from https://portals.broadinstitute.org/ccle/data (Pharmacologic profiling) (signin required).
# 
# A processed `.tsv` file is output to be visualized in `scripts/viz/ras_ccle_pharmacology.R`.

# Load in pharmacological results
pharm_file = args.phar_data or os.path.join('..', 'data', 'CCLE_NP24.2009_Drug_data_2015.02.24.csv')
pharm_df = pd.read_csv(pharm_file, index_col=0)
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

gg.ggplot(compound_heatmap, gg.aes('factor(tissue)', 'factor(Compound)', fill='count')) + gg.geom_tile(gg.aes(width=.95, height=.95)) + gg.theme(axis_text_x=gg.element_text(rotation=90),
         panel_background=gg.element_rect(fill='white'))
# p = os.path.join(classifier, 'figures', 'pharmacology_compound_heatmap.png')
# plt.savefig(p)
# Write out pharm_full_df to file to plot in ggplot2
# plotnine does not include all the functionality required to create the plot
#pharm_file = os.path.join('..', 'data', 'pharmacology_predictions_ERBB2_PIK3CA_KRAS_AKT1_ccle.tsv')

pharm_file = os.path.join(classifier, 'tables', 'pharmacology_predictions_targene_ccle.tsv')
pharm_full_df.to_csv(pharm_file, sep='\t')



