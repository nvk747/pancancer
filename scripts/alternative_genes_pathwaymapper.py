#!/usr/bin/env python3
# Pancancer_Aberrant_Pathway_Activity_Analysis scripts/alternative_genes_pathwaymapper.py

import os
import sys
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from sklearn.metrics import roc_auc_score, average_precision_score

sys.path.insert(0, os.path.join(os.path.split(os.path.realpath(__file__))[0], 'util'))

from tcga_util import add_version_argument


def get_gene_auroc(x, w):
    score = roc_auc_score(x, w, average='weighted')
    return(score)

def get_gene_auprc(x, w):
    score = average_precision_score(x, w, average='weighted')
    return(score)

# argument passing 
parser = argparse.ArgumentParser()
add_version_argument(parser)
parser.add_argument('-s', '--scores',
                    help='string of the location of classifier scores/alt_folder')
parser.add_argument('-g', '--genes', default= 'KRAS,NRAS,HRAS',
                    help='string of the genes to extract or genelist file')
parser.add_argument('-p', '--path_genes',
                    help='pathway gene list file')
parser.add_argument( '--filename_mut', default=None,
                        help='Filename of sample/gene mutations to use in model')
parser.add_argument( '--filename_sample', default=None,
                        help='Filename of patient/samples to use in model')
parser.add_argument('-c', '--copy_number', action='store_true',
                    help='optional flag to include copy number info in pathway map')
parser.add_argument( '--filename_copy_loss', default=None,
                    help='Filename of copy number loss')
parser.add_argument( '--filename_copy_gain', default=None,
                    help='Filename of copy number gain')
args = parser.parse_args()

scores = args.scores
path_genes = args.path_genes
copy_number = args.copy_number

# if list of the genes provided by file or comma seperated values:
try:
    genes = args.genes
    genes_df = pd.read_table(genes)
    genes = genes_df['genes'].tolist()
except:
    genes = args.genes.split(',')

# if list of pathway genes are provided in a file
try:
    genes_df = pd.read_table(path_genes)
    path_genes = genes_df['genes'].tolist()
except:
    path_genes = path_genes.split(',')

mut_file = args.filename_mut or os.path.join('data', 'pancan_mutation_freeze.tsv')
sample_freeze_file = args.filename_sample or os.path.join('data', 'sample_freeze.tsv')
copy_loss_file = args.filename_copy_loss or os.path.join('data', 'copy_number_loss_status.tsv')
copy_gain_file = args.filename_copy_gain or os.path.join('data', 'copy_number_gain_status.tsv')

mutation_df = pd.read_table(mut_file, index_col=0)
sample_freeze = pd.read_table(sample_freeze_file, index_col=0)
copy_loss_df = pd.read_table(copy_loss_file, index_col=0)
copy_gain_df = pd.read_table(copy_gain_file, index_col=0)

# Load Pathway Genes
pathway_genes_file = args.path_genes or os.path.join('data', 'ras_genes.csv')
pathway_genes_df = pd.read_table(pathway_genes_file)

# Load classifier weights
targene_decision_file = os.path.join(scores, 'classifier_decisions.tsv')
targene_decisions_df = pd.read_table(targene_decision_file)

pathway_mutations_df = mutation_df[pathway_genes_df['genes']]

# Add status to the Y matrix depending on if the gene is a tumor suppressor
# or an oncogene. An oncogene can be activated with copy number gains, but
# a tumor suppressor is inactivated with copy number loss

oncogene = pathway_genes_df[pathway_genes_df['og_tsg'] == 'OG']
tumor_suppressor = pathway_genes_df[pathway_genes_df['og_tsg'] == 'TSG']

# Subset copy number information
pathway_copy_gain_sub_df = copy_gain_df[oncogene['genes']]
pathway_copy_loss_sub_df = copy_loss_df[tumor_suppressor['genes']]

# Combine Copy Number data
pathway_copy_df = pd.concat([pathway_copy_gain_sub_df, pathway_copy_loss_sub_df], axis=1)

pathway_status_df = pathway_mutations_df + pathway_copy_df
pathway_status_df[pathway_status_df == 2] = 1

subset_columns = ['SAMPLE_BARCODE', 'DISEASE', 'weight', 'total_status', 'log10_mut',
                  'hypermutated', 'include']
targene_decisions_subset_df = targene_decisions_df[subset_columns]
pathway_full_status_df = pathway_status_df.merge(targene_decisions_subset_df, left_index=True,
                                         right_on='SAMPLE_BARCODE')
pathway_full_status_df.index = pathway_full_status_df['SAMPLE_BARCODE']

# Remove hyper mutated samples
burden_filter = pathway_full_status_df['hypermutated'] == 0
burden_filter = burden_filter & pathway_full_status_df['log10_mut'] < 5 * pathway_full_status_df['log10_mut'].std()
pathway_full_status_df = pathway_full_status_df[burden_filter]

full_auroc = (
    pathway_full_status_df[pathway_genes_df['genes']]
    .apply(lambda x: get_gene_auroc(x, pathway_full_status_df['weight']))
    )

full_auprc = (
    pathway_full_status_df[pathway_genes_df['genes']]
    .apply(lambda x: get_gene_auprc(x, pathway_full_status_df['weight']))
    )

# Remove targene pathway positive samples, and recalculate metrics
#drop targene genes:
remove_targene_status = pathway_full_status_df[pathway_full_status_df['total_status'] == 0]
remove_targene_status_df = remove_targene_status[pathway_genes_df['genes']]
remove_targene_status_df = remove_targene_status_df.drop(genes, axis=1)
full_auroc_remove = remove_targene_status_df.apply(lambda x: get_gene_auroc(x, w=remove_targene_status['weight']))
full_auprc_remove = remove_targene_status_df.apply(lambda x: get_gene_auprc(x, w=remove_targene_status['weight']))

# Get output metrics for targene classification
output_pathway_metrics = pd.concat([full_auroc, full_auroc_remove], axis=1, sort=False)
output_pathway_metrics = output_pathway_metrics * 100  # To get percent
output_pathway_metrics = output_pathway_metrics - 50  # Subtract 50 from AUROC only

# Combine with AUPRC
output_pathway_metrics = pd.concat([output_pathway_metrics, full_auprc * 100,
                                full_auprc_remove * 100], axis=1, sort=False)
output_pathway_metrics.columns = ['pathway_auroc', 'no_targene_auroc', 'pathway_auprc', 'no_targene_auprc']

# Fill removed targene metrics with included metrics
output_pathway_metrics['no_targene_auroc'] = (
    output_pathway_metrics['no_targene_auroc'].fillna(output_pathway_metrics['pathway_auroc'])
    )
output_pathway_metrics['no_targene_auprc'] = (
    output_pathway_metrics['no_targene_auprc'].fillna(output_pathway_metrics['pathway_auprc'])
    )

# Write results to file
tables_folder = os.path.join(scores, 'tables')

if not os.path.exists(tables_folder):
    os.makedirs(tables_folder)

pathway_metric_file = os.path.join(scores, 'tables', 'pathway_metrics_pathwaymapper.txt')
output_pathway_metrics.to_csv(pathway_metric_file, sep='\t')

# Display targene pathway metrics
all_samples_targene_pathway_status = pathway_full_status_df[pathway_genes_df['genes']].max(axis=1)
print('targene Pathway Performance Summary: All pathway Genes')
print('AUROC:')
print(roc_auc_score(all_samples_targene_pathway_status,
                    pathway_full_status_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(all_samples_targene_pathway_status,
                              pathway_full_status_df['weight'], average='weighted'))

print('targene Pathway Performance Summary:', genes)
print('AUROC:')
print(roc_auc_score(pathway_full_status_df['total_status'],
                    pathway_full_status_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(pathway_full_status_df['total_status'],
                              pathway_full_status_df['weight'], average='weighted'))

print('targene Pathway Performance Summary: Held Out Samples')
held_out_pathway_df = pathway_full_status_df[pathway_full_status_df['include'] == 0]
print('AUROC:')
print(roc_auc_score(held_out_pathway_df['total_status'],
                    held_out_pathway_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(held_out_pathway_df['total_status'],
                              held_out_pathway_df['weight'], average='weighted'))

# # Visualize Distribution of AUROC and AUPRC for all genes

# Subset mutation file by samples
sub_full_mutation_df = mutation_df[burden_filter]
low_mutation_count_filter = (
    sub_full_mutation_df.sum()
    [sub_full_mutation_df.sum() >= 10].sort_values(ascending=False).index
    )
sub_full_mutation_df = sub_full_mutation_df[low_mutation_count_filter]
sub_full_mutation_df.head()

# Get Metrics for All Genes
all_auprc = sub_full_mutation_df.apply(lambda x: get_gene_auprc(x, w = pathway_full_status_df['weight']))
all_auroc = sub_full_mutation_df.apply(lambda x: get_gene_auroc(x, w = pathway_full_status_df['weight']))

# Process file and save results
all_gene_metrics_file = os.path.join(scores, 'tables', 'all_gene_metric_ranks.tsv')

all_genes_auprc_df = pd.DataFrame(all_auprc.sort_values(ascending=False), columns=['auprc'])
all_genes_auroc_df = pd.DataFrame(all_auroc.sort_values(ascending=False), columns=['auroc'])

all_genes_auprc_df = all_genes_auprc_df.assign(auprc_rank = list(range(0, all_genes_auprc_df.shape[0])))
all_genes_auroc_df = all_genes_auroc_df.assign(auroc_rank = list(range(0, all_genes_auprc_df.shape[0])))

all_genes_auprc_df = all_genes_auprc_df.assign(targene = 0)
all_genes_auprc_df.loc[all_genes_auprc_df.index.isin(pathway_genes_df['genes']), 'targene'] = 1

all_genes_metrics_df = all_genes_auprc_df.reset_index().merge(all_genes_auroc_df,
                                                              left_on='index', right_index=True)

all_genes_metrics_df.columns = ['Gene', 'AUPRC', 'AUPRC Rank', 'targene', 'AUROC', 'AUROC Rank']
all_genes_metrics_df.to_csv(all_gene_metrics_file, sep='\t', index=False)
