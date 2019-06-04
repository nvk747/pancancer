
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from sklearn.metrics import roc_auc_score, average_precision_score
sys.path.insert(0, os.path.join('scripts', 'util'))
from tcga_util import get_args


# In[2]:


# Get the current working directory
cwd = os.getcwd()

# Ensure that the path is starting in the scripts directory
if not cwd.split('/')[-1] == 'scripts':
    sys.path.append(os.path.join(cwd, 'scripts'))


# In[3]:


def get_gene_auroc(x, w):
    score = roc_auc_score(x, w, average='weighted')
    return(score)

def get_gene_auprc(x, w):
    score = average_precision_score(x, w, average='weighted')
    return(score)


# In[4]:

# setting abc directory [abc_dir='classifiers/PI3K']

args = get_args()
alt_folder = args.alt_folder
if alt_folder != 'Auto':
    abc_folder = alt_folder

#abc_folder = os.path.join('..', 'classifiers', 'PI3K_LOSS')


# In[5]:


# Load Datasets

#mut_file = os.path.join('..', 'data', 'pancan_mutation_freeze.tsv')
#sample_freeze_file = os.path.join('..', 'data', 'sample_freeze.tsv')
#copy_loss_file = os.path.join('..', 'data', 'copy_number_loss_status.tsv')
#copy_gain_file = os.path.join('..', 'data', 'copy_number_gain_status.tsv')

mut_file = os.path.join('data', 'pancan_mutation_freeze.tsv')
sample_freeze_file = os.path.join('data', 'sample_freeze.tsv')
copy_loss_file = os.path.join('data', 'copy_number_loss_status.tsv')
copy_gain_file = os.path.join('data', 'copy_number_gain_status.tsv')


mutation_df = pd.read_table(mut_file, index_col=0)
sample_freeze = pd.read_table(sample_freeze_file, index_col=0)
copy_loss_df = pd.read_table(copy_loss_file, index_col=0)
copy_gain_df = pd.read_table(copy_gain_file, index_col=0)


# In[6]:


# Load pi3k Pathway Genes
pathway_genes_file = os.path.join('data', 'pi3k_genes.csv')
pathway_genes_df = pd.read_csv(pathway_genes_file)

# In[7]:


# Load classifier weights
abc_decision_file = os.path.join(abc_folder, 'classifier_decisions.tsv')
abc_decisions_df = pd.read_table(abc_decision_file)

# In[8]:


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


# In[9]:

pathway_status_df = pathway_mutations_df + pathway_copy_df
pathway_status_df[pathway_status_df == 2] = 1


# In[10]:


subset_columns = ['SAMPLE_BARCODE', 'DISEASE', 'weight', 'total_status', 'log10_mut',
                  'hypermutated', 'include']
abc_decisions_subset_df = abc_decisions_df[subset_columns]
pathway_full_status_df = pathway_status_df.merge(abc_decisions_subset_df, left_index=True,
                                         right_on='SAMPLE_BARCODE')
pathway_full_status_df.index = pathway_full_status_df['SAMPLE_BARCODE']


# In[11]:


# Remove hyper mutated samples
burden_filter = pathway_full_status_df['hypermutated'] == 0
burden_filter = burden_filter & pathway_full_status_df['log10_mut'] < 5 * pathway_full_status_df['log10_mut'].std()
pathway_full_status_df = pathway_full_status_df[burden_filter]

# In[12]:


full_auroc = (
    pathway_full_status_df[pathway_genes_df['genes']]
    .apply(lambda x: get_gene_auroc(x, pathway_full_status_df['weight']))
    )

full_auprc = (
    pathway_full_status_df[pathway_genes_df['genes']]
    .apply(lambda x: get_gene_auprc(x, pathway_full_status_df['weight']))
    )


# In[13]:


# Remove pi3k positive samples, and recalculate metrics
#drop abc genes:
genes_df = pd.read_csv("/data/vijay/git/pancancer/data/abc_genes.csv")
genes = genes_df['genes'].tolist()
remove_abc_status = pathway_full_status_df[pathway_full_status_df['total_status'] == 0]
remove_abc_status_df = remove_abc_status[pathway_genes_df['genes']]
remove_abc_status_df = remove_abc_status_df.drop(genes, axis=1)
full_auroc_remove = remove_abc_status_df.apply(lambda x: get_gene_auroc(x, w=remove_abc_status['weight']))
full_auprc_remove = remove_abc_status_df.apply(lambda x: get_gene_auprc(x, w=remove_abc_status['weight']))


# In[16]:


# Get output metrics for pi3k classification
output_pathway_metrics = pd.concat([full_auroc, full_auroc_remove], axis=1, sort=False)
output_pathway_metrics = output_pathway_metrics * 100  # To get percent
output_pathway_metrics = output_pathway_metrics - 50  # Subtract 50 from AUROC only

# Combine with AUPRC
output_pathway_metrics = pd.concat([output_pathway_metrics, full_auprc * 100,
                                full_auprc_remove * 100], axis=1, sort=False)
output_pathway_metrics.columns = ['pathway_auroc', 'no_abc_auroc', 'pathway_auprc', 'no_abc_auprc']

# Fill removed pi3k metrics with included metrics
output_pathway_metrics['no_abc_auroc'] = (
    output_pathway_metrics['no_abc_auroc'].fillna(output_pathway_metrics['pathway_auroc'])
    )
output_pathway_metrics['no_abc_auprc'] = (
    output_pathway_metrics['no_abc_auprc'].fillna(output_pathway_metrics['pathway_auprc'])
    )

# Write results to file
tables_folder = os.path.join(abc_folder, 'tables')

if not os.path.exists(tables_folder):
    os.makedirs(tables_folder)

pathway_metric_file = os.path.join(abc_folder, 'tables', 'pathway_metrics_pathwaymapper.txt')
output_pathway_metrics.to_csv(pathway_metric_file, sep='\t')


# In[17]:


# Display pi3k pathway metrics
all_samples_abc_pathway_status = pathway_full_status_df[pathway_genes_df['genes']].max(axis=1)
print('abc Pathway Performance Summary: All pathway Genes')
print('AUROC:')
print(roc_auc_score(all_samples_abc_pathway_status,
                    pathway_full_status_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(all_samples_abc_pathway_status,
                              pathway_full_status_df['weight'], average='weighted'))


# In[18]:

print('abc Pathway Performance Summary:', genes)
print('AUROC:')
print(roc_auc_score(pathway_full_status_df['total_status'],
                    pathway_full_status_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(pathway_full_status_df['total_status'],
                              pathway_full_status_df['weight'], average='weighted'))


# In[19]:

print('abc Pathway Performance Summary: Held Out Samples')
held_out_pathway_df = pathway_full_status_df[pathway_full_status_df['include'] == 0]
print('AUROC:')
print(roc_auc_score(held_out_pathway_df['total_status'],
                    held_out_pathway_df['weight'], average='weighted'))
print('AUPRC:')
print(average_precision_score(held_out_pathway_df['total_status'],
                              held_out_pathway_df['weight'], average='weighted'))


# # Visualize Distribution of AUROC and AUPRC for all genes

# In[20]:


# Subset mutation file by samples
sub_full_mutation_df = mutation_df[burden_filter]
low_mutation_count_filter = (
    sub_full_mutation_df.sum()
    [sub_full_mutation_df.sum() >= 10].sort_values(ascending=False).index
    )
sub_full_mutation_df = sub_full_mutation_df[low_mutation_count_filter]
sub_full_mutation_df.head()


# In[21]:


# Get Metrics for All Genes
all_auprc = sub_full_mutation_df.apply(lambda x: get_gene_auprc(x, w = pathway_full_status_df['weight']))
all_auroc = sub_full_mutation_df.apply(lambda x: get_gene_auroc(x, w = pathway_full_status_df['weight']))


# In[23]:


# Process file and save results
all_gene_metrics_file = os.path.join(abc_folder, 'tables', 'all_gene_metric_ranks.tsv')

all_genes_auprc_df = pd.DataFrame(all_auprc.sort_values(ascending=False), columns=['auprc'])
all_genes_auroc_df = pd.DataFrame(all_auroc.sort_values(ascending=False), columns=['auroc'])

all_genes_auprc_df = all_genes_auprc_df.assign(auprc_rank = list(range(0, all_genes_auprc_df.shape[0])))
all_genes_auroc_df = all_genes_auroc_df.assign(auroc_rank = list(range(0, all_genes_auprc_df.shape[0])))

all_genes_auprc_df = all_genes_auprc_df.assign(abc = 0)
all_genes_auprc_df.loc[all_genes_auprc_df.index.isin(pathway_genes_df['genes']), 'abc'] = 1

all_genes_metrics_df = all_genes_auprc_df.reset_index().merge(all_genes_auroc_df,
                                                              left_on='index', right_index=True)

all_genes_metrics_df.columns = ['Gene', 'AUPRC', 'AUPRC Rank', 'abc', 'AUROC', 'AUROC Rank']
all_genes_metrics_df.to_csv(all_gene_metrics_file, sep='\t', index=False)

print("all_gene_metric_ranks file generated")
print("abc_alternative_genes_pathwaymapper done")