
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser()

#parser.add_argument('-c', '--classifier', default=None,
#                    help='folder location of classifier file')
#parser.add_argument('-u', '--copy_number', action='store_true',
#                    help='Supplement Y matrix with copy number events')
parser.add_argument('-g', '--genes',
                    help='gene list file')
parser.add_argument('-p', '--path_genes',
                    help='pathway gene list file')
parser.add_argument('-o', '--alt_folder', default='Auto',
                        help='Provide an alternative folder to save results')
parser.add_argument('-x', '--x_matrix', default=None,
                    help='Filename of features to use in model')
parser.add_argument( '--filename_mut', default=None,
                    help='Filename of sample/gene mutations to use in model')
parser.add_argument( '--filename_mut_burden', default=None,
                    help='Filename of sample mutation burden to use in model')
parser.add_argument( '--filename_sample', default=None,
                    help='Filename of patient/samples to use in model')
parser.add_argument( '--filename_copy_loss', default=None,
                    help='Filename of copy number loss')
parser.add_argument( '--filename_copy_gain', default=None,
                    help='Filename of copy number gain')
parser.add_argument( '--filename_cancer_gene_classification', default=None,
                    help='Filename of cancer gene classification table')

args = parser.parse_args()

# Load command arguments
#classifier_file = os.path.join(args.classifier, 'classifier_summary.txt')
#copy_number = args.copy_number

# Load Constants
#rnaseq_file = args.x_matrix or os.path.join('data', 'pancan_rnaseq_freeze.tsv')
#mut_file = args.filename_mut or os.path.join('data', 'pancan_mutation_freeze.tsv')
#sample_freeze_file = args.filename_sample or os.path.join('data', 'sample_freeze.tsv')
#cancer_gene_file = args.filename_cancer_gene_classification or os.path.join('data', 'vogelstein_cancergenes.tsv')
#copy_loss_file = args.filename_copy_loss or os.path.join('data', 'copy_number_loss_status.tsv')
#copy_gain_file = args.filename_copy_gain or os.path.join('data', 'copy_number_gain_status.tsv')
#mutation_burden_file = args.filename_mut_burden or os.path.join('data', 'mutation_burden_freeze.tsv')

alt_folder = args.alt_folder
rnaseq_file = args.x_matrix 
mut_file = args.filename_mut 
sample_freeze_file = args.filename_sample  
cancer_gene_file = args.filename_cancer_gene_classification 
copy_loss_file = args.filename_copy_loss 
copy_gain_file = args.filename_copy_gain 
mutation_burden_file = args.filename_mut_burden 

# In[2]:

# Get the current working directory
#cwd = os.getcwd()

# Ensure that the path is starting in the scripts directory
#if not cwd.split('/')[-1] == 'scripts':
#    os.chdir(os.path.join(cwd, 'scripts'))

# In[3]:

#get_ipython().run_line_magic('matplotlib', 'inline')
#plt.style.use('seaborn-notebook')

# In[4]:

# Load Datasets
#mut_file = os.path.join('..', 'data', 'pancan_mutation_freeze.tsv.gz')
#sample_freeze_file = os.path.join('..', 'data', 'sample_freeze.tsv')
#copy_loss_file = os.path.join('..', 'data', 'copy_number_loss_status.tsv.gz')
#copy_gain_file = os.path.join('..', 'data', 'copy_number_gain_status.tsv.gz')
#cancer_genes_file = os.path.join('..', 'data', 'vogelstein_cancergenes.tsv')

mutation_df = pd.read_table(mut_file, index_col=0)
mutation_df.head()
sample_freeze = pd.read_table(sample_freeze_file, index_col=0)
copy_loss_df = pd.read_table(copy_loss_file, index_col=0)
copy_gain_df = pd.read_table(copy_gain_file, index_col=0)
cancer_genes_df = pd.read_table(cancer_gene_file)

# In[5]:

# Load Ras Pathway Genes
#results_path = os.path.join('..', 'classifiers', 'RAS')
results_path = alt_folder

# Load Ras Pathway Genes
#genes_file = os.path.join('..', 'data', 'ras_genes.csv')
genes_file = args.genes
genes_df = pd.read_table(genes_file)
genes = genes_df['genes'].tolist()

pathgenes_file = args.path_genes
pathgenes_df = pd.read_csv(pathgenes_file)
print(pathgenes_df)

n = pathgenes_df['og_tsg'].tolist()
n_OG = n.count('OG')
n_TSG = n.count('TSG')
print(n_OG)
print(n_TSG)

# In[6]:

genes_df.head()
pathgenes_df.head()

# In[7]:

# Subset mutation data
mutation_sub_df = mutation_df.loc[:, pathgenes_df['genes']]

# In[8]:

# Find if the input genes are in this master list
genes_sub = cancer_genes_df[cancer_genes_df['Gene Symbol'].isin(pathgenes_df['genes'])]

# In[9]:

# Add status to the Y matrix depending on if the gene is a tumor suppressor
# or an oncogene. An oncogene can be activated with copy number gains, but
# a tumor suppressor is inactivated with copy number loss
tumor_suppressor = pathgenes_df[pathgenes_df['og_tsg'] == 'TSG']
oncogene = pathgenes_df[pathgenes_df['og_tsg'] == 'OG']

# Subset copy number information
copy_loss_sub_df = copy_loss_df[tumor_suppressor['genes']]
copy_gain_sub_df = copy_gain_df[oncogene['genes']]

# ## Output Mutation, Copy Number, and Total Heatmap (Gene by Cancer-type)

# In[10]:

mutation_sub_total_df = mutation_sub_df.assign(Total=mutation_sub_df.max(axis=1))
mut_disease_df = mutation_sub_total_df.merge(sample_freeze, left_index=True,
                                             right_on='SAMPLE_BARCODE')
mut_heatmap_df = mut_disease_df.groupby('DISEASE').mean()

# In[11]:

gene_avg = mut_disease_df.mean()
gene_avg.name = 'Total'

# In[12]:

mut_heatmap_df = mut_heatmap_df.append(gene_avg)

# In[13]:

sns.set_style("whitegrid")
plt.figure(figsize = (10,10),dpi= 300)
sns.heatmap(mut_heatmap_df, linewidths=0.2, linecolor='black', 
            cmap='Blues_r', square=True, cbar=True)
plt.autoscale(enable=True, axis ='x', tight = True)
plt.autoscale(enable=True, axis ='y', tight = True)
plt.ylabel('Cancer Types', fontsize=16)
plt.xlabel('Pathway Genes', fontsize=16)
plt.savefig(os.path.join(results_path, 'mut_df.pdf'))
plt.savefig(os.path.join(results_path, 'mut_df.jpeg'))

# In[14]:

copy_df = pd.concat([copy_gain_sub_df, copy_loss_sub_df], axis=1)
copy_total_df = copy_df.assign(Total=copy_df.max(axis=1))
copy_disease_df = copy_total_df.merge(sample_freeze, left_index=True,
                                      right_on='SAMPLE_BARCODE')
copy_heatmap_df = copy_disease_df.groupby('DISEASE').mean()

# In[15]:

copy_avg = copy_disease_df.mean()
copy_avg.name = 'Total'

# In[16]:

copy_heatmap_df = copy_heatmap_df.append(copy_avg)

# In[17]:

sns.set_style("whitegrid")
plt.figure(figsize = (10,10),dpi= 300)
sns.heatmap(copy_heatmap_df, linewidths=0.2, linecolor='black', 
            cmap='Blues_r', square=True)
plt.ylabel('Cancer Types', fontsize=16)
plt.xlabel('Pathway Genes', fontsize=16)
plt.autoscale(enable=True, axis ='x', tight = True)
plt.autoscale(enable=True, axis ='y', tight = True)
plt.savefig(os.path.join(results_path, 'copy_df.pdf'))
plt.savefig(os.path.join(results_path, 'copy_df.jpeg'))

# In[18]:

# Combined heatmap
comb_heat = mutation_sub_df + copy_df
comb_heat[comb_heat == 2] = 1  # Replace duplicates with just one

# In[19]:

comb_heat_df = comb_heat.merge(sample_freeze, left_index=True, right_on='SAMPLE_BARCODE')
comb_heat_total_df = comb_heat_df.assign(Total=comb_heat_df.max(axis=1))
comb_heatmap_df = comb_heat_total_df.groupby('DISEASE').mean()

# In[20]:

comb_avg = comb_heat_total_df.mean()
comb_avg.name = 'Total'

# In[21]:

comb_heatmap_plot = comb_heatmap_df.append(comb_avg)

# In[22]:

sns.set_style("whitegrid")
plt.figure(figsize = (10,10),dpi= 300)
sns.heatmap(comb_heatmap_plot, linewidths=0.2, linecolor='black', 
            cmap='Blues_r', square=True)
plt.ylabel('Cancer Types', fontsize=16)
plt.xlabel('Pathway Genes', fontsize=16)
plt.autoscale(enable=True, axis ='x', tight = True)
plt.autoscale(enable=True, axis ='y', tight = True)
plt.tight_layout()
plt.savefig(os.path.join(results_path, 'combined_df.pdf'))
plt.savefig(os.path.join(results_path, 'combined_df.jpeg'))

# ## Generating Pathway Mapper Text Files

# In[23]:

summary_score_df = (
    pd.DataFrame(
        [mut_heatmap_df.loc['Total', :], copy_heatmap_df.loc['Total', :]]
    )
    .transpose()
)
summary_score_df.columns = ['mutation', 'copy_number']
summary_score_df = summary_score_df * 100
summary_score_df = summary_score_df.round(decimals = 1)

# In[24]:
# Create negative percentages for tumor suppressors in the Pathway
tum_sup_mult = pd.Series([1] * n_OG + [-1] * n_TSG + [1])
tum_sup_mult.index = summary_score_df.index

# In[25]:

summary_score_df = summary_score_df.mul(tum_sup_mult, axis=0)
pathway_mapper_file = os.path.join(results_path, 'tables',
                                   'pathwaymapper_percentages.txt')
summary_score_df.to_csv(pathway_mapper_file, sep='\t')

# ## Output number of Ras events per sample

# In[26]:
decision_file = os.path.join(results_path, 'classifier_decisions.tsv')
decisions_df = pd.read_table(decision_file)
decisions_df.head()

# In[27]:
#other_ras_df = mutation_sub_df.drop(['KRAS', 'HRAS', 'NRAS'], axis=1)
#other_ras_copy_df = copy_df.drop(['KRAS', 'HRAS', 'NRAS'], axis=1)
#other_ras_all_df = comb_heat_df.drop(['KRAS', 'HRAS', 'NRAS'], axis=1)

other_genes_df = mutation_sub_df.drop(genes, axis=1)
other_genes_copy_df = copy_df.drop(genes, axis=1)
other_genes_all_df = comb_heat_df.drop(genes, axis=1)

# In[28]:

total_genes_mutations = pd.DataFrame(other_genes_df.sum(axis=1), columns=['mutation_count'])
total_genes_copy_events = pd.DataFrame(other_genes_copy_df.sum(axis=1), columns=['copy_count'])
total_genes_all = pd.DataFrame(other_genes_all_df.sum(axis=1), columns=['all_count'])
total_genes_all.index = comb_heat_df['SAMPLE_BARCODE']

# In[29]:

# Define output summary of mutation, copy, and total counts per sample by Ras pathway
count_summary = (
    decisions_df[['SAMPLE_BARCODE', 'DISEASE', 'weight']]
    .merge(total_genes_mutations, left_on='SAMPLE_BARCODE', right_index=True)
    )
hyper_samples = decisions_df[decisions_df['hypermutated'] == 1]['SAMPLE_BARCODE']
count_summary.loc[count_summary['SAMPLE_BARCODE'].isin(hyper_samples),
                  'mutation_count'] = 'hyper'
count_summary.head()

# In[30]:

count_summary['mutation_count'].value_counts()

# In[31]:

count_summary = total_genes_copy_events.merge(count_summary, left_index=True,
                                            right_on='SAMPLE_BARCODE')
count_summary = total_genes_all.merge(count_summary, left_index=True,
                                    right_on='SAMPLE_BARCODE')
count_summary = (
    decisions_df[['SAMPLE_BARCODE', 'total_status']]
    .merge(count_summary, left_on='SAMPLE_BARCODE', right_on='SAMPLE_BARCODE')
    )
count_summary.head()

# In[32]:

count_summary_file = os.path.join(results_path, 'tables',
                                  'path_events_per_sample.tsv')
count_summary.to_csv(count_summary_file, sep='\t', index=False)

print("abc_pathway_count_heatmaps_generated and done")