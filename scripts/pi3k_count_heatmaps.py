
# coding: utf-8

# In[1]:


import os
import sys
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


# Get the current working directory
cwd = os.getcwd()

# Ensure that the path is starting in the scripts directory
if not cwd.split('/')[-1] == 'scripts':
    os.chdir(os.path.join(cwd, 'scripts'))


# In[3]:

# get_ipython().run_line_magic('matplotlib', 'inline')
plt.style.use('seaborn-notebook')


# In[4]:


# Load Datasets
print("executing pi3k_count_heatmaps and loading datasets")
mut_file = os.path.join('..', 'data', 'pancan_mutation_freeze.tsv.gz')
sample_freeze_file = os.path.join('..', 'data', 'sampleset_freeze.csv')
copy_loss_file = os.path.join('..', 'data', 'copy_number_loss_status.tsv.gz')
copy_gain_file = os.path.join('..', 'data', 'copy_number_gain_status.tsv.gz')
cancer_genes_file = os.path.join('..', 'data', 'vogelstein_cancergenes.tsv')

mutation_df = pd.read_table(mut_file, index_col=0)
sample_freeze = (pd.read_csv(sample_freeze_file, index_col=0)
                   .drop(['SQUAMOUS', 'COMPLETE'], axis='columns')
                )
copy_loss_df = pd.read_table(copy_loss_file, index_col=0
  )
copy_gain_df = pd.read_table(copy_gain_file, index_col=0)
cancer_genes_df = pd.read_table(cancer_genes_file)


# In[5]:


# Load pi3k Pathway Genes
results_path= os.path.join('..', 'classifiers', 'PI3K')

# Load pi3k Pathway Genes
genes_file = os.path.join('..', 'data', 'pi3k_genes_2.csv')
genes_df = pd.read_table(genes_file)


# In[6]:


genes_df.head()


# In[7]:


# Subset mutation data
print("mutation_df_from line 74")
print(mutation_df)
mutation_sub_df = mutation_df.loc[:, genes_df['genes']]
print("mutation_sub_df_line 74 only to pi3k pathway genes")
print(mutation_sub_df)

# In[8]:


# Find if the input genes are in this master list
print("cancer_genes_df from line 49")
print(cancer_genes_df)
genes_sub = cancer_genes_df[cancer_genes_df['Gene Symbol'].isin(genes_df['genes'])]
print("genes_sub_line 86")
print(genes_sub)

# In[9]:


# Add status to the Y matrix depending on if the gene is a tumor suppressor
# or an oncogene. An oncogene can be activated with copy number gains, but
# a tumor suppressor is inactivated with copy number loss
tumor_suppressor = genes_df[genes_df['og_tsg'] == 'TSG']
oncogene = genes_df[genes_df['og_tsg'] == 'OG']

# Subset copy number information
copy_loss_sub_df = copy_loss_df[tumor_suppressor['genes']]
copy_gain_sub_df = copy_gain_df[oncogene['genes']]

print("copy_loss_sub_df_line 102")
print(copy_loss_sub_df)

print("copy_gain_sub_df_line 105")
print(copy_gain_sub_df)

# ## Output Mutation, Copy Number, and Total Heatmap (Gene by Cancer-type)

# In[10]:


mutation_sub_total_df = mutation_sub_df.assign(Total=mutation_sub_df.max(axis=1))
mut_disease_df = mutation_sub_total_df.merge(sample_freeze, left_index=True,
                                             right_on='SAMPLE_BARCODE')
mut_heatmap_df = mut_disease_df.groupby('DISEASE').mean()

print("mut_heatmap_df _line 119")
print(mut_heatmap_df)
# In[11]:

print("mut_disease_df from line 114, merge of mutations_sub_total with sample_freeze")
print(mut_disease_df)
gene_avg = mut_disease_df.mean()
gene_avg.name = 'Total'
print("gene_avg _line 127, mean of mutations from disease")
print(gene_avg)

# In[12]:


mut_heatmap_df = mut_heatmap_df.append(gene_avg)
print("mut_heatmap_df_line 134 adjusted to gene_avg")
print(mut_heatmap_df)

# In[13]:


sns.set_style("whitegrid")
sns.heatmap(mut_heatmap_df, linewidths=0.2, linecolor='black',
            cmap='Blues_r', square=True, cbar=True)
plt.ylabel('Cancer Types', fontsize=16)
plt.xlabel('pi3k Pathway Genes', fontsize=16)
plt.savefig(os.path.join(results_path, 'mut_df.svg'))


# In[14]:


copy_df = pd.concat([copy_gain_sub_df, copy_loss_sub_df], axis=1)
copy_total_df = copy_df.assign(Total=copy_df.max(axis=1))
copy_disease_df = copy_total_df.merge(sample_freeze, left_index=True,
                                      right_on='SAMPLE_BARCODE')
copy_heatmap_df = copy_disease_df.groupby('DISEASE').mean()

print("copy_heatmap_df _line 156 disease mean")
print(copy_heatmap_df)
# In[15]:


copy_avg = copy_disease_df.mean()
copy_avg.name = 'Total'


# In[16]:


copy_heatmap_df = copy_heatmap_df.append(copy_avg)
print("copy_heatmap_df _line 169 copy avg")
print(copy_heatmap_df)

# In[17]:


sns.set_style("whitegrid")
sns.heatmap(copy_heatmap_df, linewidths=0.2, linecolor='black',
            cmap='Blues_r', square=True)
plt.ylabel('Cancer Types', fontsize=16)
plt.xlabel('pi3k Pathway Genes', fontsize=16)
plt.savefig(os.path.join(results_path, 'copy_df.svg'))


# In[18]:


# Combined heatmap
print("mutation_sub_df from line 75")
print(mutation_sub_df)
print("copy_df from line 149")
print(copy_df)
comb_heat = mutation_sub_df + copy_df
comb_heat[comb_heat == 2] = 1  # Replace duplicates with just one
print("comb_heat line 193 adding mutation and copy no for pi3k genes")
print(comb_heat)

# In[19]:


comb_heat_df = comb_heat.merge(sample_freeze, left_index=True, right_on='SAMPLE_BARCODE')
comb_heat_total_df = comb_heat_df.assign(Total=comb_heat_df.max(axis=1))
comb_heatmap_df = comb_heat_total_df.groupby('DISEASE').mean()

print("comb_heatmap_df _line 203 disease mean" )
print(comb_heatmap_df)
# In[20]:


comb_avg = comb_heat_total_df.mean()
comb_avg.name = 'Total'


# In[21]:


comb_heatmap_plot = comb_heatmap_df.append(comb_avg)
print("comb_heatmap_df _line 215 appened with comb_avg" )
print(comb_heatmap_df)

# In[22]:


sns.set_style("whitegrid")
sns.heatmap(comb_heatmap_plot, linewidths=0.2, linecolor='black',
            cmap='Blues_r', square=True)
plt.ylabel('Cancer Types', fontsize=16)
plt.xlabel('pi3k Pathway Genes', fontsize=16)
plt.tight_layout()
plt.savefig(os.path.join(results_path, 'combined_df.svg'))


# ## Generating Pathway Mapper Text Files

# In[23]:


summary_score = pd.DataFrame([mut_heatmap_df.ix['Total', :], copy_heatmap_df.ix['Total', :]])
print("summary_score _line 237 combining mut_heatmap_df and copy_heatmap_df")
print(summary_score)
# transposing the summaryscore
summary_score=summary_score.T
summary_score.columns=['mutation', 'copy_number']
summary_score=summary_score * 100
summary_score = summary_score.round(decimals = 1)
print("summary_score _line 244 after modification")
print(summary_score)

# In[24]:


mut_heatmap_df


# In[25]:

# ras pathway has 38 genes and pi3k pathway has 18 genes so this dimentions mismatches
# tum_sup_mult = pd.Series([1] * 34 + [-1] * 4 + [1])
tum_sup_mult = pd.Series([1] * 14 + [-1] * 4 + [1])
print("tum_sup_mult _line 258 creating a pandas Series")
print(tum_sup_mult)
tum_sup_mult.index = summary_score.index


# In[26]:


summary_score = summary_score.mul(tum_sup_mult, axis=0)
pathway_mapper_file = os.path.join(results_path, 'tables',
                                   'pathwaymapper_percentages.txt')
summary_score.to_csv(pathway_mapper_file, sep='\t')


# ## Output number of Ras events per sample

# In[27]:


decision_file = os.path.join(results_path, 'classifier_decisions.tsv')
decisions_df = pd.read_table(decision_file)
decisions_df.head()
print("decisions_df_line 280")
print(decisions_df)

# In[28]:

other_pi3k_df = mutation_sub_df.drop(['PIK3R1', 'PIK3CA', 'AKT1'], axis=1)
other_pi3k_copy_df = copy_df.drop(['PIK3R1', 'PIK3CA', 'AKT1'], axis=1)
other_pi3k_all_df = comb_heat_df.drop(['PIK3R1', 'PIK3CA', 'AKT1'], axis=1)
print("other_pi3k_all_df _line 288")
print(other_pi3k_all_df)

# In[29]:


total_pi3k_mutations = pd.DataFrame(other_pi3k_df.sum(axis=1), columns=['mutation_count'])
total_pi3k_copy_events = pd.DataFrame(other_pi3k_copy_df.sum(axis=1), columns=['copy_count'])
total_pi3k_all = pd.DataFrame(other_pi3k_all_df.sum(axis=1), columns=['all_count'])
total_pi3k_all.index = comb_heat_df['SAMPLE_BARCODE']


# In[30]:


# Define output summary of mutation, copy, and total counts per sample by Ras pathway
count_summary = (
    decisions_df[['SAMPLE_BARCODE', 'DISEASE', 'weight']]
    .merge(total_pi3k_mutations, left_on='SAMPLE_BARCODE', right_index=True)
    )
hyper_samples = decisions_df[decisions_df['hypermutated'] == 1]['SAMPLE_BARCODE']
count_summary.ix[count_summary['SAMPLE_BARCODE'].isin(hyper_samples),
                 'mutation_count'] = 'hyper'
count_summary.head()
print("count_summary _line 312")
print(count_summary)

# In[31]:


count_summary['mutation_count'].value_counts()


# In[32]:


count_summary = total_pi3k_copy_events.merge(count_summary, left_index=True,
    right_on='SAMPLE_BARCODE')
count_summary = total_pi3k_all.merge(count_summary, left_index=True,
    right_on='SAMPLE_BARCODE')
count_summary = (
    decisions_df[['SAMPLE_BARCODE', 'total_status']]
    .merge(count_summary, left_on='SAMPLE_BARCODE', right_on='SAMPLE_BARCODE')
    )
count_summary.head()
print("count_summary _line 333")
print(count_summary)

# In[33]:

count_summary_file = os.path.join(results_path, 'tables',
                                  'pi3k_events_per_sample.tsv')
count_summary.to_csv(count_summary_file, sep='\t', index=False)
print("pi3k_count_heatmaps_done")