#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
from decimal import Decimal
from scipy.stats import ttest_ind
from statsmodels.stats.proportion import proportions_chisquare
from sklearn.preprocessing import StandardScaler
from Bio.SeqUtils import IUPACData
# from openpyxl import load_workbook
import matplotlib.pyplot as plt
import seaborn as sns
import plotnine as gg
import argparse

# Store protein change dictionary
aa = IUPACData.protein_letters_1to3_extended

#get_ipython().run_line_magic('matplotlib', 'inline')

parser = argparse.ArgumentParser()

parser.add_argument('-c', '--classifier', default= None,
                    help='location of classifier_summary file')
parser.add_argument('-e', '--ex_vlog',default= None,
                    help='path for external sample expression data file[fpkm/rlog/vlog')
parser.add_argument('-s', '--sign',
                    help='assigned tumor [1] or normal sample status[-1]')
args = parser.parse_args()

# load targene classifier summary file

classifier = args.classifier
classifier_file = os.path.join( classifier , "classifier_summary.txt")
all_coef_df = pd.read_table(os.path.join( classifier , "classifier_coefficients.tsv"), index_col=0)
# with open(classifier_file) as class_fh:
#    for line in class_fh:
#        line = line.strip().split('\t')
#        if line[0] == 'Coefficients:':
#            all_coef_df = pd.read_table(os.path.join(line[1]), index_col=0)

# Only non-zero coefficients contribute to model performance
coef_df = all_coef_df[all_coef_df['abs'] > 0]

# load external sample gene expression data: vlog or rlog or fpkm values
vlog_file = args.ex_vlog
vlog_df = pd.read_csv(vlog_file, index_col= 0)

# Determine the extent of coefficient overlap
common_genes = list(set(coef_df['feature']) & set(vlog_df.index))
common_coef = coef_df[coef_df['feature'].isin(common_genes)]
print('There are a total of {} out of {} genes in common between the datasets'
      .format(common_coef.shape[0], coef_df.shape[0]))

vlog_df = vlog_df.loc[common_coef['feature'], vlog_df.columns[0:]]

pd.set_option('display.max_rows', 500)
vlog_df = vlog_df[~vlog_df.index.duplicated(keep='first')]

# Which Genes are Missing?
missing_genes = list(set(coef_df['feature']).difference(set(vlog_df.index)))
all_coef_df[all_coef_df['feature'].isin(missing_genes)]

# Transform the cell line data by z-score
scaled_fit = StandardScaler().fit(vlog_df.T)
vlog_df = pd.DataFrame(scaled_fit.transform(vlog_df.T),
                                index=vlog_df.columns,
                                columns=vlog_df.index)
# Get the weights ready for applying the classifier
apply_weights = pd.DataFrame(common_coef['weight'])
apply_weights.index = common_coef.feature

# Apply a logit transform [y = 1/(1+e^(-wX))] to output probabilities
result = apply_weights.T.dot(vlog_df.T)
result = 1 / (1 + np.exp(-1 * result))

result2 = result.T.sort_values(by='weight')

result = result2.assign(name=result2.index)
result = result.sort_values(by='name')
print(result) # printing the result table

# load status of the external-sample tumors :+1 normal : -1
from csv import reader
opened_file = open(args.sign)
s = reader(opened_file)
status = list(s)
f_status = []
for i in status:
    n = int(i[0])
    f_status.append(n)
f_status[0]

# Tumor or normal status from RNAseq
output = result.assign(mut_status = f_status)
output = output.assign(sample_name = output.index)
output = output.assign(dummy_y = 0)
output

# Perform a t-test to determine if weights are significantly different
targene_geo_mutant = output[output['mut_status'] == 1]
targene_geo_wt = output[output['mut_status'] == -1]

# Output t-test results
t_results_geo_targene = ttest_ind(a = targene_geo_mutant['weight'],
                              b = targene_geo_wt['weight'], equal_var = False)
print('Statistic = {:.2f}, p = {:.2E}'.format(t_results_geo_targene[0],
                                              Decimal(t_results_geo_targene[1])))

# graphical output for predictions
p = (gg.ggplot(output,
               gg.aes(x='weight', y='dummy_y', color='factor(mut_status)')) +
     gg.geom_hline(gg.aes(yintercept=0), linetype='solid') +
     gg.geom_point(size=4) +
     gg.scale_color_manual(values=["#377eb8", "#ff7f00"], labels=['WT', 'Mutant']) +
     gg.ylim([-0.1, 0.1]) +
     gg.xlim([-0.001, 1.001]) +
     gg.theme_seaborn(style='whitegrid') +
     gg.xlab('Targene Classifier Score') +
     gg.ylab('') +
     gg.labs(color='Sample_status') +
     gg.ggtitle('Mutant vs WT \n') +
     gg.theme(
        plot_title=gg.element_text(size=22),
        axis_title_x=gg.element_text(size=16),
        axis_text_x=gg.element_text(size=16),
        axis_text_y=gg.element_blank(),
        axis_ticks_length=4,
        axis_ticks_major_y=gg.element_blank(),
        axis_ticks_minor_y=gg.element_blank(),
        axis_ticks_minor_x=gg.element_blank(),
        legend_position=(1.02, 0.8),
        legend_background=gg.element_blank(),
        legend_key=gg.element_rect(fill='white'),
        legend_text=gg.element_text(size=9),
        legend_title=gg.element_text(size=12),
        panel_border=gg.element_blank(),
        panel_grid_major=gg.element_blank(),
        panel_grid_minor=gg.element_blank()))
# targene_fig_file = os.path.join('..', 'figures', 'cell_line', 'targene_external_sample_predictions.pdf')
targene_fig_file = os.path.join(classifier, 'figures','targene_external_sample_predictions.pdf')
p.save(targene_fig_file, width=6, height=0.5)
p

# graphical output for predictions

from matplotlib.pyplot import figure
figure(num=None, figsize=(4, 4), dpi=300, facecolor='w', edgecolor='k')
x = targene_geo_mutant['weight']
y = targene_geo_wt['weight']

plt.title('Mutant vs WT')

sns.distplot(x, hist = False, kde = True,  rug=True, rug_kws={"color": "darkblue"},
                 kde_kws = {'shade': True, 'linewidth': 2, 'clip': (0.0, 1.0) }, 
                  label = 'Mutant', color = 'blue')
sns.distplot(y, hist = False, kde = True, rug=True, rug_kws={"color": "darkorange"},
                 kde_kws = {'shade': True, 'linewidth': 2, 'clip': (0.0, 1.0)},
                  label = 'WT', axlabel = 'Classifier Score', color = 'orange')

plt.xlim(left = -0.1)
plt.xlim(right = 1.1)
locs_x, labels_x = plt.xticks(np.arange(0,1.25,0.25))
plt.axvline(0.5, color='black', linestyle='dashed', linewidth=1)
plt.savefig(os.path.join(classifier, 'figures','targene_external_sample_predictions_1.pdf'))
plt.close()

l_x = len(x)
l_y = len(y)
xscore = 0
for j in x:
    if j > 0.5:
        xscore = xscore + 1
yscore = 0
for j in y:
    if j < 0.5:
        yscore = yscore +1
x_per = xscore/l_x * 100
y_per = yscore/l_y * 100

print('Stated as Tumor:',l_x)
print('Stated as control:',l_y)
print('no of tumor predicted as Mutant:',xscore)
print('no of controls predicted as WT:',yscore)
print('Accuracy for tumor samples:',x_per)
print('Accuracy for control samples:',y_per)
