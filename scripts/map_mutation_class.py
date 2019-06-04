"""
Gregory Way 2017
PanCancer Classifier
scripts/map_mutation_class.py

Merge per sample classifier scores with mutation types present in each sample

Usage: Run in command line:

    python scripts/map_mutation_class.py

with required command arguments:

    --scores     string indicating the location of all sample predictions
    --genes      comma separated string of gene symbols or file name

and optional command argument:

    --copy_number   decision to include copy number information into output

Output:
.tsv file of classifier scores according to mutation class
"""

import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--scores',
                    help='string of the location of classifier scores')
#parser.add_argument('-g', '--genes',
#                    help='string of the genes to extract or genelist file')
parser.add_argument('-p', '--path_genes',
                    help='pathway gene list file')
parser.add_argument('-c', '--copy_number', action='store_true',
                    help='optional flag to include copy number info in map')

parser.add_argument( '--filename_copy_loss', default=None,
                    help='Filename of copy number loss')
parser.add_argument( '--filename_copy_gain', default=None,
                    help='Filename of copy number gain')
parser.add_argument( '--filename_raw_mut', default=None,
                    help='Filename of raw mut MAF')

args = parser.parse_args()

scores = args.scores
path_genes = args.path_genes
copy_number = args.copy_number

# Load command arguments
prediction_file = os.path.join(scores, 'classifier_decisions.tsv')

try:
    genes_df = pd.read_csv(path_genes)
    print("path_genes_file")
    print(genes_df)
    path_genes = genes_df['genes'].tolist()
    print(path_genes)
except:
    path_genes = path_genes.split(',')

out_dir = os.path.join(os.path.dirname(prediction_file), 'tables')

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

out_file = os.path.join(out_dir, 'mutation_classification_scores.tsv')
raw_mut_file = args.filename_raw_mut

pred_df = pd.read_table(prediction_file, index_col=0)
mut_df = pd.read_table(raw_mut_file)

# Process mutation file
mut_df = mut_df.assign(ID=mut_df.Tumor_Sample_Barcode.str.slice(start=0,
                                                                stop=15))
sub_mut_df = mut_df[mut_df['Hugo_Symbol'].isin(path_genes)]
sub_mut_df = sub_mut_df[['ID', 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'HGVSc',
                         'HGVSp', 'Variant_Classification']]

map_df = pred_df.merge(sub_mut_df, left_index=True, right_on='ID', how='outer')
map_df.index = map_df['ID']
map_df['Variant_Classification'] = map_df['Variant_Classification']\
                                         .fillna('Wild-Type')
if copy_number:
    # Load Copy Number info
    copy_loss_file = args.filename_copy_loss 
    copy_gain_file = args.filename_copy_gain 

    copy_loss_df = pd.read_table(copy_loss_file, index_col=0)
    copy_gain_df = pd.read_table(copy_gain_file, index_col=0)

    tumor_suppressor = genes_df[genes_df['og_tsg'] == 'TSG']
    oncogene = genes_df[genes_df['og_tsg'] == 'OG']

    # Subset copy number information
    copy_loss_sub_df = copy_loss_df[tumor_suppressor['genes']]
    copy_gain_sub_df = copy_gain_df[oncogene['genes']]

    # Melt dataframes to merge with output file
    copy_loss_melt = pd.melt(copy_loss_sub_df.reset_index(), id_vars='index',
                             var_name='copy_loss', value_name='indicator')
    copy_loss_melt = copy_loss_melt[copy_loss_melt['indicator'] != 0]
    copy_loss_melt.drop('indicator', axis=1, inplace=True)

    copy_gain_melt = pd.melt(copy_gain_sub_df.reset_index(), id_vars='index',
                             var_name='copy_gain', value_name='indicator')
    copy_gain_melt = copy_gain_melt[copy_gain_melt['indicator'] != 0]
    copy_gain_melt.drop('indicator', axis=1, inplace=True)
    copy_df = copy_loss_melt.merge(copy_gain_melt, left_on='index',
                                   right_on='index', how='outer')

    # Merge with mapped dataframe
    map_df = map_df.merge(copy_df, left_on='ID', right_on='index', how='outer')

map_df.to_csv(out_file, sep='\t')
print("map_mutation class done")