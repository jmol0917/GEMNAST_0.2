"""
Created by Juan M.
on 30/04/2020
"""

'''
Compiles .csv files generated in combinatorial_analysis.py (script 1) from individual microbes into single .csv
tables. Then generates charts and other relevant files for further analysis.
'''
import pandas as pd
import seaborn as sns
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt
import scipy.stats as ss
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

# Input (directory where individual strain files are saved) and Output directory (where Boolean table will be saved)
path_out = ''

# makes a list of models/microbes that are in the consumption folder (models generated in script 1)
# this avoids any errors if for some reason an output file wasn't generated while running script 1
model_names = [f for f in listdir(path_out + 'growth_no-growth/') if isfile(join(path_out + 'growth_no-growth/', f))]
model_names = [os.path.splitext(f)[0] for f in model_names]

final_compilation_table_growth = pd.DataFrame()
final_compilation_table_growth_rate = pd.DataFrame()
final_compilation_table_genome_growth = pd.DataFrame()

for model in model_names:
    # model results are added to the yes/no growth compilation table
    microbe_growth = pd.read_csv(path_out + 'growth_no-growth/' + model + '.csv', index_col=[0])
    microbe_growth = microbe_growth.transpose()
    final_compilation_table_growth = pd.concat([final_compilation_table_growth, microbe_growth])

# Growth compilation table
final_compilation_table_growth.to_csv(path_out + 'compilation/growth_compilation.csv')
sns.set(font_scale=0.25)
fig1 = sns.clustermap(final_compilation_table_growth, cmap="YlGnBu", cbar_pos=None, col_cluster=False,
                      figsize=(60, 160), dendrogram_ratio=(.01, .02))
fig1.savefig(path_out + '/Graphs/growth_compilation.jpeg', pil_kwargs={'optimize': True})
reordered_index = fig1.dendrogram_row.reordered_ind
reordered = fig1.data2d
reordered.insert(0, 'Original_Index', reordered_index, True)
reordered.to_csv(path_out + 'compilation/growth_compilation_reordered.csv')
