"""
Created by Juan M.
on 30/04/2020
"""

'''
Compiles .csv files generated in essential_nutrient_assessment.py from individual microbes into a single .csv
tables. The new file is saved in the same folder were individual files are read from.
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
path = ''

# makes a list of models/microbes that are in the consumption folder (models generated in script 1)
# this avoids any errors if for some reason an output file wasn't generated while running script 1
model_names = [f for f in listdir(path + 'growth_no-growth/') if isfile(join(path + 'growth_no-growth/', f))]
model_names = [os.path.splitext(f)[0] for f in model_names]

final_compilation_table_growth = pd.DataFrame()

for model in model_names:
    # model results are added to the yes/no growth compilation table
    microbe_growth = pd.read_csv(path + 'growth_no-growth/' + model + '.csv', index_col=[0])
    microbe_growth = microbe_growth.transpose()
    final_compilation_table_growth = pd.concat([final_compilation_table_growth, microbe_growth])

# Growth compilation table
final_compilation_table_growth.to_csv(path + 'compilation/growth_compilation.csv')
