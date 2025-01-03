"""
Created by Juan M.
on 30/04/2020
"""
'''
This script compiles .csv files from individual microbes into single .csv table.
'''
import pandas as pd
import seaborn as sns
import numpy as np
import os
from os import listdir
from os.path import isfile, join

# Local path to individual .csv files to be read
path = 'C:/'

# make a list of model names (.csv files to be compiled)
model_names = [f for f in listdir(path + 'growth_no-growth/') if isfile(join(path + 'growth_no-growth/', f))]
model_names = [os.path.splitext(f)[0] for f in model_names]

final_compilation_table_growth = pd.DataFrame()

for model in model_names:
    # model results are added to the yes/no growth compilation table
    microbe_growth = pd.read_csv(path + 'growth_no-growth/' + model + '.csv', index_col=[0])
    microbe_growth = microbe_growth.transpose()
    final_compilation_table_growth = pd.concat([final_compilation_table_growth, microbe_growth])

# Growth compilation table
final_compilation_table_growth.to_csv(path_out + 'compilation/growth_compilation.csv')
