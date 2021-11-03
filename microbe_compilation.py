"""
Created by Juan M.
on 30/04/2020
"""
'''
Script 2 out of 3: 
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

# rename this variable based on the dimension being assessed
dimension = 'requirements_no_ubiquinone/consumption/'

# to run a comparison analysis between genome size and number of instances (growth/production)
regression_instances_vs_genome = False

# to run a comparison analysis between genome size and growth rate (same result no matter the dimension/aspect)
regression_growth_rate_vs_genome = False
growth_column = -1        # column location in relation to index from a boolean table that has recorded growth rates

# Path where simulation results and media are saved
# path_out = 'results/AGORA/'
# path_in = 'sbml/'

# Local path options
path_out = 'C:/Users/jpmo_/Dropbox (Sydney Uni)/gsm_analysis/results/AGORA/' + dimension

# make a list of models/microbes that are in the consumption folder (models generated in script 1)
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

    # model results are added to the growth rate at each media formulation compilation table
    if regression_growth_rate_vs_genome:
        microbe_growth_rate = pd.read_csv(path_out + 'growth_rate/' + model + '_growth_rate.csv')
        microbe_growth_rate = microbe_growth_rate.transpose()
        final_compilation_table_growth_rate = pd.concat([final_compilation_table_growth_rate, microbe_growth_rate])

    # model results are added to the growth instances vs genome size compilation table
    if regression_instances_vs_genome:
        microbe_genome_vs_growth = pd.read_csv(path_out + 'growth_instances/' + model + '_genome_vs_growth.csv')
        microbe_genome_vs_growth = microbe_genome_vs_growth.transpose()
        final_compilation_table_genome_growth = pd.concat([final_compilation_table_genome_growth,
                                                           microbe_genome_vs_growth])

'''
Compilation tables for growth yes/no, growth rate at every instance and genome vs total number of growth occurrences
are exported as csv. files and then corresponding charts are created
 '''
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

# Growth rates ( instead of 1s (ones)) compilation table
if regression_growth_rate_vs_genome:
    final_compilation_table_growth_rate.to_csv(path_out + 'compilation/growth_rate_compilation.csv')

    fig2 = sns.clustermap(final_compilation_table_growth_rate, cmap="YlGnBu", cbar_pos=(0, .2, .03, .4), col_cluster=False,
                          figsize=(150, 400), dendrogram_ratio=(.1, .2))
    fig2.savefig(path_out + 'consumption/Graphs/growth_rate_compilation.jpeg', optimize=True)

    reordered_index2 = fig2.dendrogram_row.reordered_ind
    reordered2 = fig2.data2d
    reordered2.insert(0, 'Original_Index', reordered_index2, True)
    reordered2.to_csv(path_out + 'compilation/growth_rate_compilation_reordered.csv')

# Growth vs growth compilation table
if regression_instances_vs_genome:
    final_compilation_table_genome_growth.to_csv(path_out + 'compilation/genome_growth_compilation.csv')

if regression_instances_vs_genome:
    X1 = final_compilation_table_genome_growth.iloc[:, 1].to_numpy()
    Y1 = final_compilation_table_genome_growth.iloc[:, 0].to_numpy()

    f, ax = plt.subplots(constrained_layout=True)
    xy = np.vstack([X1, Y1])
    z = ss.gaussian_kde(xy)(xy)
    ax.scatter(X1, Y1, c=z, s=40)
    plt.tight_layout()
    f.savefig(path_out + dimension + 'Graphs/genome_vs_growth_KDE.jpeg')

    X1 = X1.reshape(-1, 1)
    Y1 = Y1.reshape(-1, 1)

    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(X1, Y1)  # perform linear regression
    Y_pred = linear_regressor.predict(X1)  # make predictions

    fig3 = plt.scatter(X1, Y1)
    plt.plot(X1, Y_pred, color='red')
    f = fig3.get_figure()
    f.savefig(path_out + '/Graphs/genome_vs_growth_reg.jpeg')

    pcorr, _ = pearsonr(X1, Y1)
    scorr, _ = spearmanr(X1, Y1)

    pearson = 'Pearson: ', str(pcorr)
    spearman = 'Spearman:', str(scorr)
    coef = 'Coefficients:', str(linear_regressor.coef_)
    r2 = 'R2 Score:', str(linear_regressor.score(X1, Y1))
    intercept = 'Intercept:', str(linear_regressor.intercept_)
    f = open(path_out + 'regression_scores.txt', 'w')
    f.write("%s \n %s \n %s \n %s \n %s \n" % (pearson, spearman, coef, r2, intercept))
    f.close()

if regression_growth_rate_vs_genome:
    # This is saved in a general folder (not dimension specific) since the output will be the same for every dimension
    # The growth rate used in this case is the one a microbe would achieve under UDM with all its ingredients
    X2 = final_compilation_table_genome_growth.iloc[:, 1]
    y2 = final_compilation_table_growth.iloc[:, growth_column]
    g = sns.regplot(X2, y2)
    g.set(ylabel='growth rate (mmol/gDW/h)', xlabel='genome size (Mb)')
    fig4 = g.get_figure()
    fig4.savefig(path_out + 'Graphs/genome_vs_g_rate_reg.jpeg')

    t = sns.scatterplot(X2, y2)
    t.set(ylabel='growth rate (mmol/gDW/h)', xlabel='genome size (Mb)')
    fig5 = t.get_figure()
    fig5.savefig(path_out + 'Graphs/genome_vs_g_rate.jpeg')

    linear_regressor2 = LinearRegression()  # create object for the class

    pcorr, _ = pearsonr(X2, y2)
    scorr, _ = spearmanr(X2, y2)

    pearson2 = 'Pearson: ', str(pcorr)
    spearman2 = 'Spearman:', str(scorr)
    coef2 = 'Coefficients:', str(linear_regressor2.coef_)
    r22 = 'R2 Score:', str(linear_regressor2.score(X2, y2))
    intercept2 = 'Intercept:', str(linear_regressor2.intercept_)
    f = open(path_out + 'regression_scores.txt', 'w')
    f.write("%s \n %s \n %s \n %s \n %s \n" % (pearson2, spearman2, coef2, r22, intercept2))
    f.close()

# a compilation table with nutrient combination necessary for growth is generated from final_compilation_table_growth
# in other words, 1s (ones) found in final_compilation_table_growth are replaced by their corresponding index headers
#  NOTE THAT THERE MIGHT BE AN EASIER WAY TO DO THIS WITH LESS STEPS AND USING THE SAME final_compilation_table_growth
outcome_dict = {}

for index, row in enumerate(final_compilation_table_growth.iterrows()):
    microbe = row[0]
    outcome_list = []
    for item in range(len(final_compilation_table_growth.columns)):
        result = final_compilation_table_growth.iloc[index][final_compilation_table_growth.columns[item]]
        if result > 0:
            product = final_compilation_table_growth.columns[item]
            outcome_list.append(product)
        else:
            outcome_list.append('0')
    outcome_dict.update({microbe: outcome_list})

outcome_table = pd.DataFrame.from_dict(outcome_dict, orient='index')
outcome_table.columns = final_compilation_table_growth.columns
outcome_table.to_csv(path_out + 'compilation/growth_compilation_by_nutrient.csv')


