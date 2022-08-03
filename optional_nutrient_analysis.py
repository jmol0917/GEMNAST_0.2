"""
Created by Juan M.
on 26/03/2021
"""

"""
This script reads individual microbial GSMs from a file location. Microbes metabolism is then inspected to check if a 
nutrient that has been previously removed from the media has been synthesized intra-cellularly. The set of 
nutrients to be inspected has to be specified.
"""

import cobra
from cobra.exceptions import OptimizationError
import pandas as pd
import seaborn as sns
import os
import warnings
from os import listdir
from os.path import isfile, join

warnings.filterwarnings("error")

path_in = ''
path_out = ''

output_folder = ''

if not os.path.exists(path_out + output_folder):
    os.makedirs(path_out + output_folder)
    os.makedirs(path_out + output_folder + 'graphs/')
    os.makedirs(path_out + output_folder + 'compilation/')
    os.makedirs(path_out + output_folder + 'cluster/')

# Dictionary with experimental energy sources
simple_sugars = {
    # Example
    'D-glucose': "EX_glc_D(e)"
                 }

amino_acids = {
               }

cations = {
}

anions = {
}

metals = {
}

main_cofactors = {
}

secondary_cofactors = {
}

dipeptide = {
}

fatty_acids = {
}

bile_acids = {
}

other = {
}

# nutrients from every group are added to the media. Nutrients inspected for production are removed from the media
# below
rich_media_no_explored_n = {}
rich_media_no_explored_n.update(simple_sugars)
rich_media_no_explored_n.update(amino_acids)
rich_media_no_explored_n.update(main_cofactors)
rich_media_no_explored_n.update(other)
rich_media_no_explored_n.update(bile_acids)
rich_media_no_explored_n.update(fatty_acids)
rich_media_no_explored_n.update(dipeptide)
rich_media_no_explored_n.update(main_cofactors)
rich_media_no_explored_n.update(secondary_cofactors)
rich_media_no_explored_n.update(metals)
rich_media_no_explored_n.update(anions)
rich_media_no_explored_n.update(cations)

explored_groups = {

    # Example
    'B1': {'thiamine': "thm", 'thiamine monophosphate': "thmmp"},

}

# Creates a list of bacteria names (models) located in the path_in directory when running several microbes at once
models_in = [f for f in listdir(path_in) if isfile(join(path_in, f))]
models_in = [os.path.splitext(f)[0] for f in models_in]

# Testing individual models
# models_in = ['Enterococcus_faecalis_V583']

rich_media_df = pd.DataFrame()

for ingredient in rich_media_no_explored_n:
    code = rich_media_no_explored_n[ingredient]
    new_ingredient = pd.DataFrame([100], index=[code])
    rich_media_df = pd.concat([rich_media_df, new_ingredient])

production_boolean_table = pd.DataFrame()

for name in models_in:

    print(name)
    microbe_boolean_table = pd.DataFrame()

    for explored_group in explored_groups:

        # Most recent update 18/10/2021: I have moved the following line from line 227 to this location to load the
        # model every time a new metabolite is examined. This prevents crashes.
        model = cobra.io.read_sbml_model(path_in + name + '.xml')

        media_dict = rich_media_df.to_dict()
        uptakes = media_dict[0]

        group_of_reactions = explored_groups[explored_group]

        # remove reactions that belong to the same group from the media above
        for metabolite in group_of_reactions:
            reaction = group_of_reactions[metabolite]
            ex_reaction = 'EX_' + reaction + '(e)'
            if ex_reaction in uptakes:
                del uptakes[ex_reaction]

        # value is out of the lower loop, so if value changes for one of the reactions in the current group it
        # conserves a value of 1 even if the later reactions in the group don't return a positive outcome.
        value = 0

        for metabolite in group_of_reactions:
            reaction = group_of_reactions[metabolite]
            explored_metabolite = reaction + '[c]'
            with model:
                medium = model.medium

                for ingredient in medium:
                    if ingredient not in uptakes:
                        medium[ingredient] = 0.0

                model.medium = medium
                if explored_metabolite in model.metabolites:
                    try:
                        solution = model.optimize()
                        if solution.objective_value is not None and solution.objective_value > 0.09 and solution.status != 'Infeasible':
                            metabolite_reactions = model.metabolites.get_by_id(explored_metabolite).summary().to_string()
                            if 'Empty DataFrame' not in metabolite_reactions:
                                print(explored_metabolite, 'HAS BEEN used in one or more reactions')
                                value = 1
                    except (UserWarning, OptimizationError):
                        value = 0

        group_test = pd.DataFrame([value], index=[explored_group])
        group_test.columns = [name]
        microbe_boolean_table = pd.concat([microbe_boolean_table, group_test])

    microbe_boolean_table = microbe_boolean_table.transpose()
    production_boolean_table = pd.concat([production_boolean_table, microbe_boolean_table])

production_boolean_table.to_csv(path_out + output_folder + 'compilation/synthesis_compilation.csv')

with open(path_out + output_folder + 'experimental_design.txt', 'w') as file:
    file.write('This results were generated using the nutrient_biosynthesis_aa_vit.py script')
    line = '\n\n'
    file.write(str(rich_media_no_explored_n))
