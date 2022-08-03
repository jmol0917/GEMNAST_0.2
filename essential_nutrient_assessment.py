"""
Created by Juan M.
on 01/07/2020
"""

"""
This script was designed to test every possible combination of a set of nutrients that are added to a 
Universally Defined Media (UDM).
The components of both, UDM media and the list of nutrients to be tested, should be specified as
dictionaries below. 

It has been expanded to work with other nutrient dimensions.
"""

import cobra
from cobra.exceptions import OptimizationError
import pandas as pd
import itertools
import os
from os import listdir
from os.path import isfile, join
import warnings

warnings.filterwarnings("error")

# Local path options
path_out = ''
path_in = ''

output_folder = ''

if not os.path.exists(path_out + output_folder):
    os.makedirs(path_out + output_folder)
    os.makedirs(path_out + output_folder + 'consumption/growth_no-growth/')

simple_sugars = {
    # example
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

explored_nutrients = {
    # example
    'B1': {'thiamine': "EX_thm(e)", 'thiamine monophosphate': "EX_thmmp(e)"},
}


nutrient_uptake = 20

rich_media_no_exp_source = {}
rich_media_no_exp_source.update(amino_acids)
rich_media_no_exp_source.update(simple_sugars)
rich_media_no_exp_source.update(other)
rich_media_no_exp_source.update(bile_acids)
rich_media_no_exp_source.update(fatty_acids)
rich_media_no_exp_source.update(dipeptide)
rich_media_no_exp_source.update(main_cofactors)
rich_media_no_exp_source.update(secondary_cofactors)
rich_media_no_exp_source.update(main_cofactors)
rich_media_no_exp_source.update(metals)
rich_media_no_exp_source.update(anions)
rich_media_no_exp_source.update(cations)

# Creates a list of bacteria names (models) located in the path_in directory when running several microbes at once
models_in = [f for f in listdir(path_in) if isfile(join(path_in, f))]
models_in = [os.path.splitext(f)[0] for f in models_in]

final_compilation_table = pd.DataFrame()

for file in models_in:
    print(file)

    # an SBML model is loaded to be simulated
    model = cobra.io.read_sbml_model(path_in + file + '.xml')  # Models with the same name as 'file' are loaded
    
    rich_media_no_exp_source_df = pd.DataFrame()

    final_boolean_table = pd.DataFrame()
    growth_rate_table = pd.DataFrame()

    for ingredient in rich_media_no_exp_source:
        code = rich_media_no_exp_source[ingredient]
        new_ingredient = pd.DataFrame([nutrient_uptake], index=[code])
        rich_media_no_exp_source_df = pd.concat([rich_media_no_exp_source_df, new_ingredient])

    rich_media = rich_media_no_exp_source_df
    code_list = []

    for ingredient in explored_nutrients:
        code_list.append(ingredient)

    for r in range(0, len(code_list) + 1):
        for subset in itertools.combinations(code_list, r):
            combination = list(subset)

            rich_media = rich_media_no_exp_source_df
            grouped = {}

            for group in combination:
                grouped.update(explored_nutrients[group])

            list_of_reactions = []

            for code in grouped:
                list_of_reactions.append(grouped[code])

            for ingredient in list_of_reactions:
                new_ingredient = pd.DataFrame([nutrient_uptake], index=[ingredient])
                rich_media = pd.concat([rich_media, new_ingredient])

            media_dict = rich_media.to_dict()
            uptakes = media_dict[0]
            growth_rate = 0.0
            value = 0
            with model:
                medium = model.medium

                for reaction in medium:
                    if reaction not in uptakes:
                        medium[reaction] = 0.0
                    else:
                        medium[reaction] = uptakes[reaction]

                model.medium = medium
                try:
                    solution = model.optimize()
                    if solution.objective_value < 0.09:
                        print(solution.objective_value)
                        value = 0
                    else:
                        if solution.objective_value is not None and solution.objective_value > 0.09:
                            growth_rate = solution.objective_value
                            count += 1
                            print('Growth')
                            value = 1
                            for reaction in grouped:
                                if grouped[reaction] in solution.fluxes:
                                    if solution.fluxes[grouped[reaction]] < 0.0:
                                        print(growth_rate)
                                        print(reaction, 'was consumed at the following rate:',
                                              solution.fluxes[grouped[reaction]])
                except (UserWarning, OptimizationError):
                    value = 0

                nutrient_test = pd.DataFrame([value], index=[str(subset)])
                growth_rate_log = pd.DataFrame([growth_rate], index=[str(subset)])

                nutrient_test.columns = [file]
                growth_rate_log.columns = [file]

                final_boolean_table = pd.concat([final_boolean_table, nutrient_test])
                growth_rate_table = pd.concat([growth_rate_table, growth_rate_log])

            solution = model.optimize()

    final_boolean_table.to_csv(path_out + output_folder + 'consumption/growth_no-growth/' + file + '.csv')

with open(path_out + output_folder + 'experimental_design.txt', 'w') as file:
    file.write('This results were generated using the essential_nutrient_assessment.py script\n\n')
    file.write('Using the following media:\n', str(rich_media_no_exp_source))
