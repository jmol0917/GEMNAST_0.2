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

path_in = 'C:/Users/jpmo_/Dropbox (Sydney Uni)/gsm_analysis/sbml/'
path_out = 'C:/Users/jpmo_/Dropbox (Sydney Uni)/gsm_analysis/results/AGORA/'

output_folder = 'no_bio_b1b6_test/cofactor_synthesis/'

if not os.path.exists(path_out + output_folder):
    os.makedirs(path_out + output_folder)
    os.makedirs(path_out + output_folder + 'graphs/')
    os.makedirs(path_out + output_folder + 'compilation/')
    os.makedirs(path_out + output_folder + 'cluster/')

# Dictionary with experimental energy sources
simple_sugars = {'D-glucose': "EX_glc_D(e)", 'fructose': "EX_fru(e)",
                 'galactose': "EX_gal(e)", 'mannose': "EX_man(e)",
                 'ribose': "EX_rib_D(e)", 'lactose': "EX_lcts(e)",
                 'L-fucose': "EX_fuc_L(e)", 'inulin': "EX_inulin(e)",
                 'maltose': "EX_malt(e)", 'D-xylose': "EX_xyl_D(e)",
                 'sucrose': "EX_sucr(e)", 'arabinose': "EX_arab_D(e)",
                 'ac_glucosamine': "EX_acgam(e)", 'chitobiose': "EX_chtbs(e)",
                 'glucosamine': "EX_gam(e)"
                 }

amino_acids = {'D-alanine': "EX_ala_D(e)", 'alanine': "EX_ala_L(e)",
               'asparagine': "EX_asn_L(e)", 'aspartate': "EX_asp_L(e)",
               'arginine': "EX_arg_L(e)", 'cysteine': "EX_cys_L(e)",
               'glutamine': "EX_gln_L(e)", 'glycine': "EX_gly(e)",
               'glutamate': "EX_glu_L(e)", 'histidine': "EX_his_L(e)",
               'isoleucine': "EX_ile_L(e)", 'leucine': "EX_leu_L(e)",
               'lysine': "EX_lys_L(e)", 'D-methionine': "EX_met_D(e)",
               'L-methionine': "EX_met_L(e)", 'phenylalanine': "EX_phe_L(e)",
               'proline': "EX_pro_L(e)", 'D-serine': "EX_ser_D(e)",
               'L-serine': "EX_ser_L(e)", 'threonine': "EX_thr_L(e)",
               'tryptophan': "EX_trp_L(e)", 'tyrosine': "EX_tyr_L(e)",
               'valine': "EX_val_L(e)", 'aspactic acid': "EX_asp_L(e)",
               'L-cystine': "EX_Lcystin(e)", 'L-methionine sulfoxide': "EX_metsox_S_L(e)"
               }

cations = {
    'calcium': "EX_ca2(e)", 'cadmium': "EX_cd2(e)",
    'mercury': "EX_hg2(e)", 'magnesium': "EX_mg2(e)",
    'sodium': "EX_na1(e)", 'ammonia': "EX_nh4(e)",
    'potassium': "EX_k(e)", 'hydrogen ion': "EX_h(e)",
    'nitrogen': "EX_n2(e)"
}

anions = {
    'chloride ion': "EX_cl(e)", 'phosphate': "EX_pi(e)", 'sulfate': "EX_so4(e)",
    'sulfite': "EX_so3(e)", 'hydrogen sulfide': "EX_h2s(e)", 'hydrogen': "EX_h2(e)",
    'thiosulfate': "EX_tsul(e)", 'nitrite': "EX_no2(e)", 'nitrate': "EX_no3(e)",
}

metals = {
    'copper': "EX_cu2(e)", 'fe2': "EX_fe2(e)", 'cobalt': "EX_cobalt2(e)",
    'fe3': "EX_fe3(e)", 'manganese': "EX_mn2(e)", 'nickel': "EX_ni2(e)",
    'zinc': "EX_zn2(e)"
}

main_cofactors = {
    'biotin': "EX_btn(e)",
    'menaquionine-7': "EX_mqn7(e)",
    'cobalamin I': "EX_cbl1(e)",
    'menaquionine-8': "EX_mqn8(e)",
    'cobalamin II': "EX_cbl2(e)", 'nicotinic acid': "EX_nac(e)",
    'adenosylcobalamin': "EX_adpcbl(e)",
    'folic acid': "EX_fol(e)", 'niacinamide': "EX_ncam(e)",
    'nicotinamide ribotide': "EX_nmn(e)", 'pantothenic acid': "EX_pnto_R(e)",
    'pyridoxine': "EX_pydxn(e)",
    'reduced riboflavin': "EX_rbflvrd(e)",
    'riboflavin': "EX_ribflv(e)",
    'tetrahydrofolic acid': "EX_thf(e)",
    'thiamine': "EX_thm(e)", 'thiamine monophosphate': "EX_thmmp(e)",
    'demethylmenaquinone': "EX_2dmmq8(e)",
    'pyridoxal': "EX_pydx(e)",
    'pyridoxamine': "EX_pydam(e)",
    'ubiquinone-8': "EX_q8(e)"
}

secondary_cofactors = {
    'heme': "EX_pheme(e)", 'siroheme': "EX_sheme(e)",
    'thymidine': "EX_thymd(e)", 'cytosine': "EX_csn(e)",
    'uracil': "EX_ura(e)", 'adenosine': "EX_adn(e)",
    'adenine': "EX_ade(e)", 'guanine': "EX_gua(e)",
    'deoxyadenosine': "EX_dad_2(e)", 'deoxyguanosine': "EX_dgsn(e)",
    'guanosine': "EX_gsn(e)", 'guanosine triphosphate': "EX_gtp(e)",
    'Methylthioadenosine': "EX_5mta(e)", 'adenosine monophosphate': "EX_amp(e)",
    'S-adenosylmethionine': "EX_amet(e)", 'deoxyadenosine triphosphate': "EX_datp(e)",
    '5-Thymidylic acid': "EX_dtmp(e)", 'hypoxanthine': "EX_hxan(e)",
    'cytidine': "EX_cytd(e)", 'inosine': "EX_ins(e)",
    'xanthine': "EX_xan(e)", 'deoxycytidine': "EX_dcyt(e)",
    'uridine': "EX_uri(e)", 'deoxyinosine': "EX_din(e)",
    'cytidine monophosphate': "EX_cmp(e)", 'xanthosine': "EX_xtsn(e)"
}

dipeptide = {
    'Alanyl-glutamine': 'EX_alagln(e)', 'Carnosine': 'EX_alahis(e)',
    'Cysteinylglycine': 'EX_cgly(e)', 'Glycyl-L-asparagine': 'EX_glyasn(e)',
    'Glycyl-L-glutamine': 'EX_glygln(e)', 'Glycylleucine': 'EX_glyleu(e)',
    'Glycyl-L-methionine': 'EX_glymet(e)',
    'Gly-Cys': 'EX_glycys(e)', 'Glycyl-L-tyrosine': 'EX_glytyr(e)',
    'Glycyl-Phenylalanine': 'EX_glyphe(e)', 'L-alanyl-L-threonine': 'EX_alathr(e)',
    'L-methionyl-L-alanine': 'EX_metala(e)', 'L-alanyl-L-leucine': 'EX_alaleu(e)',
    'Glycylproline': 'EX_glypro(e)', 'L-alanyl-L-aspartate': 'EX_alaasp(e)',
    'L-alanylglycine': 'EX_alagly(e)', 'Alanyl-glutamate': 'EX_alaglu(e)',
    'Glycyl-L-aspartate': 'EX_glyasp(e)', 'Glycyl-L-glutamate': 'EX_glyglu(e)'
}

fatty_acids = {
    'Stearic acid': 'EX_ocdca(e)', 'Myristic acid': 'EX_ttdca(e)',
    'Dodecanoic acid': 'EX_ddca(e)', 'Oleic acid': 'EX_ocdcea(e)'
}

bile_acids = {
    'Chenodeoxycholic acid-glycine': 'EX_dgchol(e)',
    'Glycocholic acid': 'EX_gchola(e)',
    'Taurocholic acid': 'EX_tchola(e)'
}

other = {
    'Glycerol 3-phosphate': 'EX_glyc3p(e)',
    '4-Aminobenzoate': 'EX_4abz(e)',
    'Glutathione': 'EX_gthrd(e)', 'Diaminoheptanedioate': 'EX_26dap_M(e)',
    'Dephospho-CoA': 'EX_dpcoa(e)', '1,2-Diacyl-sn-glycerol': 'EX_12dgr180(e)',
    'Methyl-Oxovaleric Acid': 'EX_3mop(e)',
    'Chorismate': 'EX_chor(e)',
    '4-Hydroxybenzoic acid': 'EX_4hbz(e)', 'Oxidized glutathione': 'EX_gthox(e)',
    'Putrescine': 'EX_ptrc(e)', 'Indole': 'EX_indole(e)',
    'Lanosterin': 'EX_lanost(e)',
    'Choline sulfate': 'EX_chols(e)',
    'Ketobutyric acid': 'EX_2obut(e)',
    'Glycolaldehyde': 'EX_gcald(e)',
    'Trimethylamine': 'EX_tma(e)', 'NADP': 'EX_nadp(e)',
    'Acetic acid': 'EX_ac(e)',
    'Formic acid': 'EX_for(e)',
    'Gamma-butyrobetaine': 'EX_gbbtn(e)',
    'Acetoacetic acid': 'EX_acac(e)',
    'Coenzyme A': 'EX_coa(e)',
    'Ethanolamine': 'EX_etha(e)',
    'Tetrathionate': 'EX_tet(e)',
    'Dehydro-deoxy-gluconate': 'EX_2ddglcn(e)',
    'Carbon dioxide': 'EX_co2(e)', 'Allantoin': 'EX_alltn(e)',
    'Cholesterol': 'EX_chsterol(e)', 'Formaldehyde': 'EX_fald(e)',
    'Water': 'EX_h2o(e)',
    'Phenylpyruvic acid': 'EX_phpyr(e)',
    'Urea': 'EX_urea(e)',
    'L-Lactic acid': 'EX_lac_L(e)',
    'D-Galacturonate': 'EX_galur(e)',
    'Citric acid': 'EX_cit(e)',
    'Malic acid': 'EX_mal_L(e)',
    'Acetylmannosamine': 'EX_acmana(e)', 'Glycerol': 'EX_glyc(e)',
    'Carnitine': "EX_crn(e)", 'Ornithine': "EX_orn(e)", 'Spermidine': 'EX_spmd(e)'
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
rich_media_no_explored_n.update(secondary_cofactors)
rich_media_no_explored_n.update(metals)
rich_media_no_explored_n.update(anions)
rich_media_no_explored_n.update(cations)

explored_groups = {


    'B1': {'thiamine': "EX_thm(e)", 'thiamine monophosphate': "EX_thmmp(e)", 'thiamine pyrophosphate': "EX_thmpp[c]",
           'Hydroxy ethyl Thiamin diphosphate': "EX_hethmpp[c]"},
    'B2': {'riboflavin': "EX_ribflv(e)", 'reduced riboflavin': "EX_rbflvrd(e)", 'FAD': "EX_fad[c]", 'FMN': 'EX_fmn[c]'},
    'B3': {'nicotinic acid': "EX_nac(e)", 'niacinamide': "EX_ncam(e)", 'nicotinamide ribotide': "EX_nmn(e)",
           'NAD': "EX_nad[c]"},
    'B5': {'pantothenic acid': "EX_pnto_R(e)"},
    # Here I don't use 'pyridoxamine 5-phosphate': "EX_pyam5p[c]" beacuse it is not the active form and
    # it cannot be imported if released to the media to be converted into the active form. The rest of metabolites
    # can be converted into pydx5p[c]
    'B6': {'pyridoxine': "EX_pydxn(e)", 'pyridoxal': "EX_pydx(e)", 'pyridoxamine': "EX_pydam(e)",
           'pyridoxal 5-phosphate': "EX_pydx5p[c]"},
    # 'B7': {'biotin': "EX_btn(e)"},
    'B9': {'folic acid': "EX_fol(e)", 'tetrahydrofolic acid': "EX_thf(e)", '5-MTHF': "EX_5mthf[c]"},
    'B12': {'cobalamin I': "EX_cbl1(e)", 'cobalamin II': "EX_cbl2(e)", 'adenosylcobalamin': "EX_adocbl(e)",
            'adenosyl cobinamide phosphate': "EX_adocbip[c]", 'adenosyl cobinamide': "EX_adocbip[c]",
            'cobinamide': "EX_cbi[c]"},
    'K': {'menaquionine-7': "EX_mqn7(e)", 'menaquionine-8': "EX_mqn8(e)", 'demethylmenaquinone': "EX_2dmmq8(e)",
          'ubiquinone': "EX_q8(e)"}

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
            if reaction in uptakes:
                del uptakes[reaction]

        # value is out of the lower loop, so if value changes for one of the reactions in the current group it
        # conserves a value of 1 even if the later reactions in the group don't return a positive outcome.
        value = 0

        for metabolite in group_of_reactions:
            reaction = group_of_reactions[metabolite]
            explored_metabolite = reaction[3:-3] + '[c]'
            with model:
                medium = model.medium

                for ingredient in medium:
                    if ingredient not in uptakes:
                        medium[ingredient] = 0.0

                model.medium = medium
                if explored_metabolite in model.metabolites:
                    if reaction in model.reactions:
                        # print(model, reaction)
                        constraint = model.problem.Constraint(
                            model.reactions.get_by_id(reaction).flux_expression, lb=0.001, ub=50)
                        model.add_cons_vars(constraint)
                    # the try/catch loop shouldn't be nested in the if statement, there are some microbes that can
                    # synthesise a vitamin but don't have the export reaction that I use to constrain them above
                    try:
                        solution = model.optimize()
                        if solution.objective_value is not None and solution.objective_value > 0.09 and \
                                solution.status != 'Infeasible':
                            metabolite_reactions = model.metabolites.get_by_id(explored_metabolite)\
                                .summary().to_string()
                            if 'Empty DataFrame' not in metabolite_reactions:
                                print(explored_metabolite, 'HAS BEEN used in one or more reactions')
                                # print(metabolite_reactions)
                                value = 1
                    except (UserWarning, OptimizationError):
                        value = 0
                # else:
                #     print(explored_metabolite, "CANNOT be processed by this microbe.")

        group_test = pd.DataFrame([value], index=[explored_group])
        group_test.columns = [name]
        microbe_boolean_table = pd.concat([microbe_boolean_table, group_test])

    microbe_boolean_table = microbe_boolean_table.transpose()
    production_boolean_table = pd.concat([production_boolean_table, microbe_boolean_table])

production_boolean_table.to_csv(path_out + output_folder + 'compilation/synthesis_compilation.csv')

g = sns.clustermap(production_boolean_table, cmap="YlGnBu", cbar_pos=None, col_cluster=False, figsize=(75, 200),
                   dendrogram_ratio=(.1, .2))
g.savefig(path_out + output_folder + 'graphs/synthesis_compilation.jpeg')

reordered = g.data2d
reordered.to_csv(path_out + output_folder + 'compilation/synthesis_compilation_reordered.csv')

outcome_dict = {}

for index, row in enumerate(reordered.iterrows()):
    microbe = row[0]
    outcome_list = []
    for number in range(len(reordered.columns)):
        result = reordered.iloc[index][reordered.columns[number]]
        if result > 0:
            product = reordered.columns[number]
            outcome_list.append(product)
        else:
            outcome_list.append(0)
    outcome_dict.update({microbe: outcome_list})

outcome_table = pd.DataFrame.from_dict(outcome_dict, orient='index')
outcome_table.columns = reordered.columns
outcome_table.to_csv(
    path_out + output_folder + 'compilation/synthesis_compilation_reordered_nutrient_name.csv')

final_production_ruleset = {}

for index, row in outcome_table.iterrows():
    produced_nutrients = []
    for value in row:
        if value != 0 and value != '0':
            produced_nutrients.append(value)
    nutrient_set = ''
    for nutrient in produced_nutrients:
        nutrient_set = nutrient_set + nutrient + ','
    nutrient_set = '(' + nutrient_set[:-1] + ')'
    if nutrient_set not in final_production_ruleset:
        final_production_ruleset.update({nutrient_set: [index]})
    else:
        final_production_ruleset[nutrient_set].append(index)

with open(path_out + output_folder + 'cluster/production_by_nutrient_set.txt', 'w') as file:
    for combination in final_production_ruleset:
        inner_cluster = final_production_ruleset[combination]
        line = combination + ': ' + str(inner_cluster) + '\n\n'
        file.write(line)

with open(path_out + output_folder + 'cluster/production_by_nutrient_set_len.txt', 'w') as file:
    for combination in final_production_ruleset:
        inner_cluster = final_production_ruleset[combination]
        line = combination + ': ' + str(len(inner_cluster)) + '\n\n'
        file.write(line)
        line = '\n\n'

with open(path_out + output_folder + 'cluster/_experimental_design.txt', 'w') as file:
    file.write('This results were generated using the nutrient_biosynthesis_aa_vit.py script')
    line = '\n\n'
    file.write(str(rich_media_no_explored_n))
