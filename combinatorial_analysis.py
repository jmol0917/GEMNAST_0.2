"""
Created by Juan M.
on 01/07/2020
"""

"""
This script was designed to test every possible combination of a set of nutrients that are added to a 
Universally Defined Media (UDM).
The components of both, UDM media and the list of nutrients to be tested, should be specified as
dictionaries below. 

It has been expanded to work with other nutrient dimensions; minor adjustments are needed.

Follow comments to run this either on a personal computer or on the HPC. For this script to run on the HPC, the 
count_files.sh document will be needed. After uploading count_files to the HPC it has to be converted to unix with the
'dos2unix count_files.sh' command. Submit a job as any other regular job using the './count_files.sh' command.
"""

import cobra
from cobra.exceptions import OptimizationError
import pandas as pd
import itertools
import os
from os import listdir
from os.path import isfile, join
import warnings

# Path where simulation results are saved
# path_out = 'results/AGORA/'
# path_in = 'sbml/'

# Uncomment the following two lines to read the name of a model from the command line
# file = str(sys.argv[1])
# file = os.path.splitext(file)[0]

warnings.filterwarnings("error")

# Local path options
path_out = 'C:/Users/jpmo_/Dropbox (Sydney Uni)/gsm_analysis/results/AGORA/'
path_in = 'C:/Users/jpmo_/Dropbox (Sydney Uni)/gsm_analysis/sbml/'

output_folder = 'requirements_no_ubiquinone/'

if not os.path.exists(path_out + output_folder):
    os.makedirs(path_out + output_folder)
    os.makedirs(path_out + output_folder + 'consumption/growth_no-growth/')
    os.makedirs(path_out + output_folder + 'consumption/growth_rate/')
    os.makedirs(path_out + output_folder + 'consumption/growth_instances/')

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
               'L-cystine': "EX_Lcystin(e)", 'L-methionine sulfoxide': "EX_metsox_S_L(e)",
               'carnitine': "EX_crn(e)", 'ornithine': "EX_orn(e)"
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
    'biotin': "EX_btn(e)", 'menaquionine-7': "EX_mqn7(e)",
    'cobalamin I': "EX_cbl1(e)", 'menaquionine-8': "EX_mqn8(e)",
    'cobalamin II': "EX_cbl2(e)", 'nicotinic acid': "EX_nac(e)",
    'adenosylcobalamin': "EX_adpcbl(e)",
    'folic acid': "EX_fol(e)", 'niacinamide': "EX_ncam(e)",
    'nicotinamide ribotide': "EX_nmn(e)", 'pantothenic acid': "EX_pnto_R(e)",
    'pyridoxine': "EX_pydxn(e)", 'reduced riboflavin': "EX_rbflvrd(e)",
    'riboflavin': "EX_ribflv(e)", 'tetrahydrofolic acid': "EX_thf(e)",
    'thiamine': "EX_thm(e)", 'thiamine monophosphate': "EX_thmmp(e)",
    'demethylmenaquinone': "EX_2dmmq8(e)", 'pyridoxal': "EX_pydx(e)",
    'pyridoxamine': "EX_pydam(e)", 'ubiquinone-8': "EX_q8(e)"
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
    'Glycyl-L-methionine': 'EX_glymet(e)', 'Spermidine': 'EX_spmd(e)',
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
    'Glycerol 3-phosphate': 'EX_glyc3p(e)', '4-Aminobenzoate': 'EX_4abz(e)',
    'Glutathione': 'EX_gthrd(e)', 'Diaminoheptanedioate': 'EX_26dap_M(e)',
    'Dephospho-CoA': 'EX_dpcoa(e)', '1,2-Diacyl-sn-glycerol': 'EX_12dgr180(e)',
    'Methyl-Oxovaleric Acid': 'EX_3mop(e)', 'Chorismate': 'EX_chor(e)',
    '4-Hydroxybenzoic acid': 'EX_4hbz(e)', 'Oxidized glutathione': 'EX_gthox(e)',
    'Putrescine': 'EX_ptrc(e)', 'Indole': 'EX_indole(e)',
    'Lanosterin': 'EX_lanost(e)', 'Choline sulfate': 'EX_chols(e)',
    'Ketobutyric acid': 'EX_2obut(e)', 'Glycolaldehyde': 'EX_gcald(e)',
    'Trimethylamine': 'EX_tma(e)', 'NADP': 'EX_nadp(e)',
    'Acetic acid': 'EX_ac(e)', 'Formic acid': 'EX_for(e)',
    'Gamma-butyrobetaine': 'EX_gbbtn(e)', 'Acetoacetic acid': 'EX_acac(e)',
    'Coenzyme A': 'EX_coa(e)', 'Ethanolamine': 'EX_etha(e)',
    'Tetrathionate': 'EX_tet(e)', 'Dehydro-deoxy-gluconate': 'EX_2ddglcn(e)',
    'Carbon dioxide': 'EX_co2(e)', 'Allantoin': 'EX_alltn(e)',
    'Cholesterol': 'EX_chsterol(e)', 'Formaldehyde': 'EX_fald(e)',
    'Water': 'EX_h2o(e)', 'Phenylpyruvic acid': 'EX_phpyr(e)',
    'Urea': 'EX_urea(e)', 'L-Lactic acid': 'EX_lac_L(e)',
    'D-Galacturonate': 'EX_galur(e)',
    'Citric acid': 'EX_cit(e)', 'Glycerol': 'EX_glyc(e)',
    'Malic acid': 'EX_mal_L(e)',
    'Acetylmannosamine': 'EX_acmana(e)',
    'ubiquinone-8': "EX_q8(e)",
    'biotin': "EX_btn(e)"
}

main_cofactors_grouped = {

    'B1': {'thiamine': "EX_thm(e)", 'thiamine monophosphate': "EX_thmmp(e)"},
    'B2': {'riboflavin': "EX_ribflv(e)", 'reduced riboflavin': "EX_rbflvrd(e)"},
    'B3': {'nicotinic acid': "EX_nac(e)", 'niacinamide': "EX_ncam(e)", 'nicotinamide ribotide': "EX_nmn(e)"},
    'B5': {'pantothenic acid': "EX_pnto_R(e)"},
    'B6': {'pyridoxine': "EX_pydxn(e)", 'pyridoxal': "EX_pydx(e)", 'pyridoxamine': "EX_pydam(e)"},
    # 'B7': {'biotin': "EX_btn(e)"},
    'B9': {'folic acid': "EX_fol(e)", 'tetrahydrofolic acid': "EX_thf(e)"},
    'B12': {'cobalamin I': "EX_cbl1(e)", 'cobalamin II': "EX_cbl2(e)", 'adenosylcobalamin': "EX_adocbl(e)"},
    'K': {'menaquionine-7': "EX_mqn7(e)", 'menaquionine-8': "EX_mqn8(e)", 'demethylmenaquinone': "EX_2dmmq8(e)"}
          # 'ubiquinone-8': "EX_q8(e)"}

}


nutrient_uptake = 20

rich_media_no_exp_source = {}
rich_media_no_exp_source.update(amino_acids)
rich_media_no_exp_source.update(simple_sugars)
rich_media_no_exp_source.update(other)
rich_media_no_exp_source.update(bile_acids)
rich_media_no_exp_source.update(fatty_acids)
rich_media_no_exp_source.update(dipeptide)
rich_media_no_exp_source.update(secondary_cofactors)
rich_media_no_exp_source.update(metals)
rich_media_no_exp_source.update(anions)
rich_media_no_exp_source.update(cations)

# Creates a list of bacteria names (models) located in the path_in directory when running several microbes at once
models_in = [f for f in listdir(path_in) if isfile(join(path_in, f))]
models_in = [os.path.splitext(f)[0] for f in models_in]

# Uncomment this line to test a single microbe
# models_in = ['Lactobacillus_fermentum_ATCC_14931']

# models_in2 = [f for f in listdir(path_out + 'consumption/') if isfile(join(path_out + 'consumption/', f))]
# models_in2 = [os.path.splitext(f)[0] for f in models_in2]

final_compilation_table = pd.DataFrame()

for file in models_in:
    print(file)

    # an SBML model is loaded to be simulated
    model = cobra.io.read_sbml_model(path_in + file + '.xml')  # Models with the same name as 'file' are loaded

    genome_size = 0

    # genome size (if available) is recorded for later comparison
    # GENOME SIZE WASN'T DETERMINED FOR TWO MICROBES IN AGORA
    every_genome_size = pd.read_csv('../microbe_and_genome.csv', index_col=[0])
    if every_genome_size.loc[file, 'Genome size'] > 0.0:
        genome_size = every_genome_size.loc[file, 'Genome size']

    count = 0

    rich_media_no_exp_source_df = pd.DataFrame()

    final_boolean_table = pd.DataFrame()
    growth_rate_table = pd.DataFrame()

    for ingredient in rich_media_no_exp_source:
        code = rich_media_no_exp_source[ingredient]
        new_ingredient = pd.DataFrame([nutrient_uptake], index=[code])
        rich_media_no_exp_source_df = pd.concat([rich_media_no_exp_source_df, new_ingredient])

    rich_media = rich_media_no_exp_source_df
    code_list = []

    for ingredient in main_cofactors_grouped:
        code_list.append(ingredient)

    for r in range(0, len(code_list) + 1):
        for subset in itertools.combinations(code_list, r):
            combination = list(subset)
            print(combination)

            rich_media = rich_media_no_exp_source_df
            grouped = {}

            for group in combination:
                grouped.update(main_cofactors_grouped[group])

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

    #     Uncomment the line below to write one file per model to the directory specified on path_out
    final_boolean_table.to_csv(path_out + output_folder + 'consumption/growth_no-growth/' + file + '.csv')

    # Since the file generated in the following lines will be compared to genome size, there is no reason to create it
    # if the genome size of the microbe in question is not known
    if genome_size > 0.0:
        growth_rate_table.to_csv(path_out + output_folder + 'consumption/growth_rate/' + file + '_growth_rate.csv')
    #
    # The following lines of code register genome size and number of time a microbe grew
    genome_growth_info = pd.DataFrame([count, genome_size], index=['growth instances', 'genome size'])
    genome_growth_info.columns = [file]
    genome_growth_info.to_csv(path_out + output_folder +
                              'consumption/growth_instances/' + file + '_genome_vs_growth.csv')
