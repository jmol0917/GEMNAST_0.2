"""
Created by Juan M.
on 26/03/2021
"""

"""
Script to analyse intermediate metabolite export and monosaccharide degradation \.
This is an independent analysis that extracts data directly from GSM growth modelling without reading data from other 
sources.
"""

import cobra
from cobra.exceptions import Infeasible
import pandas as pd
import seaborn as sns
import os
from os import listdir
from os.path import isfile, join

path_in = 'C:/'
path_out = 'C:/'

"""
Create necessary paths if they don't exist already
"""

output_folder = 'im_export/'

if not os.path.exists(path_out + output_folder):
  os.makedirs(path_out + output_folder)
  os.makedirs(path_out + output_folder + 'compilation/')

# Modified UDM by categories
simple_sugars = {'D-glucose': "EX_glc_D(e)", 'fructose': "EX_fru(e)",
                 'galactose': "EX_gal(e)", 'mannose': "EX_man(e)",
                 'ribose': "EX_rib_D(e)", 'lactose': "EX_lcts(e)",
                 'L-fucose': "EX_fuc_L(e)", 'inulin': "EX_inulin(e)",
                 'maltose': "EX_malt(e)", 'D-xylose': "EX_xyl_D(e)",
                 'sucrose': "EX_sucr(e)", 'arabinose': "EX_arab_D(e)",
                 'ac_glucosamine': "EX_acgam(e)", 'chitobiose': "EX_chtbs(e)",
                 'glucosamine': "EX_gam(e)", 'D-Galacturonate': 'EX_galur(e)'
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
    '4-Aminobenzoate': 'EX_4abz(e)',
    'Glutathione': 'EX_gthrd(e)', 'Diaminoheptanedioate': 'EX_26dap_M(e)',
    'Dephospho-CoA': 'EX_dpcoa(e)', '1,2-Diacyl-sn-glycerol': 'EX_12dgr180(e)',
    'Methyl-Oxovaleric Acid': 'EX_3mop(e)',
    'Chorismate': 'EX_chor(e)',
    '4-Hydroxybenzoic acid': 'EX_4hbz(e)', 'Oxidized glutathione': 'EX_gthox(e)',
    'Putrescine': 'EX_ptrc(e)', 'Indole': 'EX_indole(e)',
    'Lanosterin': 'EX_lanost(e)',
    'Choline sulfate': 'EX_chols(e)',
    'Trimethylamine': 'EX_tma(e)', 'NADP': 'EX_nadp(e)',
    'Gamma-butyrobetaine': 'EX_gbbtn(e)',
    'Ethanolamine': 'EX_etha(e)',
    'Tetrathionate': 'EX_tet(e)',
    'Dehydro-deoxy-gluconate': 'EX_2ddglcn(e)',
    'Carbon dioxide': 'EX_co2(e)', 'Allantoin': 'EX_alltn(e)',
    'Cholesterol': 'EX_chsterol(e)', 'Formaldehyde': 'EX_fald(e)',
    'Water': 'EX_h2o(e)',
    'Phenylpyruvic acid': 'EX_phpyr(e)',
    'Urea': 'EX_urea(e)',
    }

rich_media_no_vit_k = {}
# rich_media_no_vit_k.update(simple_sugars)  Simple sugars are not part of modified UDM
rich_media_no_vit_k.update(amino_acids)
rich_media_no_vit_k.update(main_cofactors)
rich_media_no_vit_k.update(other)
rich_media_no_vit_k.update(bile_acids)
rich_media_no_vit_k.update(fatty_acids)
rich_media_no_vit_k.update(dipeptide)
rich_media_no_vit_k.update(secondary_cofactors)
rich_media_no_vit_k.update(metals)
rich_media_no_vit_k.update(anions)
rich_media_no_vit_k.update(cations)

explored_ch = {
    
    'N-Acetylgalactosamine': 'EX_acgal(e)',
    'N-Acetyl-D-glucosamine': 'EX_acgam(e)',
    'N-Acetylneuraminic acid': 'EX_acnam(e)',
    'D-Arabinose': 'EX_arab_D(e)',
    'L-Arabinose': 'EX_arab_L(e)',
    'Deoxyribose': 'EX_drib(e)',
    'D-Fructose': 'EX_fru(e)',
    'L-Fucose': 'EX_fuc_L(e)',
    'Glucosamine': 'EX_gam(e)',
    'D-Glucose': 'EX_glc_D(e)',
    'L-lyxose': 'EX_lyx_L(e)',
    'D-Mannose': 'EX_man(e)',
    'D-Ribose': 'EX_rib_D(e)',
    'L-Rhamnose': 'EX_rmn(e)',
    'Salicin': 'EX_salcn(e)',
    'D-Xylose': 'EX_xyl_D(e)'
}

intermediate_metabolites = {
    'Acetic acid': 'EX_ac(e)',
    'Acetaldehyde': 'EX_acald(e)',
    'Formic acid': 'EX_for(e)',
    'L-Lactic acid': 'EX_lac_L(e)',
    'Malic acid': 'EX_mal_L(e)',
    'Propionate': 'EX_ppa(e)',
    'Pyruvic acid': 'EX_pyr(e)',
    'Butyrate': 'EX_but(e)',
    'Methylbutyrate': 'EX_2mbut(e)',
    'Succinate': 'EX_succ(e)',
    'Fumarate': 'EX_fum(e)'
}

# Creates a list of bacteria names (models) located in the path_in directory when running several microbes at once
models_in = [f for f in listdir(path_in) if isfile(join(path_in, f))]
models_in = [os.path.splitext(f)[0] for f in models_in]

rich_media_df = pd.DataFrame()

for ingredient in rich_media_no_vit_k:
    code = rich_media_no_vit_k[ingredient]
    new_ingredient = pd.DataFrame([100], index=[code])
    rich_media_df = pd.concat([rich_media_df, new_ingredient])

production_boolean_table = pd.DataFrame()

for name in models_in:
    metabolites_generated_pre_exploration = []
    microbe_boolean_table = pd.DataFrame()
    model = cobra.io.read_sbml_model(path_in + name + '.xml')
    print(name)

    for ch in explored_ch:
        media_dict = rich_media_df.to_dict()
        uptakes = media_dict[0]
        ch_reaction = explored_ch[ch]
        uptakes[ch_reaction] = 1000

        value = 0

        with model:
            medium = model.medium

            for ingredient in medium:
                if ingredient not in uptakes:
                    medium[ingredient] = 0.0

            model.medium = medium
            if ch_reaction in model.reactions:
                for r in intermediate_metabolites:
                    metabolite_reaction = intermediate_metabolites[r]
                    if metabolite_reaction in model.reactions:
                    #     molecule constraint to be degraded
                        constraint = model.problem.Constraint(model.reactions.get_by_id(ch_reaction).flux_expression,
                                                              lb=-1000, ub=-0.0001)
                        model.add_cons_vars(constraint)
                    #     molecule constraint to be exported
                        constraint = model.problem.Constraint(model.reactions.get_by_id(metabolite_reaction).
                                                              flux_expression, lb=0.0001, ub=1000)
                        model.add_cons_vars(constraint)
                        try:
                            solution = model.optimize()
                            if solution.fluxes[ch_reaction] < 0.0 and solution.fluxes[metabolite_reaction] > 0.0 and \
                                    solution.objective_value is not None and solution.objective_value > 0.09:
                                value = 1
                               
                            else:
                                value = 0
                        except (UserWarning, Infeasible):
                            value = 0
                            print('error')
                    else:
                        value = 0
                    production_test = pd.DataFrame([int(value)], index=[ch + ', ' + r])
                    production_test.columns = [name]
                    microbe_boolean_table = pd.concat([microbe_boolean_table, production_test])

            else:
                value = 0
                for r in intermediate_metabolites:
                    production_test = pd.DataFrame([int(value)], index=[ch + ', ' + r])
                    production_test.columns = [name]
                    microbe_boolean_table = pd.concat([microbe_boolean_table, production_test])

    microbe_boolean_table = microbe_boolean_table.transpose()
    production_boolean_table = pd.concat([production_boolean_table, microbe_boolean_table])

production_boolean_table.to_csv(path_out + output_folder + 'compilation/production_compilation.csv')

