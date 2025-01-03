"""
Created by Juan M.
on 26/03/2021
"""

"""
This script was designed to test strains' metabolite degradation capabilities from a list of previously specified metabolites.

This does not assess growth by itself but growth and metabolite degradation.
"""

import cobra
import pandas as pd
import itertools
import os
from os import listdir
from os.path import isfile, join

# Paths for inputs and outputs to be specified by the user
path_out = 'C:/'
path_in = 'C:/'

# Name of the output folder containing Boolean table
output_folder = '/'

if not os.path.exists(path_out + output_folder):
    os.makedirs(path_out + output_folder)
    os.makedirs(path_out + output_folder + 'growth_no-growth/')
    os.makedirs(path_out + output_folder + 'compilation/')

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
               'valine': "EX_val_L(e)",
               'L-cystine': "EX_Lcystin(e)", 'L-methionine sulfoxide': "EX_metsox_S_L(e)",
               'carnitine': "EX_crn(e)", 'ornithine': "EX_orn(e)",
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
    'D-Galacturonate': 'EX_galur(e)'
}

experimental_carbon_sources = {
  # Proteoglycans
    'Core 2': {'Core 2': 'EX_core2(e)'},
    'Core 3': {'Core 3': 'EX_core3(e)'},
    'Core 4': {'Core 4': 'EX_core4(e)'},
    'Core 5': {'Core 5': 'EX_core5(e)'},
    'Core 6': {'Core 6': 'EX_core6(e)'},
    'Core 7': {'Core 7': 'EX_core7(e)'},
    'Core 8': {'Core 8': 'EX_core8(e)'},
    'Chondroitin Sulfate A Proteoglycan': {'Chondroitin Sulfate A Proteoglycan': 'EX_cspg_a(e)'},
    'Dimer GalNAc4S-GlcA': {'Dimer GalNAc4S-GlcA': 'EX_cspg_a_degr(e)'},
    'Chondroitin Sulfate B Proteoglycan': {'Chondroitin Sulfate B Proteoglycan': 'EX_cspg_b(e)'},
    'Dimer GalNAc4S-IdoA2S': {'Dimer GalNAc4S-IdoA2S': 'EX_cspg_b_degr(e)'},
    'Chondroitin Sulfate C Proteoglycan': {'Chondroitin Sulfate C Proteoglycan': 'EX_cspg_c(e)'},
    'Dimer GalNAc6S-GlcA': {'Dimer GalNAc6S-GlcA': 'EX_cspg_c_degr(e)'},
    'Disialyl-T Antigen': {'Disialyl-T Antigen': 'EX_dsT_antigen(e)'},
    'GlcNac-Alpha-1,4-Core 1': {'GlcNac-Alpha-1,4-Core 1': 'EX_gncore1(e)'},
    'GlcNac-Alpha-1,4-Core 2': {'GlcNac-Alpha-1,4-Core 2': 'EX_gncore2(e)'},
    'Hyaluronic acid': {'Hyaluronic acid': 'EX_ha(e)'},
    'Heparan Sulfate Proteoglycan': {'Heparan Sulfate Proteoglycan': 'EX_hspg(e)'},
    'Product 1 HEPARL1_e': {'Product 1 HEPARL1_e': 'EX_hspg_degr_1(e)'},
    'Product AS3TASE_HS2': {'Product AS3TASE_HS2': 'EX_hspg_degr_10(e)'},
    'Product GAM2STASE_HS3': {'Product GAM2STASE_HS3': 'EX_hspg_degr_11(e)'},
    'Product GLCNACASE_HS2': {'Product GLCNACASE_HS2': 'EX_hspg_degr_12(e)'},
    'Product GLCAASE_HS': {'Product GLCAASE_HS': 'EX_hspg_degr_13(e)'},
    'Product GLCNACASE_HS3': {'Product GLCNACASE_HS3': 'EX_hspg_degr_14(e)'},
    'Product IDOURASE_HS3': {'Product IDOURASE_HS3': 'EX_hspg_degr_15(e)'},
    'Product 2 HEPARL1_e': {'Product 2 HEPARL1_e': 'EX_hspg_degr_2(e)'},
    'Product IDOURASE_HS1': {'Product IDOURASE_HS1': 'EX_hspg_degr_3(e)'},
    'Product AS6TASE_HS1': {'Product AS6TASE_HS1': 'EX_hspg_degr_4(e)'},
    'Product AS3TASE_HS1': {'Product AS3TASE_HS1': 'EX_hspg_degr_5(e)'},
    'Product GAM2STASE_HS2': {'Product GAM2STASE_HS2': 'EX_hspg_degr_6(e)'},
    'Product IS2TASE_HS1': {'Product IS2TASE_HS1': 'EX_hspg_degr_7(e)'},
    'Product IDOURASE_HS2': {'Product IDOURASE_HS2': 'EX_hspg_degr_8(e)'},
    'Product AS6TASE_HS2': {'Product AS6TASE_HS2': 'EX_hspg_degr_9(e)'},
    'Sialyl-T Antigen': {'Sialyl-T Antigen': 'EX_sT_antigen(e)'},
    'Sialyl-Tn Antigen': {'Sialyl-Tn Antigen': 'EX_sTn_antigen(e)'},
    'T-Antigen (Core 1)': {'T-Antigen (Core 1)': 'EX_T_antigen(e)'},
    'Tn Antigen': {'Tn Antigen': 'EX_Tn_antigen(e)'},
  # Polysaccharides
    'Amylopectin': {'Amylopectin': 'EX_amylopect900(e)'},
    'Amylose': {'Amylose': 'EX_amylose300(e)'},
    'Arabinan': {'Arabinan': 'EX_arabinan101(e)'},
    'Larch arabinogalactan': {'Larch arabinogalactan': 'EX_arabinogal(e)'},
    'Arabinoxylan': {'Arabinoxylan': 'EX_arabinoxyl(e)'},
    'Beta-Glucans': {'Beta-Glucans': 'EX_bglc(e)'},
    'Cellulose': {'Cellulose': 'EX_cellul(e)'},
    'Dextran': {'Dextran': 'EX_dextran40(e)'},
    'Carob galactomannan': {'Carob galactomannan': 'EX_galmannan(e)'},
    'Glycogen, Structure 2': {'Glycogen, Structure 2': 'EX_glygn2(e)'},
    'Glycogen, Structure 4': {'Glycogen, Structure 4': 'EX_glygn4(e)'},
    'Glycogen, Structure 5': {'Glycogen, Structure 5': 'EX_glygn5(e)'},
    'Homogalacturonan': {'Homogalacturonan': 'EX_homogal(e)'},
    'Inulin': {'Inulin': 'EX_inulin(e)'},
    'Levan': {'Levan': 'EX_levan1000(e)'},
    'Lichenin': {'Lichenin': 'EX_lichn(e)'},
    'Laminarin': {'Laminarin': 'EX_lmn30(e)'},
    'Alpha-mannan, yeast': {'Alpha-mannan, yeast': 'EX_mannan(e)'},
    'Pectin': {'Pectin': 'EX_pect(e)'},
    'Pectic galactan (potato)': {'Pectic galactan (potato)': 'EX_pecticgal(e)'},
    'Pullulan': {'Pullulan': 'EX_pullulan1200(e)'},
    'Potato rhamnogalacturonan I': {'Potato rhamnogalacturonan I': 'EX_rhamnogalurI(e)'},
    'Wine rhamnogalacturonan II': {'Wine rhamnogalacturonan II': 'EX_rhamnogalurII(e)'},
    'Starch': {'Starch': 'EX_starch1200(e)'},
    'Starch, Structure 1': {'Starch, Structure 1': 'EX_strch1(e)'},
    'Xylan': {'Xylan': 'EX_xylan(e)'},
    'Xyluglucan': {'Xyluglucan': 'EX_xyluglc(e)'},
  # Oligosaccharides
    'Arabinotriose': {'Arabinotriose': 'EX_arabttr(e)'},
    'Cellobiose':{'Cellobiose': 'EX_cellb(e)'},
    'N,N-diacetylchitobiose': {'N,N-diacetylchitobiose': 'EX_chtbs(e)'},
    'Kestopentaose': {'Kestopentaose': 'EX_kestopt(e)'},
    'Kestotetraose': {'Kestotetraose': 'EX_kestottr(e)'},
    'D-Maltose': {'D-Maltose': 'EX_malt(e)'},
    'Maltohexaose': {'Maltohexaose': 'EX_malthx(e)'},
    'Mannotriose (beta-1,4)': {'Mannotriose (beta-1,4)': 'EX_mantr(e)'},
    'Melibiose': {'Melibiose': 'EX_melib(e)'},
    'Raffinose': {'Raffinose': 'EX_raffin(e)'},
    'Stachyose': {'Stachyose': 'EX_stys(e)'},
    'Sucrose': {'Sucrose': 'EX_sucr(e)'},
    'Trehalose': {'Trehalose': 'EX_tre(e)'},
    'Alpha-Lactose': {'Alpha-Lactose': 'EX_lcts(e)'},
  # Monosaccharides
    'N-Acetylgalactosamine': {'N-Acetylgalactosamine': 'EX_acgal(e)'},
    'N-Acetyl-D-glucosamine': {'N-Acetyl-D-glucosamine': 'EX_acgam(e)'},
    'N-Acetylneuraminic acid': {'N-Acetylneuraminic acid': 'EX_acnam(e)'},
    'D-Arabinose': {'D-Arabinose': 'EX_arab_D(e)'},
    'L-Arabinose': {'L-Arabinose': 'EX_arab_L(e)'},
    'Deoxyribose': {'Deoxyribose': 'EX_drib(e)'},
    'D-Fructose': {'D-Fructose': 'EX_fru(e)'},
    'L-Fucose': {'L-Fucose': 'EX_fuc_L(e)'},
    'Glucosamine': {'Glucosamine': 'EX_gam(e)'},
    'D-Glucose': {'D-Glucose': 'EX_glc_D(e)'},
    'L-lyxose': {'L-lyxose': 'EX_lyx_L(e)'},
    'D-Mannose': {'D-Mannose': 'EX_man(e)'},
    'D-Ribose': {'D-Ribose': 'EX_rib_D(e)'},
    'L-Rhamnose': {'L-Rhamnose': 'EX_rmn(e)'},
    'Salicin': {'Salicin': 'EX_salcn(e)'},
    'D-Xylose': {'D-Xylose': 'EX_xyl_D(e)'},
  # Intermediate metabolites
    'Acetic acid': {'Acetic acid': 'EX_ac(e)'},
    'Acetaldehyde': {'Acetaldehyde': 'EX_acald(e)'},
    'Formic acid': {'Formic acid': 'EX_for(e)'},
    'L-Lactic acid': {'L-Lactic acid': 'EX_lac_L(e)'},
    'Malic acid': {'Malic acid': 'EX_mal_L(e)'},
    'Propionate': {'Propionate': 'EX_ppa(e)'},
    'Pyruvic acid': {'Pyruvic acid': 'EX_pyr(e)'},
    'Butyrate': {'Butyrate': 'EX_but(e)'},
    'Succinate': {'Succinate': 'EX_succ(e)'},
    'Fumarate': {'Fumarate': 'EX_fum(e)'},
    'Ethanol': {'Ethanol': 'EX_etoh(e)'}
}

nutrient_qtt = 100

rich_media_no_exp_source = {}
rich_media_no_exp_source.update(amino_acids)
rich_media_no_exp_source.update(main_cofactors)
rich_media_no_exp_source.update(other)
rich_media_no_exp_source.update(bile_acids)
rich_media_no_exp_source.update(fatty_acids)
rich_media_no_exp_source.update(dipeptide)
rich_media_no_exp_source.update(secondary_cofactors)
rich_media_no_exp_source.update(metals)
rich_media_no_exp_source.update(anions)
rich_media_no_exp_source.update(cations)

# Creates a list of strain names (model names) located in the path_in directory when running several microbes at once
models_in = [f for f in listdir(path_in) if isfile(join(path_in, f))]
models_in = [os.path.splitext(f)[0] for f in models_in]

final_compilation_table = pd.DataFrame()

for file in models_in:
    print(file)

    # an SBML model is loaded to be simulated
    model = cobra.io.read_sbml_model(path_in + file + '.xml')  # Models with the same name as 'file' are loaded
    
    rich_media_no_exp_source_df = pd.DataFrame()

    final_boolean_table = pd.DataFrame()
    
    for ingredient in rich_media_no_exp_source:
        code = rich_media_no_exp_source[ingredient]
        new_ingredient = pd.DataFrame([nutrient_qtt], index=[code])
        rich_media_no_exp_source_df = pd.concat([rich_media_no_exp_source_df, new_ingredient])

    rich_media = rich_media_no_exp_source_df
    code_list = []

    for ingredient in experimental_carbon_sources:
        code_list.append(ingredient)

    for r in range(0, len(code_list) + 1):
        for subset in itertools.combinations(code_list, r):
            combination = list(subset)
           
            if len(combination) > 1:
                break

            rich_media = rich_media_no_exp_source_df
            grouped = {}
            for group in combination:
                grouped.update(experimental_carbon_sources[group])

            list_of_reactions = []

            for code in grouped:
                list_of_reactions.append(grouped[code])

            for ingredient in list_of_reactions:
                new_ingredient = pd.DataFrame([nutrient_qtt], index=[ingredient])
                rich_media = pd.concat([rich_media, new_ingredient])

            media_dict = rich_media.to_dict()
            uptakes = media_dict[0]

            with model:
                medium = model.medium

                for reaction in medium:
                    if reaction not in uptakes:
                        medium[reaction] = 0.0
                    else:
                        medium[reaction] = uptakes[reaction]

                model.medium = medium
                value = 0

                for added_reaction in list_of_reactions:
                    if added_reaction in model.reactions:
                        constraint = model.problem.Constraint(model.reactions.get_by_id(added_reaction).flux_expression,
                                                              lb=-1000, ub=-0.0001)
                        model.add_cons_vars(constraint)
                        solution = model.optimize()

                        if added_reaction in solution.fluxes and solution.fluxes[
                            added_reaction] != 0.0 and \
                                solution.objective_value is not None and solution.objective_value > 0.09:
                            value = 1
                             
                nutrient_test = pd.DataFrame([value], index=[str(subset)])
                
                nutrient_test.columns = [file]
                
                final_boolean_table = pd.concat([final_boolean_table, nutrient_test])
                
            solution = model.slim_optimize()

    
    final_boolean_table.to_csv(path_out + output_folder + 'growth_no-growth/' + file + '.csv')
