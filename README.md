# GEMNAST
Computational high throughput platform for the assessment of nutritional and metabolic attributes of microorganisms/cells as represented in genome scale models. 
GEMNAST was written in the Python coding language (https://www.python.org/downloads/) and requires a compatible interface or IDE to be ran, additional to the instalation of the Python libraries specified within each of the corresponding scripts (e.g. COBRA: https://github.com/opencobra/cobrapy/blob/stable/INSTALL.rst). Installation of Anaconda Distribution is recommended: https://www.anaconda.com/products/distribution

The scripts are ready to be ran once the correct directory inputs and outputs are specified. Further, download of the corresponding GEM models in SBML format is required and the directory where these files are located should be provided as an input.

To download AGORA SBML files:
1. Go to https://www.vmh.life/#home
2. Click on 'Gut Microbiota'
3. Click on 'Download' at the bar near the TOP of the webpage
4. Click on 'Microbe collections'
5. Click on the 'AGORA 1.03 without mucin' version
6. All the models within this version will be downloaded as part of a .zip file

GEMNAST is divided in three folders. The Essential, Optional and Replaceable Nutrients analyses.

The Essential Nutrient analysis:
# Objective: Infer essential nutrients in a set of GEMs from a defined set of molecules
The essential_nutrients_assessment.py script produces individual files for every examined GEM. The essential_nutrients_compilation.py script can be used 
to generate a single boolean table with information from all the original files.

The optional_nutrients_assessment.py and the optional_nutrients_export.py files already produce a single boolean table with information from all the assessed models and
therefore additional scripts are not required.

The Replaceable Nutrients Assessment is composed of three scripts:
The replaceable_nutrients_degradation_assessment.py script produces individual files for every examined GEM. The replaceable_nutrients_degradation_compilation.py
script can be used to generate a single boolean table that compiles individual model results.

The replaceable_nutrients_degradation_and_intermediate_metabolite_synthesis.py and the replaceable_nutrients_degradation_and_intermediate_metabolite_export.py scripts produce boolean tables with information from every model assessed.

For questions or comments please email
