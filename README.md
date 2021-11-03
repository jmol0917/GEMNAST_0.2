# gsm_analysis
The scripts in this repository have been designed to assess three aspects of microbial metabolism separately from corresponding Genome Scale Models (GSMs). The aspects that
are assessed are a) nutrient requirements (combinatorial_analysis.py), b) nutrient biosynthesis (nutrient_biosynthesis.py) and c) nutrient secretion (nutrient_secretion.py).
The scripts are compatible with the AGORA and vmh.life reaction and annotation format. For more information visit www.vmh.life

All of the scripts require the same setup:
1. Provide a directory in the 'path_out' variable. Output files will be saved here.
2. Provide a directory in the 'path_in' variable from where GSMs will be accessed.
3. Provide a name for the output folder (which will be automatically created in the specified directory in 'path_out').
4. Provide a list of nutrient names in at least one of the specified dictionaries. NOTE: There is no need to distribute nutrients by category (sugars, amino acids, cations, etc.) 
   and new categories (with a dictionary format) can be created and later added to the rich_media variable
5. Provide a list of nutrients to be explored within the "explored_groups" dictionary following the format provided in the example. NOTE that combinatorial_analysis and
   nutrient_secretion only work with exchange reactions (reactions with the format EX_xxx[e] where xxx should be replaced with the compound's abbreviation). Providing the
   compound's abbreviation alone is required for the nutrient_biosynthesis script to work.  

Different scripts provide different outputs.
1. Nutrient requirement analysis generates individual files per analysed GSM and the microbe compilation script is required to compile the information from such files into a
   single boolean table. Users need to provide a read/destination directory where individual files to be read are saved and where a Boolean compilation table will be saved 
   (path_out variable in line 26).
2. Nutrient biosynthesis and secretion both provide the same outputs: a Boolean table containing information from all the assessed GSMs is directly generated from both scripts. 
   A second 'reordered' boolean table and clustermap is also generated based on strain profile similarity. Finally, two files suggesting potential clusters in the data analysed 
   are also generated in its own folder.
   
For questions or comments please email juan.molina@sydney.edu.au
