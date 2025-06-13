Seasonal Duckweed microbiomes experiment data and analysis

This repository contains the data and analyses associated with manuscript:
Seasonal differences in effects of microbial communities on duckweed traits
Authored by: Emma Kinerson, Alex Trott, and Anna M. O'Brien

All data and analysis scripts are attributable to the same authors.

Script ProcessImageJOutputData.r uses the output files generated from manual use of image J from experimental images:
NewDuckweedPlateMap.xlsx
NewDuckweedPlateDat.xlsx
It also uses information about the images in file:
ImageSettings.xlsx 
The output file from this script is summaries of the duckweed data in wells:
DuckweedData.csv

Script SeasonAnalysis.r is the primary data analysis script.
This script uses the output of the previous script as input, as well as treatment and optical density information:
DuckweedTreatments.csv
OpticalDensityData.xlsx
This script outputs a concatenated file of all treatments and data:
MappedGrowth_trt_data.csv
This script also outputs all remaining files, which are main text or supplementary figures (all .png files)

