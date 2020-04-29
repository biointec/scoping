# Parental_Selection_for_long_term_profits
## Running the code
- Install the package hypred (zip file included)
- Install the package GSSimTPUpdate (zip file included)
- Open the R project Article_files.Rproj
- Run Create_directory.R
- Run MakeGenome_File.R
- Run one of the 'run_experiment' files in the main directory
- Transfer the simulation results from the 'own_results' directory to the 'data' directory
- Run the correct 'make_figure' script

## data 
This directory contains the original base population used for each simulation study and several directories where the simulation results can be stored.

## Figures
This directory contains the original figures used for the article.

## make_figures
This directory contains the R scripts used to make the figures. The simulated data should be put in the data/\<method\> directory. 

## own_results
Empty directory where simulation results will be saved. This directory is creating by running the script Create_directory.R.

## MakeGenome_File 
Makes a list of n different genomes that can be used in the run_experiment files. At each iteration, each method will use the same genome, making it possible to compare the different methods.

## Genome
Directory containing an example files of genomes that can be used in the run_experiment files.

## run_experiment files
R scripts used to simulate a population following the oracle, baseline, backcrossing, scoping or combined selection method. The results will be saved in the directory 'own_results'. 
The user will have to load the correct genome from the genome directory and choose the number of nodes that are available for calculation. 

## Supplementary data
Directory containing the supplementary figures and tables.

## R
Directory containing functions that are used in the run_experiment files.