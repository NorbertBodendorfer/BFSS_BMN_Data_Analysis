This repository contains a selection of the data analysis scripts that I used for analysing data form Monte Carlo Simulations of the BFSS and BMN matrix models. The results are published in the following papers:
https://arxiv.org/abs/1909.04592
https://arxiv.org/abs/2005.04103
https://arxiv.org/abs/2110.01312
https://arxiv.org/abs/2205.06098
https://arxiv.org/abs/2210.04881

The scripts are written in R and either take raw output data from the MC simulations (samples found in prod_broad_output) or preprocessed data (samples found in processed). Executing the scripts prompts instructions on how to use them correctly. (Preprocessed) data for the examples provided is contained in the repository. 



Filenames of the raw output data read e.g. outGN24S12T0.543M2.0D9_F1_5.txt
This means:
out: it is an output file (as opposed to e.g. configurations files containing complete field configurations, which are not included in this repository)
G: it is the gauged matrix model (as opposed to ungauged = UG)
N24: number of colours is 24
S12: number of temporal lattice points is 12
T0.543: the dimensionless temperature is 0.543
M2.0: the BMN flux parameter is 2.0
D9: the number of matrices is 9
F1_5: this is the 5th output file in the Monte Carlo stream 1 for the above parameters. Multiple output files per stream are needed due to runtime limits on the HPC clusters. Multiple streams are useful to increase parallelization or for varying other parameters not included in the file name. 



There are two types of preprocessed files included in the repository:

matrix_sizes_GN8S24T0.3M0.5D9_F11_.csv
- contains the average sizes and errors of the 9 squared matrices plus their sum for the given stream number and other parameters

observables_GN16S72T0.3M0.5D9_F11_.csv
- contains the average sizes and errors of a set of observables for the given stream number and other parameters



The R scripts contain the path to the R executable at the beginning. This may need to be adapted depending on where they are run. 

There are sample outputs of the scripts in the folder example_output


List of scripts:

Script:        r_binning_joint.R 
What is does:  Plot histograms of an observable for various temperatures in the same plot. This is useful for visualising phase transitions, see e.g. https://arxiv.org/abs/2110.01312
Usage:         ./r_binning_joint.R 1000 prod_broad_output GN24S12 9 2.0 0.542 0.543 0.544
Output sample: plot_jointbins_Polyak_GN24S12M2.0D9acl1.pdf



Script:        r_correlogram.R
What is does:  Plots a correlogram of several observables. Uses base R. 
Usage:         ./r_correlogram.R prod_broad_output GN24S12T0.543M2.0D9_F1 1 20000 
Output sample: correlogram_GN24S12T0.543M2.0D9_F1.pdf

Script:        r_correlogram_gg.R
What is does:  Plots a correlogram of several observables. Uses base ggplot2 and produces nicer output and smaller plots than base R. 
Usage:         ./r_correlogram_gg.R prod_broad_output GN24S12T0.543M2.0D9_F1 1 20000 
Output sample: correlogram_gg_GN24S12T0.543M2.0D9_F1.pdf

Script:        r_matrix_asymmetry.R
What is does:  Processes a MC stream and output a preprocessed file containing the average sizes of the nine matrices squares and their sum. 
Usage:         ./r_matrix_asymmetry.R prod_broad_output GN16S30T0.3M0.5D9_F12 2000 -1 
Output sample: see processed/matrix...

Script:        r_matrix_asymmetry_large_S.R
What is does:  Performs a large S extrapolation of matrix sizes to check whether size differences (caused by lattice artefacts) vanish in the continuum limit, see e.g. https://arxiv.org/abs/2210.04881
Usage:         ./r_matrix_asymmetry_large_S.R GN16 T0.3M0.5D9_F
Output sample: large_S_GN16T0.3M0.5D9_F.pdf

Script:        r_MC_history_single.R
What is does:  Plots the MC history of a single observable using base R. 
Usage:         ./r_MC_history_single.R prod_broad_output GN24S12T0.543M2.0D9_F1 1 80000 
Output sample: MChistory_GN24S12T0.543M2.0D9_F1_Polyak.pdf

Script:        r_observable_extractor.R
What is does:  Processes a single MC stream (or multiple if stream number is not given) and calculates averages and errors on a list of observables specified within the script. 
Usage:         ./r_observable_extractor.R prod_broad_output GN16S30T0.3M0.5D9_F12_ 2000 -1 
Output sample: see processed/observables...

Script:        r_observable_large_N.R
What is does:  Performs a large N extrapolation for given simulation parameters. 
Usage:         ./r_observable_large_N.R G S30T0.3M0.5D9_F 
Output sample: observable_large_N_GS30T0.3M0.5D9_F.pdf observable_large_N_GS30T0.3M0.5D9_F_errors.pdf

Script:        r_observable_large_S.R
What is does:  Performs a continuum extrapolation for given simulation parameters. 
Usage:         ./r_observable_large_S.R GN16 T0.3M0.5D9_F
Output sample: observable_large_S_GN16T0.3M0.5D9_F.pdf observable_large_S_GN16T0.3M0.5D9_F_errors.pdf

Script:        r_observable_large_NS.R 
What is does:  Performs a simultaneous large N and continuum extrapolation for the given simulation parameters. Performs a Kolmogorov-Smirnov test to check whether the model fits the data well enough (errors have to be standard normally distributed). 
Usage:         ./r_observable_large_NS.R G T0.3M0.5D9_F 
Output sample: observable_large_NS_GT0.3M0.5D9_F.pdf observable_large_NS_KS_test_GT0.3M0.5D9_F.pdf

Script:        r_PD_1d_bins.R 
What is does:  Fits a power law to the dependence of the average of one observable to another. Useful to establish partial deconfinement in https://arxiv.org/abs/1909.04592 
Usage:         ./r_PD_1d_bins.R prod_broad_output GN24S12T0.543M2.0D9_F1 1 20000 0 0.5 0.7
Output sample: plot_1d_bins_Polyak_s_trx2_GN24S12T0.543M2.0D9_F1.pdf

Script:        r_PD_2d_bins.R
What is does:  Makes a 2d density plot of two observables. Useful to establish partial deconfinement in https://arxiv.org/abs/1909.04592 
Usage:         ./r_PD_2d_bins.R prod_broad_output GN24S12T0.543M2.0D9_F1 1 20000 
Output sample: plot_2d_bins_Polyak_energy_GN24S12T0.543M2.0D9_F1.pdf

Script:        
What is does:  
Usage:         
Output sample: 


There are also some additional R scripts that contain functions called by the above scripts:
r_GeneralFunctions.R
r_GeneralSettings.R
r_ReadOutFilesToDataTable.R
r_ErrorBars.R