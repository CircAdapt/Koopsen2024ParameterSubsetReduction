# Koopsen2024ParameterSubsetReduction

This folder contains Matlab files to perform Morris Screening Method (MSM), Sobol simulations, diaphony calculation, and Dynamic Multi-Swarm Particle Swarm Optimization (DMS-PSO).

The required CircAdapt P-structs are provided:
PRefHF.mat --> reference heart failure P-struct
PRef_AHA18.mat --> reference P-struct with 18 LV segments (according to AHA segmentation)
PRef_AHA18_CPP.mat --> reference P-struct with 18 LV segments (according to AHA segmentation) for use in C++ model solver

Performance of MSM requires the following files:
MorrisScreening_input.m
MorrisScreening_analysis.m
CalculateOutputs_Morris.m

Performance of Sobol simulations requires the following files:
SobolSamples.m
RunSobolSimulation.m

Performance of diaphony calculation (using the best Sobol simulations as determined from calculation of the objective function value) requires the following files:
ParSobolObjectiveFunction.m
SobolObjectiveFunction.m
CalculateObjectiveFunction.m
CalculateDiaphony.m
CalculateMaxDiaphony.m

Performance of DMS-PSO requires the following files:
DMSPSO_main.m
CalculateObjectiveFunction_18.m

Further instructions are provided in the respective Matlab files
