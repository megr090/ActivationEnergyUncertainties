# ActivationEnergyUncertainties

This repository includes code associated with the preprint "The effect of uncertainties in creep activation energies on modeling ice flow and deformation", which computes
the dominant deformation mechanisms as a function of deformation rate and activation energies.

The main files here begin with "main". 

main_computeGamma_Idealized.m computes the fraction of dislocation creep in a probabilistic setup (using ensembles) in an idealized setup, where you can prescribe strain 
rate, ice thickness, surface mass balance, and compute a posterior ensemble of fraction of dislocation creep. This replicates Figure 2 of the preprint.

main_ComputeGamma_Antarctica.m computes the fraction of dislocation creep in a deterministic setup over the Antarctic Ice Sheet (only including fast-flowing regions, with
velocity > 30 m/yr). This code requires running a few codes first:
  main_ComputeIceTemperature_Antarctica.m estimates ice temperature in these regions of the Antarctic Ice Sheet following Meyer and Minchew 2018 (EPSL).
  main_computeDeformationMap.m computes the fraction of dislocation creep for varying strain rates and ice temperatures
Then, running this file replicates Figure 5 of the preprint.

The remaining figures can be recreated with the available code (please feel free to reach out at meghanar@ucar.edu if you would like the specific code for the 
remaining figures).
