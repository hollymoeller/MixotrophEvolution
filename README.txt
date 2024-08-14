Archibald, K. M., S. Dutkiewicz, C. Laufkötter, and H. V. Moeller (submitted). Evolved reductions in respiration in marine mixotrophic microbes under thermal stress. The American Naturalist.

Code Author: Kevin Archibald, karchibald@ucsb.edu

This directory contains MATLAB code and data needed to reproduce the model simulations contained in the manuscript. This is a mechanistic model of mixotroph growth and respiration based on a set of differential equations that describe short-term, plastic changes to metabolic strategy and long-term, evolutionary changes to respiration under variable temperature conditions. The model was used to simulate changes to metabolic activity at different temperatures under evolutionary and non-evolutionary scenarios.

Directory Contents:

mixo_model.m -- MATLAB script containing model equations
params.m -- MATLAB script containing model parameter values
simulation.m -- MATLAB script to run model simulations used in the manuscript
figures.m -- MATLAB script to produce models contained in the manuscript

fig2_data.mat -- MATLAB data file containing simulation output for reciprocal transplant simulations
	T -- nx1 vector of temperature values (degrees C)
	Tevo -- mx1 vector of evolutionary temperatures for mixotrophs (degrees C)
	dn -- nxm matrix of the following model run diagnostics
		qc -- quality control flag (1 = simulation converged to constant growth rate, 0 = no convergence)
		ngr -- net growth rate (per day)
		strat -- metabolic strategy (growth, photosynthesis, and grazing investment)
		cellc -- carbon content per cell (pg C)
		celln -- nitrogen content per cell (pg N)
		photo -- photosynthetic rate (pg C per cell per day)
		graz -- grazing rate (prey per cell per day)
		resp -- respiration rate (per day)
		cue -- carbon use efficiency
fig3_data*.mat -- MATLAB data file containing simulation output for long-term time-series simulation. Separate files for each simulation (temperature increase and decrease, control and evolving populations).
	Tinc or Tdec -- time-series of temperature values (degrees C)
	dn* -- data structure containing the time-series diagnostics (same as above)

figs1.mat -- MATLAB data file containing sensitivity analysis output data.
	T -- nx1 vector of temperature (degrees C)
	u -- nxn matrix of growth rates using baseline parameters (columns = evolutionary temperature, rows = environmental temperature)
	uminus -- nxnx22 matrix of growth rates using parameters - 20%
	uplus -- nxnx22 matrix of growth rates using parameters + 20%

Experimental Data -- directory containing experimental data files compiled from Lepori-Bui et al. (2022).
	cue_highlight.csv -- text file containing carbon use efficiency data with the following data columns
		Strain -- ID number for mixotroph strain
		Ev -- temperature each strain was evolved at (degrees C)
		Acc -- temperature each strain was acclimated at (degrees C)
		N -- number of replicates for assay
		CUE -- carbon use efficiency measurement 
		sd -- standard deviation
		se -- standard error
		ci -- width of confidence interval
	graze_highlight.csv -- grazing rate measurements (prey per cell per day)
	growth_highlight.csv -- growth rate measurements (per day)
	photo_highlight.csv -- photosynthesis measurements (gC per gC per day)
	resp_highlight.csv -- respiration measurements (gC per gC per day)

All analysis conducted in MATLAB: '9.9.0.1538559 (R2020b) Update 3'

