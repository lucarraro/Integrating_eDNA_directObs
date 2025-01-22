# Integrating_eDNA_directObs

Scripts supporting Carraro, L., "Integrating direct observation and environmental DNA data to enhance species distribution models in riverine environments"

 - `MAIN.R` executes the simulations for the in-silico experiment.
 - `analyze_data.R` analyzes the simulation output and produces the manuscript figures.
 - `eval_wilcox_paired.R` stores utility functions performing statistical tests and other graphical operations.
 - `test_multiplicity.R` produces Fig. S1 of the manuscript.
 - `OCN` contains the OCNs used in the in-silico experiment.
 - `results` contains simulation output in compact form (`results_df.rda`) and can store the output of `MAIN.R` (one file per simulation when `MAIN.R` is run). It also contains simulation output for the multiplicity test both in compact form (`results_multiplicity.rda`) and in extensive form (one file per simulation within subfolders `multiplicity_BT` and `multiplicity_optim`).
 

