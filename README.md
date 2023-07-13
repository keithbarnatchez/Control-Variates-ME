## README.md

This file describes the code layout of the control variates project

#### MAIN CODE FOR RUNNING SIMULATIONS
`ms-exposure-sim.R`: This file is used for carrying out the main simulation 
exercise as described in the manuscript. Data-generating parameters to be considered
can be specified in the file. This simulation code calls functions from the 
following files:

- `ms-exposure-funcs.R`
This file contains code for 1) generating data according to the DGP in the 
manuscript, 2) implementing the C.V. and MIME estimators (as well as the naive
and oracle estimators), and 3) implementing a fixed iteration of the simulation
procedure

- `generalizability_funcs.R`
This file includes functions for implementing the **generalizability** version of 
the control variates estimator (e.g. the estimator that takes into account 
systematic differences in the covariate distribution between the validation data
and overall data)

#### CLUSTER CODE
- `ms-causal-exp.sh`
This file submits a batch job request to run `ms-exposure-sim.R`. Edit it if you'd
like to vary the amount of memory/cores used (the simulation is carried out in
parallel)

#### OUTPUT FILES
 - The `output` folder contains the results of different simulation runs, where the 
filename corresponds to the date/time of the simulation run. Each row
in the output dataset corresponds to a **single iteration** of a simulation (I
don't aggregate upfront so that we can calculate different stats of interest 
that may arise down-stream), with different operating characteristics

- `plotting_res.R` is the main file used for plotting results. The user specifies
a file from the `output` folder. Then, the results are averaged within different
combinations of parameter values, after which plots of %bias, RMSE and CI coverage
are produced

#### DATA EXAMPLE
The file `clean_data.R` is used for cleaning the raw NHANES data files stored in 
the `data` folder

#### OLD FILES
The `archive` folder contains old code from previous iterations of the project. 
None of the files stored here are currently being used


