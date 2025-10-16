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
and overall data). We generally recommend the generalizability implementation, as 
the resulting estimator will tend to be more efficient 

Additional extension exercises which consider (i) more complex validation schemes, and (ii) model misspecification can be found in the `extensions` folder.


#### OUTPUT FILES
 - The `sim_results` folder contains the results of different simulation runs, where the 
filename corresponds to the date/time of the simulation run. Each row
in the output dataset corresponds to a **single simulation iteration**.

- `plotting-control-variates.R` is the main file used for plotting results. The user specifies
a file from the `sim_results` folder. Then, the results are averaged within different
combinations of parameter values, after which plots of %bias, RMSE and CI coverage
are produced


