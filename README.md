# tb-rdm-method
This repository contains code for paper on random directed graph model for TB household contact study. 

## Bayesian Implementation 
- `simulation.cpp`: code of MCMC implemention for the model without powerful predictors, i.e., model 2 in the paper. 
- `simulation_powerful.cpp`: code of MCMC implementation for the model with powerful predictors, i.e., model 1 in the paper. 
- `weight_comparison.cpp`: code of computing DIC for the purpose of model comparison (with different weighting schemes). 

## Data & Output
- `out1.rds`: DIC output for the first scenario, i.e., where extra-household transmission is predominant (accounts for about 90% of the total infections).
- `out2.rds`: DIC output for the second scenario, i.e., where extra-household transmission is stronger than household transmission (accounts for about 50%-60% of the total infections).
- `out3.rds`: DIC output for the third scenario, i.e., where household transmission is stronger than extra-household transmission (accounts for about 50%-60% of the total infections).
- `out4.rds`: DIC output for the fourth scenario, i.e., where household transmission is predominant (accounts for about 90% of the total infections).
- `output.R`: Proprocessing code.
- `prediction.rds`: Simulation results based on the Brazilian household contact study. Used as the final interpretable results given by the random directed graph model.
- `si_hhc.rds`: Posterior sample obtained based on the household contact study data only (without community controls). i = 1, 2, 3 or 4 signaling the scenario ID. 
- `si_hw.rds`:  95% credible intervals obtained based on the household contact study data under optimal weighting schemes. 
- `si_whole.rds`: Posterior sample obtained based on the whole data (with perfect community controls). 
- `simulation_prediction.R`: The R code used to generate `prediction.rds`. 
- `unweight_est.rds`: Posterior sample obtained based on the Brazilian household contact study without weights.
- `weighted_est.rds`: Posterior samples obtained based on the Brazilian household contact study under different weighting schemes.
- `weight_DIC.rds`: DIC for posterior samples in `weighted_est.rds`. 

## Simulation
- `powersim.R`: data simulation function for the model with powerful predictors, i.e., model 1 in the paper.
- `simulated_data.R`: data simulation function for the model without powerful predictors, i.e., model 2 in the paper. 
