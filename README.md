Use of common control arms across sub-studies in a master protocol via ridge estimation
====

Program codes for each estimator used in simulation and application are available. 

# Description of each file

- methods_SBR.R: Function to compute each estimator, which was used for the simulation study
- hbm_sum_stats.stan: Stan code to compute hierarchical Bayesian model, which was used for the simulation study
- application_CPMP.R: Program codes to generate the application results
- hbm_sum_stats_case_study.stan: Stan code to compute hierarchical Bayesian model, which was used for the application
- CPMP_study_result_asof_16Jun2024.csv: Aggregated dataset of CPMP study results which was disclosed in clinicaltrials.gov as of 16Jun2024
