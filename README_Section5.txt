README - R Covid Paper
----------------------

This is a Readme file for the analysis of Covid-19 in Europe in: "Fast likelihood calculations for emerging epidemics" Ball and Neal (2024), Section 5.

Simulation Folder
- Contains code to simulate daily counts of deaths (simulation_P_daily) 
- simulate_data - code to simulate a Covid-19 example with multiple countries and control measures
- SimControl2 - Dates 5 control measures are introduced within 5 countries (Simulation Example)
- DailyDataS2 - Daily counts of observed deaths in 5 countries (Simulation Example)

Datasets Folder
- Contains datasets for Covid-19 Worldwide and Europe along with simulated dataset
- Description file included

MCMC Output 
- Stores output of MCMC runs along with a description file

Main Folder
Setup and MCMC
- approxBD_functions : likelihood functions for approximating time-inhomogeneous Birth-Death process using Gaussian.
- parameter_calc : code for calculating (birth, death, detection) parameters during the course of the epidemic 
- daily_convert: code for generating uniformly distributed death times on a day given aggregated daily count data and for updating a randomly selected proportion of death times
- Covid_MCMC (Covid_MCMC2) 
  - Random walk Metropolis algorithm for sampling from the posterior distribution of the parameters 
  - Each iteration updates death times and model parameters
  - Initial covariance matrix (sigma * I - Covid_MCMC; general - Covid_MCMC2) 
  - Multiple runs of different length allowed - covariance matrix updated after each run and only samples from the final MCMC run are stored. 
  - Parameters can be fixed as can death times.

Rt 
- Rt_Calc : calculates how Rt varies over time
- Rt_Plot : plots Rt quantiles over time

European Covid-19 data
- Euro_Death.csv - daily counts (death) in 11 European countries
- Intervention_March.csv - Dates in March 2020 of introduction of NPIs in 11 European countries
- EuroCovidPrep : Prepares the data
- EuroAnalysis : Runs the MCMC code for the European Covid-19 data
- Covid_R0_Calc : Calculates Rt over time for the Covid-19 data

Simulation example
- DailyDataS2.csv - daily counts (removals) in 5 countries 
- SimControl2.csv - days on which control measures were implemented. 
- SimAnalysis : Runs the MCMC code for the Simulated data
- Sim_R0_Calc : Calculates Rt over time for the Simulated data and allows for comparison with the true Rt value

