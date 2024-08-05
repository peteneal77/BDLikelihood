README - R Covid Paper
----------------------

This is a Readme file for the R code used in Supplementary Material Section 8 of: "Fast likelihood calculations for emerging epidemics" Ball and Neal (2024).

Output Folder - Stores output from analysis especially related to the exact likelihood calculations.

BD_sim - Simulation birth-death process

Likelihood functions
BD_Exact_Like - Exact birth-death process likelihood
BD_Approx_Like - Approximate birth-death process likelihood 
BD_Hybrid_Like - Hybrid birth-death process likelihood where M (transition from exact to approximate likelihood calculations) can be chosen.
BD_Hybrid_Like_Fix - Hybrid birth-death process likelihood where M (transition from exact to approximate likelihood calculations) is fixed at M=50.

MCMC_BD - MCMC algorithm for obtaining samples from the posterior distribution of the parameters of a time-inhomogeneous birth-death process. (Parameters piecewise constant between event times.) General models for how parameters define birth and death rates and detection probabilities are allowed.

BD_para - Code for transforming model parameters into birth and death rates and detection probabilities.

BD_prior - Code for calculating prior densities.

Simultation_Study_BD - Study of contour plots of the likelihood surface for a time-homogeneous birth-death process using the exact, approximate and hybrid likelihoods.

Simulation_Scenario_BD - Simulation study of four simulated data sets from different time-inhomogeneous birth-death processes. MCMC output from the posterior distribution are stored in Output.

Note in Simultation_Study_BD and Simultation_Scenario_BD the set.seed() command is used to allow reproducible analysis of the reported simulation study. Altering the seed value will enable further robustness checks.

Early_Diff - Code used to compute possible differences between the exact and approximate likelihood for k=4 and k=5 in a time-homogeneous birth-death process.
Early_Diff_Sim - Code to simulate multiple realisations of the first 5 death times in a time-homogeneous birth-death process and study the difference between the exact and approximate likelihoods.
