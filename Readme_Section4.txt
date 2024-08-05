README - R Covid Paper
----------------------

This is a Readme file for the analysis of Abakiliki data in: "Fast likelihood calculations for emerging epidemics" Ball and Neal (2024), Section 4.

Analysis folder - Contains pre-calculated results and figures. More details provided in separate readme file.

MCMC folder - MCMC output from MCMC runs detailed MCMC_Run_Abakiliki.R and MCMC_Run_Complete_Abakiliki.R

Calculation of log-likelihood
- EGSEI_functions.R - Exact GSE log-likelihood function calculation
- GSEI_functions.R - Birth-death process log-likelihood function calculations
- aGSEI_functions.R - approximate Birth-death process log-likelihood function calculation (Gaussian approximation)

Grid_paper.R - Code for evaluating the likelihood at grid points with T=90 and T=46.99 which feature in the paper. 
Grid_paper_V2.R - Additional code to Grid_paper.R for colour contour plots.

MCMC_IncGSE.R - Random walk Metropolis algorithm for epidemic in progress. Allows for choice of likelihood function (exact GSE, birth-death, approximate birth-death) 

MCMC_daily_IncGSE.R - Random walk Metropolis algorithm for epidemic in progress with daily aggregated count data. Allows for choice of likelihood function (exact GSE, birth-death, approximate birth-death) 

daily_convert.R - Program for converting aggregated daily count data to inter-removal data and to update the inter-removal times.

MCMC_Run_Abakiliki_V2.R - Implementation of MCMC algorithm for examples mentioned in the paper.

Abak_grid_posterior - Uses numerical integration and the grid of likelihood values to compute posterior quantities with T=90 and T=46.99. 
