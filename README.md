# DON-LCL
 DON LCL Population analysis

The whole scripts are performed under rstan package.

The population concentration-response raw data of each cell line is first input in the stan formula and fitted with the stan model. The Rhat distribution is used to check the convergence. The fitting curves are then visualized in plots to examine the curve trend with 95% confidence interval. A 10% decreased change in concentration-response curve is derived to represent the EC10 distribution for each cell line. The geometric standard deviation (GSD) distribution of EC10 is also used to derive toxicodynamic variability factor (TDVF).

## stan fitting

input file:
- DON-LCLdata.csv
- DONstan_dat.R

script:
  DON-LCL-fitting.R
  
output: 
- DON-LCLdata_test.csv
- DONstan_dat.R
- DON_stanfit.RData
- DON-Rhat.pdf

## Concentration-response fitting

input file:
- DON-LCLdata_test.csv
- DON_stanfit.RData

original source from script:
- DON-LCL-fitting.R

script:
  DON-LCL-plotfits.R

output: 
- Figure-DON_fit.pdf (Supplementary Figure S6)
  
## Effective concentration (EC10) derivation

input file:
- DON-LCLdata_test.csv
- DON_stanfit.RData

original source from script:
- DON-LCL-fitting.R

script: 
  DON-LCL-plot_n_and_EC10.R

output:
- Table-DON_n_CellLine_summary.csv
- Figure-DON_n.pdf
- DON_EC10_cellline_samples.csv (Manuscript  Figure 5)
- DON-EC10_pop_samples.csv
- Table-EC10.pop_posteriors.csv
- Figure-DON-EC10_pop.GM.pdf
- Figure-DON-EC10_pop.GSD.pdf
- Figure-DON-EC10_pop.GSD_compare_prior.pdf
- Figure-DON-EC10_pop.TDVF05_compare_prior.pdf
- Figure-DON-EC10_pop.TDVF01_compare_prior.pdf
- TDVF_figure.pdf (Supplementary Figure S7)
- Table-EC10.pop_posteriors.csv
  