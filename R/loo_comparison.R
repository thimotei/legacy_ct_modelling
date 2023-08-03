# performing a LOO-CV using the PSIS approximation implemented in the loo
# package, designed to be used with Stan

# For the LOO-CV to work, the likelihood must be evaluated for each 
# observation. We do this in the generated quantities block of the Stan model.
# Once this is implemented and the log likelihood is returned in the fit object,
# the $loo() function which is included in the cmdstanr package can be directly
# run, returning the elpd and loo-ic values for model comparison

# we load and run a loo comparison for the 3 models of interest:

# semi-informative priors, with VOC as a covariate (main model)
# semi-informative priors, without VOC as a covariate
# semi-informative priors, wider individual-level variation
# uninformative priors

library(loo)

fit_main <- readRDS("outputs/fits/fit_main.rds")
fit_no_voc <- readRDS("outputs/fits/fit_no_voc.rds")

loo_main <- fit_main$loo()
loo_no_voc <- fit_no_voc$loo()

# printing the psis-loo results to add them to the supplementary table
loo_main
loo_no_voc

loo_compare(loo_main, loo_no_voc)
