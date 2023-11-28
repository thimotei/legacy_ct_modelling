# performing a LOO-CV using the PSIS approximation implemented in the loo
# package, designed to be used with Stan

# For the LOO-CV to work, the likelihood must be evaluated for each 
# observation. We do this in the generated quantities block of the Stan model.
# Once this is implemented and the log likelihood is returned in the fit object,
# the $loo() function which is included in the cmdstanr package can be directly
# run, returning the elpd and loo-ic values for model comparison

library(loo)

fit_main <- readRDS("outputs/fits/main.rds")
fit_voc <- readRDS("outputs/fits/voc.rds")
fit_symptoms <- readRDS("outputs/fits/symptoms.rds")
fit_age <- readRDS("outputs/fits/age.rds")
fit_exposures <- readRDS("outputs/fits/exposures.rds")
fit_no_voc <- readRDS("outputs/fits/no_voc.rds")
fit_no_covariates <- readRDS("outputs/fits/no_covariates.rds")
fit_uninformative <- readRDS("outputs/fits/uninformative.rds")

loo_main <- fit_main$loo(cores = 8)
loo_voc <- fit_voc$loo(cores = 8)
loo_symptoms <- fit_symptoms$loo(cores = 8)
loo_exposures <- fit_exposures$loo(cores = 8)
loo_age <- fit_age$loo(cores = 8)
loo_no_voc <- fit_no_voc$loo(cores = 8)
loo_no_covariates <- fit_no_covariates$loo(cores = 8)
loo_uninformative <- fit_uninformative$loo(cores = 8)

# printing the psis-loo results to add them to the supplementary table
loo_compare(loo_main, 
            loo_voc, 
            loo_symptoms, 
            loo_exposures,
            loo_age,
            loo_no_voc, 
            loo_no_covariates,
            loo_uninformative)

# printing all results to write up the table reporting the values
loo_main
loo_voc
loo_symptoms
loo_exposures
loo_age
loo_no_voc
loo_no_covariates
loo_uninformative
