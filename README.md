# Individual-level modelling of cycle threshold (Ct) dynamics
 
**Summary:** 
Understanding how different factors, such as novel variants of concern (VOCs) or prior exposures, influence SARS-CoV-2 viral shedding kinetics is crucial for interpreting past epidemic trends and developing plans to mitigate future transmission risk. However, as population immunity has accumulated via infection and vaccination, alongside emergence of multiple VOCs, it has become harder to generalise insights from one population to another. 

This repository outlines a Bayesian modelling pipeline to reconstruct unobserved PCR cycle threshold (Ct) trajectories for each individual, using detailed data from a UK prospective cohort undergoing weekly occupational health PCR screening for SARS-CoV-2. Full details of the model and analysis are given in the below manuscript.

*Citation*

Russell TW, Townsley H, Abbott S, Hellewell J, Carr EJ et al. Within-host viral kinetics of SARS-CoV-2 in a cohort with complex life course exposures reveals different intrinsic properties of Omicron sub-variants BA.1 and BA.2 compared to Delta.

## Repository organisation

### Folder structure

This repository follows the structure of an R package, with functions defined in `R` folder, study data in `data`, and stan model code defined in `stan`. Data loading and main analysis scripts are in `scripts`.

### Key script files

To reproduce the main analysis in the accompanying paper, use the following script files:

> `scripts/setup.R` - Load and visualise data, and define design matrix and adjustment model for swab type and the gene target

> `scripts/figure_1.R` - Generate schematic of model structure

> `scripts/figure_2.R` - Plot summary of data

> `scripts/figure_3.R` - Fit model and plot estimates of Ct dynamics and effect sizes for different covariates. Supplementary models variants are defined in `scripts/figure_3_[VARIANT].R`.

> `scripts/figure_4.R` - Fit model and plot estimates of incubation period against Ct growth dynamics for different covariates. Supplementary model variants are defined in `scripts/figure_4_[VARIANT].R`.

### Dependencies

Required libraries are listed in the script files. In particular, the inference framework requires [cmdstanr](https://mc-stan.org/cmdstanr/) to be installed, which in turn requires a suitable C++ toolchain. 

