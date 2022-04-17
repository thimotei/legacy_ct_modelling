#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(stringr)
library(purrr)
library(lubridate)
library(tidyr)
library(here)
library(forcats)

# loading all functions in package directory
devtools::load_all()

dt_raw <- fread("data/raw_data.csv")

# processing data, adding machine readable dates, moving dates
# back for certain swabs, based on data collection advice and adjusting
# all wet swabs upwards according to the fitted dry vs wet adjustment
# factor
dt_clean <- process_data(dt_raw)

# adding time between last vaccine dose and first positive test calculation
# for each individual
dt_clean <- t_last_dose(dt_clean, imm_delay = 14)

# choosing whether to pool all VOC subvariants together 
# (e.g. B.1.617.2-like and AY.4-like are both considered "Delta") or
# whether we consider them as separate groups. Main analysis considers
# them as grouped. Supplementary analysis considers them separately
dt_clean <- voc_status_attribution(dt_clean,
                                   group_all = TRUE)

summary(dt_clean)

saveRDS(dt_clean, here("data", "processed-data.rds"))
