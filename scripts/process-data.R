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
dt_proc <- process_data(dt_raw, struct_arg = "long")

# adding time between last vaccine dose and first positive test calculation
# for each individual
dt_proc <- t_last_exposure(dt_proc, imm_delay = 14)

# choosing whether to pool all VOC subvariants together 
# (e.g. B.1.617.2-like and AY.4-like are both considered "Delta") or
# whether we consider them as separate groups. Main analysis considers
# them as grouped. Supplementary analysis considers them separately
dt_proc <- voc_status_attribution(dt_proc)

# Do additional processing to filter for the desired number of swabs
# per positive episode
obs <- subset_data(dt_proc,
                   no_pos_swabs = 2,
                   first_pos_min = -15)

cols_to_remove <- c("centre", "sex", "age", "barcode", "Sequenced",
                    "dose_1", "dose_2", "dose_3", "dose_4")

cols_to_keep <- setdiff(colnames(obs), cols_to_remove)

dt_clean <- obs[, ..cols_to_keep] 

# saving processed data as R object file
saveRDS(dt_clean, here("data", "processed-data.rds"))
