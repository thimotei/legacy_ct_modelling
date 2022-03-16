#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(stringr)
library(purrr)
library(lubridate)
library(tidyr)
library(here)

# loading all functions in package directory
devtools::load_all()

dt_raw <- fread("data/raw-data.csv")

# processing data, adding machine readable dates, moving dates
# back for certain swabs, based on data collection advice and adjusting
# all wet swabs upwards according to the fitted dry vs wet adjustment
# factor
dt_clean <- process_data(dt_raw)

summary(dt_clean)

saveRDS(dt_clean, here("data", "processed-data.rds"))