#--- simulating Ct trajectories, to test parameter recovery
library(data.table)
library(stringr)
library(purrr)
library(lubridate)
library(tidyr)
library(here)

# loading all functions in package directory
devtools::load_all()

dt_raw <- fread("data/ct_values_clean.csv")

# processing data, adding machine readable dates, moving dates
# back for certain swabs, based on data collection advice and adjusting
# all wet swabs upwards according to the fitted dry vs wet adjustment
# factor
dt_clean <- process_data(dt_raw)

#--- TEMPORARY, for now, as "invalid" and "inconclusive" both have many
#--- observations consistent with positive Ct values, keeping them for now
#--- will discuss with Crick partners about this
dt_clean[result == "Invalid", result := "Positive"]
dt_clean[result == "Inconclusive", result := "Positive"]
dt_clean[result == "Negative", pcr_res := 0]
dt_clean[result == "Positive", pcr_res := 1]

#--- Drop patients with gaps between positive tests of more than 60 days
ids_spurious_gaps <- dt_clean[,
  .(spurious = any(abs(time_since_first_pos) > 60) || any(abs(onset_time) > 60)
   ), by = "id"
]
ids_spurious_gaps[spurious == TRUE][]
dt_clean <- dt_clean[ids_spurious_gaps, on = "id"]
dt_clean <- dt_clean[spurious == FALSE | is.na(spurious)]

#--- Drop subject with ID 978 due to large mismatch between onset and positive
#--  tests
dt_clean <- dt_clean[id != 978]

fwrite(dt_clean, here("data", "processed-data.csv"))