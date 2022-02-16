# function to condition full dataset on VOC type and number of positive swabs
# by individual. Also reducing dataset to just the data needed for model
# inference

subset_data <- function(dt_clean_in, voc, no_pos_swabs) {
  dt_out <- dt_clean[VOC %in% voc] %>%
    .[, t_first_test := as.numeric(swab_date - min(swab_date), units = "days"),
      by = c("id", "infection_id")] %>%
    .[no_pos_results > no_pos_swabs] %>%
    .[, id := .GRP, by = id] %>%
    .[, c("id", "infection_id", "swab_date", "t_first_test", "t", "ct_value",
          "onset_time", "result", "pcr_res")]

  return(dt_out)
}
