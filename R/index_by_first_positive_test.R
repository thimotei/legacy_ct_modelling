index_by_first_positive <- function(dt) {
  pos_test <- dt[
    pcr_res == 1, .SD[t == min(t)], by = "id"][,
    .(id, t_first_pos = t)
  ]
  dt <- dt[pos_test, on = "id"]
  dt[, t := t - t_first_pos]
  return(dt[])
}
