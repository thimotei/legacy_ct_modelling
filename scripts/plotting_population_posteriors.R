# Need to run inference script currently. Once we have finished with the model
# structure, we can save the fits and just read them in
source("scripts/inference.R")

draws <- rbind(as.data.table(fit_delta$draws())[, voc := "delta"],
               as.data.table(fit_omicron$draws())[, voc := "omicron"])

# draws[variable %like% "inc_mean"] %>% 
#   ggplot(

#--- gathering and plotting population-level posteriors by VOC
pop_params_rel <- c("c_0", "c_p_mean", "c_s_mean", 
                    "t_p_mean", "t_s_mean", "t_lod_mean")

pop_posteriors_wide <- transform_pop_posteriors(pop_params_rel, draws)

pop_params_abs <- c("c_0_abs", "c_p_mean_abs", "c_s_mean_abs", 
                    "t_p_mean_abs", "t_s_mean_abs", "t_lod_mean_abs")

pop_params_abs_labs <- c("Ct value at zero & LOD", "Ct value at peak", 
                         "Ct value at switch", "Time of peak", "Time of switch", 
                         "Time of limit of detection")

pop_params_abs_melt <- c("iteration", "voc", pop_params_abs)

pop_posteriors_long <- melt(pop_posteriors_wide[, ..pop_params_abs_melt], 
                            measure.vars = pop_params_abs)

p_pop_posteriors <- plot_pop_posteriors(
  pop_posteriors_long %>% 
    setnames(., "voc", "VOC") %>% 
    .[, variable := factor(variable, 
                           levels = pop_params_abs,
                           labels = pop_params_abs_labs)]
  )


#--- gathering and plotting population-level posterior predictive curves by VOC
pop_pp_samples_dt <- pop_pp_samples(pop_posteriors_wide, t_max = 40, t_step = 0.1)
pop_pp_summary_dt <- summarise_draws(pop_pp_samples_dt, by = c("time", "voc"))

p_pred_summary <- plot_pop_pp(pop_pp_summary_dt,
                              pop_pp_samples_dt,
                              no_samples = 100) 


#--- combining posteriors and posterior predictive plots
p_final <- p_pop_posteriors / p_pred_summary + plot_layout(guides = "collect")
  
ggsave("outputs/figures/posteriors_and_predictive.png",
       p_final,
       width = 8,
       height = 8,
       bg = "white")
