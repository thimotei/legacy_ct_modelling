# library(cowplot)
# library(ggridges)
# library(lemon)
# library(ggsci)
# library(bayestestR)
# 
# # either run inference script or load in saved fit
# # fit <- readRDS("outputs/fits/fit_full.rds")
# 
# # run if hasn't been run yet (other figure plotting scripts require it too)
# draws <- extract_draws(fit)
# 
# # this takes a long time. this has already been run for figure 3, so
# # it is commented out by default. reduce the no_draws to 100 or so for
# # a rough but much faster version of the figures
# pop_ct_draws_trim <- extract_pop_ct_trajectories(fit,
#                                                  no_draws = 10000)
# 
# # list of factors in custom order
# factor_order <- c("Delta", "Omicron (BA.1)", "Omicron (BA.2)", 
#                   "Symptomatic", "Asymptomatic",
#                   "3 exposures", "4 exposures", "5+ exposures",
#                   "Age: 20-34", "Age: 35-49", "Age: 50+")
# 
# # extracting cumulative shedding draws
# shedding_draws <- extract_shedding(pop_ct_draws_trim,
#                                    vl_flag = TRUE,
#                                    trim_flag = TRUE,
#                                    pcr_pos_threshold = 40)
# 
# # summarising cumulative shedding trajectories with 95% credible interval 
# shedding_draws_plot <- summarise_shedding_trajectories(shedding_draws)
# 
# # extracting incubation period draws
# ip_draws <- extract_ip_draws(draws,
#                              factor_order = factor_order,
#                              keep_asymptomatic = FALSE,
#                              trim_flag = TRUE)
# 
# # calculating incubation period summary
# # ip_summary <- ip_draws[, .(ip_median = quantile(inc_mean_nat, 0.5)),
# #                        by = c("predictor", "regressor_category")]
# 
# shedding_plot_dt <- merge(shedding_draws_plot,
#                           ip_draws,
#                           by = c(".draw", "predictor", "regressor_category"), 
#                           all.x = TRUE)
# 
# dt_shedding_calc <- merge.data.table(shedding_draws,
#                                      ip_draws,
#                                      by = c(".draw",
#                                             "predictor",
#                                             "regressor_category"))
# 
# dt_shedding_tmp <- dt_shedding_calc[, .SD[which.min(abs(t - ip_draw))],
#                                     by = c(".draw", "predictor", "regressor_category")] 
# 
# dt_shedding_tmp %>% 
#   .[, .(me = quantile(auc_prop, 0.5), 
#         lo95 = quantile(auc_prop, 0.025),
#         hi95 = quantile(auc_prop, 0.975)), by = c("predictor", "regressor_category")] %>% 
#     ggplot() +
#   geom_pointrange(aes(x = predictor,
#                         y = me,
#                         ymin = lo95,
#                         ymax = hi95,
#                         colour = predictor)) +
#     scale_x_discrete(limits = rev) +
#     coord_flip() +
#     # scale_colour_brewer(palette = "Set1") +
#     theme_minimal() +
#     theme(axis.title.y = element_blank(),
#           legend.position = "none") +
#     scale_y_continuous(labels = scales::percent) +
#     labs(y = "Cumulative shedding (%)")
# 
# 
# # shedding_draws_sum <- summarise_shedding_trajectories(shedding_draws)
# 
# 
# # plotting the estimates in the same format as figure 3
# # p11_fig4 <- plot_cummulative_shedding_panel(shedding_draws_sum,
# #                                             ip_draws,
# #                                             "VOC",
# #                                             factor_arg = c("Delta",
# #                                                            "Omicron (BA.1)",
# #                                                            "Omicron (BA.2)"),
# #                                             title_arg = "VOC")
# # 
# # p21_fig4 <- shedding_effect_size_panel(shedding_effects,
# #                                        "VOC",
# #                                        factor_arg = c("Delta",
# #                                                       "Omicron (BA.1)",
# #                                                       "Omicron (BA.2)"),
# #                                        "Omicron (BA.1)")
# # 
# # p12_fig4 <- plot_cummulative_shedding_panel(shedding_draws_sum,
# #                                            ip_draws,
# #                                            "Symptom status",
# #                                            factor_arg = c("Symptomatic", 
# #                                                           "Asymptomatic"),
# #                                            title_arg = "Symptom Status")
# # 
# # p22_fig4 <- shedding_effect_size_panel(shedding_effects,
# #                                        "Symptom status",
# #                                        factor_arg = c("Symptomatic", 
# #                                                       "Asymptomatic"),
# #                                        "Symptomatic")
# # 
# # p13_fig4 <- plot_cummulative_shedding_panel(shedding_draws_sum,
# #                                             ip_draws,
# #                                             "Age",
# #                                             factor_arg = c("Age: 20-34",
# #                                                            "Age: 35-49",
# #                                                            "Age: 50-64",
# #                                                            "Age: 65+"),
# #                                             title_arg = "Age group")
# # 
# # p23_fig4 <- shedding_effect_size_panel(shedding_effects,
# #                                        "Age",
# #                                        factor_arg = c("Age: 20-34",
# #                                                       "Age: 35-49",
# #                                                       "Age: 50-64",
# #                                                       "Age: 65+"),
# #                                        "Age: 35-49")
# # 
# # p14_fig4 <- plot_cummulative_shedding_panel(shedding_draws_sum,
# #                                             ip_draws,
# #                                             "Number of exposures",
# #                                             factor_arg = c("3 exposures", 
# #                                                            "4 exposures",
# #                                                            "5 exposures",
# #                                                            "6 exposures"),
# #                                             title_arg = "Number of exposures")
# # 
# #                            
# # p24_fig4 <- shedding_effect_size_panel(shedding_effects,
# #                                        "Number of exposures",
# #                                        factor_arg = c("3 exposures", 
# #                                                       "4 exposures",
# #                                                       "5 exposures",
# #                                                       "6 exposures"),
# #                                        "4 exposures")
# 
# # p4_1 <- figure_4_subpanel(p11_fig4, p21_fig4, rel_widths_arg = c(2, 1))
# # p4_2 <- figure_4_subpanel(p12_fig4, p22_fig4, rel_widths_arg = c(2, 1))
# # p4_3 <- figure_4_subpanel(p13_fig4, p23_fig4, rel_widths_arg = c(2, 1))
# # p4_4 <- figure_4_subpanel(p14_fig4, p24_fig4, rel_widths_arg = c(2, 1))
# # 
# # p4 <- plot_grid(p4_1, p4_2, p4_3, p4_4, nrow = 2)
# # 
# # ggsave("outputs/figures/figure_4_tmp.png",
# #        p4,
# #        width = 15,
# #        height = 8,
# #        bg = "white")
# # 
# # ggsave("outputs/figures/figure_4_tmp.pdf",
# #        p4,
# #        width = 15,
# #        height = 8,
# #        bg = "white")
