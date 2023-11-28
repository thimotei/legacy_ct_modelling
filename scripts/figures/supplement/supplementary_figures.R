# p2_fig2 <- ct_values_interactions(obs[,
#   VOC := fct_relevel(VOC,
#                      c("Delta", "Omicron (BA.1)", "Omicron (BA.2)"))],
#   vaccine_arg = TRUE)/
#   ct_values_interactions(obs[,
#     VOC := fct_relevel(VOC, c("Delta",
#                               "Omicron (BA.1)",
#                               "Omicron (BA.2)"))][,
#   age_group := fct_relevel(age_group, c("20-34"))
#   ],
#                          vaccine_arg = FALSE)
# 
# fig_2_top <- plot_grid(p2_fig2, p1_fig2, nrow = 1)
# 
# # p3_fig2 <- time_between_last_neg_first_pos(obs)
# p3_fig2 <- ct_values_no_interactions(obs)
# p4_fig2 <- plot_raw_dynamics(obs)
# 
# fig_2_bottom <- plot_grid(p3_fig2, p4_fig2, ncol = 2)
# 
# fig_2 <- plot_grid(fig_2_top, fig_2_bottom, nrow = 2, rel_heights = c(2, 1))
# # supplementary plots
# 
# # all raw data faceted by id
# tmp1 <- obs_all_cts[VOC == "Delta"] %>% 
#   ggplot() +
#   geom_point(aes(x = t, y = ct_value, colour = ct_type)) + 
#   # geom_smooth() + 
#   facet_rep_wrap(vars(data_id)) +
#   # facet_rep_grid(vars(VOC, data_id)) +
#   scale_y_reverse() +
#   custom_plot_theme()
# 
# tmp2 <- obs_all_cts[VOC == "Omicron"] %>% 
#   ggplot() +
#   geom_point(aes(x = t, y = ct_value, colour = ct_type)) + 
#   # geom_smooth() + 
#   facet_rep_wrap(vars(data_id)) +
#   # facet_rep_grid(vars(VOC, data_id)) +
#   scale_y_reverse() +
#   custom_plot_theme()
# 
# tmp3 <- obs_all_cts[VOC == "BA.2"] %>% 
#   ggplot() +
#   geom_point(aes(x = t, y = ct_value, colour = ct_type)) + 
#   # geom_smooth() + 
#   facet_rep_wrap(vars(data_id)) +
#   # facet_rep_grid(vars(VOC, data_id)) +
#   scale_y_reverse() +
#   custom_plot_theme()