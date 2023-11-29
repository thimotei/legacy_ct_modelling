# Sourcing setup script to load in data
source("scripts/setup/main.R")

# Munging observations from long to wide format to help with plotting
dt_ct_values_correlation <- obs[, .(id, swab_date, ct_value, ct_type)][
  , sequence := seq_len(.N), by = .(id, swab_date, ct_type)] |> 
  dcast(id + swab_date + sequence~ ct_type, value.var = "ct_value")

# Plotting ORF1AB against N-gene values
p_orf1ab_n_gene <- dt_ct_values_correlation %>% 
  ggplot(aes(x = ct_value, y = ct_n_gene)) + 
  geom_point() + 
  geom_smooth() + 
  labs(x = "ORF1ab Ct value",
       y = "N gene Ct value",
       title = "ORF1ab vs N gene") + 
  theme_minimal()

# Plotting ORF1AB against S-gene values
p_orf1ab_s_gene <- dt_ct_values_correlation %>% 
  ggplot(aes(x = ct_value, y = ct_s_gene)) + 
  geom_point() + 
  geom_smooth() +
  labs(x = "ORF1ab Ct value",
       y = "S gene Ct value",
       title = "ORF1ab vs S gene") + 
  theme_minimal()

# Plotting N-gene against S-gene values
p_s_n_gene <- dt_ct_values_correlation %>% 
  ggplot(aes(x = ct_n_gene, y = ct_s_gene)) + 
  geom_point() + 
  geom_smooth() + 
  labs(x = "N gene Ct value", 
       y = "S gene Ct value",
       title = "N vs S gene") + 
  theme_minimal()

# Putting all three plots together
p_all_ct_corr <- p_orf1ab_n_gene + p_orf1ab_s_gene + p_s_n_gene

# Saving the data for the figure
saveRDS(dt_ct_values_correlation, "outputs/plot_data/supplement/figure_S3.rds")

# Saving the final figure
ggsave("outputs/figures/pdfs/supplement/figure_S3.pdf",
       p_all_ct_corr,
       width = 14,
       height = 6)

ggsave("outputs/figures/pngs/supplement/figure_S3.png",
       p_all_ct_corr,
       width = 14,
       height = 6)

