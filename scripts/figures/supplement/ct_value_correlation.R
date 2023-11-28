library(data.table)
library(stringr)
library(purrr)
library(lubridate)
library(tidyr)
library(here)
library(forcats)

source("scripts/setup/main.R")

dt_raw <- fread("data_raw/raw_data.csv")

# changing structure of cleaned data to wide, to help with correlation plots
dt_clean_ct_corr <- process_data(dt_raw, struct_arg = "wide")

# adding time between last vaccine dose and first positive test calculation
# for each individual
dt_clean_ct_corr <- t_last_exposure(dt_clean_ct_corr, imm_delay = 14)

# choosing whether to pool all VOC subvariants together 
# (e.g. B.1.617.2-like and AY.4-like are both considered "Delta") or
# whether we consider them as separate groups. Main analysis considers
# them as grouped. Supplementary analysis considers them separately
dt_clean_ct_corr <- voc_status_attribution(dt_clean_ct_corr)

obs_ct_corr <- subset_data(
  dt_clean_ct_corr, no_pos_swabs = 2, first_pos_min = -15)

p_orf1ab_n_gene <- obs_ct_corr %>% 
  ggplot(aes(x = ct_value, y = ct_n_gene)) + 
  geom_point() + 
  geom_smooth() + 
  labs(x = "ORF1ab Ct value",
       y = "N gene Ct value",
       title = "ORF1ab vs N gene") + 
  theme_minimal()

p_orf1ab_s_gene <- obs_ct_corr %>% 
  ggplot(aes(x = ct_value, y = ct_s_gene)) + 
  geom_point() + 
  geom_smooth() +
  labs(x = "ORF1ab Ct value",
       y = "S gene Ct value",
       title = "ORF1ab vs S gene") + 
  theme_minimal()

p_s_n_gene <- obs_ct_corr %>% 
  ggplot(aes(x = ct_n_gene, y = ct_s_gene)) + 
  geom_point() + 
  geom_smooth() + 
  labs(x = "N gene Ct value", 
       y = "S gene Ct value",
       title = "N vs S gene") + 
  theme_minimal()

p_all_ct_corr <- p_orf1ab_n_gene + p_orf1ab_s_gene + p_s_n_gene

ggsave("outputs/figures/pdfs/supplement/ct_correlation.pdf",
       p_all_ct_corr,
       width = 14,
       height = 6)

ggsave("outputs/figures/pngs/supplement/ct_correlation.png",
       p_all_ct_corr,
       width = 14,
       height = 6)

