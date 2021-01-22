# Base case settings -----
library(tidyverse)
library(ggpubr)
theme_set(theme_minimal(base_size = 18))
source("R/ode_2doses.R")
source("R/ode_2vaccines.R")
source("R/helpers.R")
source("R/output-helpers.R")
source("R/config.R")
source("R/config-pars.R")


# Table with BI and BD for various efficacy and delta1 values -----
# Reference columns for RI and RD are severity of epidemic without vaccination

model_i <- function(model, d1, e, comp = c("cumI", "D")) {
  if(model == "pars_le_cr")   pars <- pars_le_cr
  if(model == "pars_le_slow")   pars <- pars_le_slow
  if(model == "pars_le_fast") pars <- pars_le_fast
  sr(list_modify(pars, 
                 e1 = e,
                 delta1 = rep(1/d1, Ngroups))) %>% b_any(pop, comp)
}

# c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "40%" = 60, "49%" = 45)

df_efficacy_delta <- expand_grid(d1 = c(seq(10, 360, 10), 730, 1460, Inf),
                                 e = seq(.5, .95, .05),
                                 model = c("pars_le_cr", "pars_le_slow", "pars_le_fast")) %>%
  mutate(data = pmap(list(model, d1, e), function(x,y, z) data.frame(value = model_i(x,y, z), 
                                                                     var = c("i", "d")))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e) %>%
  mutate(ref_i = i[d1 > 1460]) %>%
  mutate(ref_d = d[d1 > 1460]) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                        labels = c("Constant risk", "Slow growth", "Fast growth"))) %>%
  mutate(ri = 1 - (i/ref_i)) %>%
  mutate(rd = 1 - (d/ref_d)) %>%
  group_by(d1, model, e)
