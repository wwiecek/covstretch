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
model_i <- function(model, d1, e) {
  if(model == "pars_le_cr")   pars <- pars_le_cr
  if(model == "pars_le_slow") pars <- pars_le_slow
  if(model == "pars_le_fast") pars <- pars_le_fast
  
  # pars <- apap_2v(pars, d1)
  pars <- list_modify(pars, delta1 = rep(1/d1, Ngroups))
  y <- sr(list_modify(pars, e1 = e))
  main_metrics(y, pop)
}

# c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "40%" = 60, "49%" = 45)

df_efficacy_delta_raw <- expand_grid(d1 = c(seq(10, 360, 10), 730, 1460, Inf),
                                     e = seq(.5, .95, .05),
                                     model = c("pars_le_cr", "pars_le_slow", "pars_le_fast")) %>%
  mutate(data = pmap(list(model, d1, e), function(x,y, z) data.frame(value = model_i(x,y, z), 
                                                                     var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e) 

df_efficacy_delta <- 
  df_efficacy_delta_raw %>%
  # mutate(benefit_vr = 1 - harm_vr) %>%
  mutate(ref_e = harm[d1 > 1460]) %>%
  mutate(ref_i = i[d1 > 1460]) %>%
  mutate(ref_d = d[d1 > 1460]) %>%
  ungroup() %>%
  mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                        labels = c("Constant risk", "Slow growth", "Fast growth"))) %>%
  mutate(re = 1 - (harm/ref_e)) %>%
  mutate(ri = 1 - (i/ref_i)) %>%
  mutate(rd = 1 - (d/ref_d)) %>%
  group_by(d1, model, e)




# Figure G1: general ilustration of the model and benefits of vaccination -----
mlist <- list(
  "Constant risk of infection" = pars_le_cr,
  "Base case (R0 = 1.5)" = pars_le_slow,
  "Fast growth (R0 = 3)" = pars_le_fast) %>%
  setNames(scenario_names)

ln2 <- ln
ln2[["cumV1"]] <- "Courses of vaccine used"
gglist <- lapply(as.list(1:3), function(i) {
  pars <- mlist[[i]]
  list(
    "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
    # "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
    "Mass vaccination (no prioritsation), delta = 1/60" = list_modify(pars, delta1 = rep(1/60, Ngroups)),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind() %>% 
    # plot_rcs(c("S", "I", "cumV1", "D"), ncol = 5, long_names = ln2,
    plot_rcs(c("S", "I", "P1", "P2", "D"), ncol = 5, long_names = ln2,
             end_date = as.Date("01-01-2021", format="%d-%m-%Y") + 300) + ylab("") +
    ggtitle(names(mlist)[i])
})

ggarrange(plotlist=gglist, common.legend = TRUE, ncol = 1, legend = "top")



# Fig G2: How many infections and deaths averted with 95% efficacious vaccine -----
d1_general <- c(60, 90, 180, 360, 730, 1460)
df_efficacy_delta  %>%
  filter(e == .95) %>%
  # filter(d1 < 340) %>%
  filter(d1 %in% d1_general) %>%
  select(ri, rd, re) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("RI (infections)", "RD (deaths)", "RE (economic harm)"))) %>%
  ggplot(aes(x = d1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = d1_general[-1]) +
  scale_color_discrete(name = "scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 12), legend.position = "top") +
  xlab("average time to vaccination, 1/ delta1 [days]") + ylab("proportion of harm averted")



# Table G3 -----
df_efficacy_delta %>% 
  filter(e == .95) %>%
  filter(d1 %in% c(d1_general, Inf)) %>% 
  # select(d1, model, i, ri, benefit_r, v1, tthi) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)


