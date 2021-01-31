# Base case settings -----
library(tidyverse)
library(ggpubr)
theme_set(theme_minimal(base_size = 18))
source("R/ode_2doses.R")
source("R/ode_2doses_v2.R")
source("R/ode_2vaccines.R")
source("R/ode_2vaccines_v2.R")
source("R/helpers.R")
source("R/output-helpers.R")
source("R/config.R")
source("R/config-pars.R")
source("R/harm_function.R")
source("R/prioritisation.R")


# Table with BI and BD for various efficacy and delta1 values -----
# Reference columns for RI and RD are severity of epidemic without vaccination
model_i <- function(model, d1, e) {
  if(model == "pars_le_cr")   pars <- pars_le_cr
  if(model == "pars_le_slow") pars <- pars_le_slow
  if(model == "pars_le_fast") pars <- pars_le_fast
  
  pars <- apap_2v(pars, d1)
  # pars <- list_modify(pars, delta1 = rep(1/d1, Ngroups))
  y <- sr(list_modify(pars, e1 = e), "2v_v2")
  main_metrics(y, pop)
}

# c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "40%" = 60, "49%" = 45)

# df_efficacy_delta_raw <- expand_grid(d1 = c(1, seq(90, 360, 10), 730, 1460, Inf),
df_efficacy_delta_raw <- expand_grid(d1 = c(90, 180, 360, 730, 1460, Inf),
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
    "Mass vaccination: 360 days" = apap_2v(pars, 360),
    # "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
    "Mass vaccination: 180 days" = apap_2v(pars, 180),
    "Mass vaccination: 90 days" = apap_2v(pars, 90),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind() %>% 
    # plot_rcs(c("S", "I", "cumV1", "D"), ncol = 5, long_names = ln2,
    plot_rcs(c("S", "I", "cumV1", "D"), ncol = 5, long_names = ln2,
             end_date = as.Date("01-01-2021", format="%d-%m-%Y") + 300) + ylab("") +
    ggtitle(names(mlist)[i])
})

ggarrange(plotlist=gglist, common.legend = TRUE, ncol = 1, legend = "top")



# Fig G2B: How many infections and deaths averted with 95% efficacious vaccine -----
d1_general <- c(60, 90, 180, 360, 730, 1460)
g2b <- df_efficacy_delta  %>%
  filter(e == .95) %>%
  # filter(d1 < 340) %>%
  filter(d1 %in% d1_general) %>%
  select(ri, rd, re) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("infections", "deaths", "economic harm"))) %>%
  ggplot(aes(x = d1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = d1_general[-1]) +
  scale_color_discrete(name = "scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 12), legend.position = "top") +
  xlab("length of vaccination campaign, 1/delta") + 
  ylab("fraction of harm averted")

# Fig G2A: absolute harm
g2a <- df_efficacy_delta  %>%
  filter(e == .95) %>%
  # filter(d1 < 340) %>%
  filter(d1 %in% d1_general) %>%
  select(i, d, harm) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("i", "d", "harm"), 
                      labels = c("infections", "deaths", "economic harm"))) %>%
  ggplot(aes(x = d1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = d1_general[-1]) +
  scale_color_discrete(name = "scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 12), legend.position = "top") +
  xlab("length of vaccination campaign, 1/delta") + 
  ylab("total harm")

ggarrange(g2a+ggtitle("A"), 
          g2b+ggtitle("B"), 
          common.legend = TRUE, ncol = 1, legend = "top")

# Table G3 -----
df_efficacy_delta %>% 
  filter(e == .95) %>%
  filter(d1 %in% c(d1_general, Inf)) %>% 
  select(d1, model, i, ri, d, rd, harm, re) %>%
  mutate(d = round(d*1e05)) %>%
  # mutate(i = i*1e05) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "harm", "ri", "rd", "re"),
                           labels = c("Infections", "Deaths per 100,000", "Economic harm",
                                      "RI", "RD", "RE"))) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)


