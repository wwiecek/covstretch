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


comp_to_display <- c("I", "D", "cumV", "cumI", "P1", "P2")
delay_default <- 18
delay_fdf <- 74
delay_hybrid <- c(rep(delay_fdf, 6), rep(delay_default, 3))



# FINDING VALUES FOR DELTA2 that fit DELTA1 -----
# d1_default <- c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "49%" = 60, "100%" = 1)
d1_default <- c(1, 60, 90, 120, 150, 180, 270, 360, 540, 730, 1460)

# d1_default <- c("4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "49%" = 45, "100%" = 1)
d_default <- sapply(d1_default, function(d){
  print(d)
  bsl <- sr(apap_2d(pars_fdf_slow, d, delay_default), f = "2d_v2") %>% 
    rescale_rcs(pop, merge = T)
  # bsl_vaccinated <- bsl[c(30, 60, 90, 180, 360),"cumV",1]
  
  # d_start <- floor(d*.25)
  d_start <- 1
  bsl_vaccinated <- bsl[d_start:min(d, 360),"cumV",1]
  
  sse <- sapply(ceiling(d/2):d, function(d2){
    fdf <- sr(apap_2d(pars_fdf_slow, d2, delay_fdf), f = "2d_v2") %>% 
      rescale_rcs(pop, merge = T)
    fdf_v <- fdf[d_start:min(d, 360),"cumV",1]
    # fdf_v <- fdf[c(30, 60, 90, 180, 360),"cumV",1]
    sum((fdf_v - bsl_vaccinated)^2)
  })
  sse_h <- sapply(ceiling(d/2):d, function(d2){
    fdf <- sr(apap_2d(pars_fdf_slow, c(rep(d2, 6), rep(d,3)), delay_hybrid), f = "2d_v2") %>% 
      # fdf <- sr(apap_2d(pars_fdf_slow, d2, delay_hybrid), f = "2d_v2") %>% 
      rescale_rcs(pop, merge = T)
    fdf_v <- fdf[d_start:min(d, 360),"cumV",1]
    sum((fdf_v - bsl_vaccinated)^2)
  })
  c(d,
    which.min(sse) + ceiling(d/2) - 1,
    which.min(sse_h) + ceiling(d/2) - 1)
})

# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# [1,]    1   60   90  120  150  180  270  360  540   730  1460
# [2,]    1   46   69   93  118  145  230  318  481   650  1300
# [3,]    1   52   77  102  128  155  239  325  496   680  1460
# > 


# FDF model function (generating main metrics for all results) -----

model_fdf <- function(model, d1, e, policy, comp = c("cumI", "D")) {
  if(model == "pars_le_cr")   pars <- pars_fdf_cr
  if(model == "pars_le_slow")   pars <- pars_fdf_slow
  if(model == "pars_le_fast") pars <- pars_fdf_fast
  d2 <- as.numeric(d_default[2, d1_default == d1])
  d3 <- as.numeric(d_default[3, d1_default == d1])
  if(is.infinite(d1)){
    d2 <- Inf; d3 <- Inf
  }
  # if(d2 == 1) browser()
  if(policy == "default")
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    d1, delay_default))
  if(policy == "fdf")
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    d2, delay_fdf))
  if(policy == "hybrid"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    d3, delay_hybrid))
  }
  if(policy == "no_vaccination")
    res <- sr(f = "2d_v2", 
              list_modify(pars, 
                          e1 = 0,
                          e2 = 0,
                          ta = rep(0, Ngroups),
                          delta1 = rep(0, Ngroups),
                          delta2 = rep(0, Ngroups)))
  main_metrics(res, pop, vat = 91)
}



# Generate a data frame with all values -----

df_fdf <- expand.grid(d1 = c(d1_default, Inf),
                      model = c("pars_le_cr", "pars_le_slow", "pars_le_fast"), 
                      e = c(.4, .5, .6, .7, .8, .9, .95), 
                      policy = c("default", "fdf", "hybrid")) %>%
  mutate(data = pmap(list(model, d1, e, policy), 
                     function(x,y,z,a) data.frame(value = model_fdf(x,y,z,a), 
                                                  var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                        labels = c("Constant risk", "Slow growth", "Fast growth"))) %>%
  group_by(d1, model, e, policy)




# Total doses used: illustration -----
lapply(as.list(1:length(d1_default)), function(i) {
  x <- d_default[,i]
  w <- rescale_and_bind(list(
    "Default (4 weeks)"      = sr(apap_2d(pars_fdf_cr, x[1], delay_default) %>% list_modify(e1 = .8), f = "2d_v2"),
    "FDF (12 weeks)"         = sr(apap_2d(pars_fdf_cr, x[2], delay_fdf) %>% list_modify(e1 = .8), f = "2d_v2"),
    "FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_cr, c(rep(x[3], 6), rep(x[1],3)), delay_hybrid) %>% list_modify(e1 = .8), f = "2d_v2")
    # "FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_cr, 156, delay_hybrid) %>% list_modify(e1 = .8), f = "2d_v2")
  ), pop)
  as.data.frame(w[,"cumV",]) %>% rownames_to_column("time")
}) %>% 
  setNames(d_default[1,]) %>%
  bind_rows(.id = "d") %>%
  mutate(time = as.numeric(time)) %>%
  mutate(d = as.numeric(d)) %>%
  gather(policy, value, -d, -time) %>%
  filter(d > 1) %>%
  ggplot(aes(x=time, y=value, color=policy)) + geom_line() + facet_wrap(~d) + ylim(0, .75)



# Fig 1: Base case example ------
rescale_and_bind(list(
  "Default (4 weeks)"        = sr(apap_2d(pars_fdf_slow, d_default[1,6], delay_default) %>% 
                                    list_modify(e1 = .8), f = "2d_v2"),
  "FDF (12 weeks)"           = sr(apap_2d(pars_fdf_slow, d_default[2,6], delay_fdf) %>% 
                                    list_modify(e1 = .8), f = "2d_v2"),
  "S-FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_slow, c(rep(d_default[3,6], 6), rep(d_default[1,6],3)), delay_hybrid) %>% 
                                    list_modify(e1 = .8), f = "2d_v2")
  # "FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_cr, 156, delay_hybrid) %>% list_modify(e1 = .8), f = "2d_v2")
), pop) %>% 
  plot_rcs(c("P1", "P2", "cumV"), long_names = ln, ncol = 3) + 
  ylab("Proportion of population") + theme(legend.position = "top") + 
  scale_color_discrete(name = "2nd dose delay policy:")

# sr(apap_2d(pars_fdf_fast, 180, delay_default) %>% list_modify(e1 = .8), f = "2d_v2") %>% main_metrics(pop)
# sr(apap_2d(pars_fdf_fast, 145, delay_fdf) %>% list_modify(e1 = .8), f = "2d_v2") %>% main_metrics(pop)



df_gg <- df_fdf %>% 
  filter(e %in% c(.5, .8)) %>% 
  filter(d1 > 60, d1 < 540) %>% 
  select(d1, model, e, policy, d, harm, i) %>%
  gather(var, value, -d1, -model, -e, -policy) %>%
  # mutate(d1 = factor(d1)) %>%
  mutate(policy = factor(policy, levels = c("default", "fdf", "hybrid"), 
                         labels = c("Default (4 wks delay)", 
                                    "FDF (12 wks for all)", 
                                    "S-FDF (4 wks for 60+, 12 for rest)"))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(e = factor(e))

df_gg %>%
  ggplot(aes(x = d1, y = value, lty = e, color = policy)) + 
  geom_line(size=1.1) +
  # geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(var~model, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45)) +
  scale_x_continuous(breaks = d1_default[c(3,6:10)]) +
  scale_linetype_manual(name = "efficacy after 1st dose", values = c("solid", "dotted")) +
  xlab("targeted vaccination speed (1/delta1 under default policy)") + ylab("")

df_gg %>%
  filter(model == "Slow growth", var %in% c("Infections", "Deaths")) %>%
  ggplot(aes(x = d1, y = value, lty = e, color = policy)) + 
  geom_line(size=1.1) +
  # geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(var~model, scales = "free") +
  theme(legend.position = "top", legend.direction = "vertical", axis.text.x = element_text(angle = 45)) +
  scale_x_continuous(breaks = d1_default[c(3,6:10)]) +
  scale_linetype_manual(name = "efficacy after 1st dose", values = c("solid", "dotted")) +
  xlab("targeted vaccination speed (1/delta1 under default policy)") + ylab("")



# Fig 4: optimal solution -----

fig4 <- df_fdf %>% 
  filter(d1 > 60, d1 < 1000) %>%
  select(d1, model, e, policy, d, harm, i) %>%
  gather(var, value, -d1, -model, -e, -policy) %>%
  group_by(d1, model, e, var) %>% 
  summarise(value_m = min(value[policy != "default"])/value[policy == "default"], 
            policy = policy[which.min(value)]
  ) %>%
  mutate(policy = factor(policy, levels = c("default", "fdf", "hybrid"), 
                         labels = c("Default (4 wks delay)", 
                                    "FDF (12 wks for all)", 
                                    "S-FDF (4 wks for 60+, 12 for rest)"))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  # mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 50%")) %>%
  mutate(speedup = factor(d1)) %>%
  mutate(e = factor(e)) %>%
  mutate(value_m = round(value_m, 2)) %>%
  ggplot(aes(x = speedup, y = e, fill = policy)) + 
  geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "") +
  theme(legend.position = "bottom") +
  facet_grid(var~model) + 
  ylab("e1 (efficacy following 1st dose)") +
  theme(axis.text.x = element_text(angle = 45)) +
  xlab("vaccination speed 1/delta1 [days]")

fig4
fig4 + geom_text(aes(label = value_m), color = "white")



# Table with all results -----
df_fdf_ratio <- df_fdf %>% select(d, i, harm) %>%
  pivot_wider(values_from = c("d", "i", "harm"), names_from = "policy") %>%
  # mutate(re = 1 - (harm_vr_fdf/harm_vr_default)) %>%
  # mutate(re = harm_vr_fdf - harm_vr_default) %>%
  # mutate(rd = d_fdf - d_default) %>%
  # mutate(ri = i_fdf - i_default) %>%
  mutate(re = (harm_fdf / harm_default)) %>%
  mutate(rd = (   d_fdf / d_default)) %>%
  mutate(ri = (   i_fdf / i_default)) %>%
  # mutate(ri_fn = 1 - (i_fdf/i_no_vaccination)) %>%
  # mutate(ri_dn = 1 - (i_default/i_no_vaccination)) %>%
  group_by(d1, model, e)


filter(df_fdf_ratio, e %in% c(.6, .8)) %>%
  filter(d1 %in% c(d1_default, Inf)) %>% 
  # select(d1, model, ri_fd, ri_fn, ri_dn, i_fdf, i_default) %>%
  # select(d1, model, i_fdf, i_default) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(e, variable, model) %>% 
  ungroup() 
# select(-e, -variable)



# Dependence between speed of vaccination and magniutude of benefits -----
df_fdf_ratio  %>%
  filter(e %in% c(.6, .8)) %>%
  filter(d1 > 60, d1 < 1000) %>%
  select(ri, rd, re) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("RRI (infections)", "RRD (deaths)", "RRE (economic harm)"))) %>%
  mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 60%")) %>%
  ggplot(aes(x = d1, y = value, group = interaction(model, efficacy), lty = efficacy, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap( ~ key, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = c(90, 180, 360, 730, 1460)) +
  scale_color_discrete(name = "scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 14), legend.position = "top") +
  xlab("average time to first dose, 1/ delta1 [days]") + 
  ylab("relative risk FDF vs default policy")




# Figure FDF4: optimal choice FDF vs default depending on speed and efficacy ------
df_fdf_ratio %>%
  filter(d1 > 1, d1 < Inf) %>%
  # mutate(rd = re) %>%
  select(d1, e, model, ri, re, rd) %>%
  gather(key, value, -d1, -e, -model) %>%
  mutate(rd = value) %>%
  mutate(fdf_better = cut(rd, c(-Inf, .95, 1.05, Inf), 
                          labels = c("FDF better by 5% or more", 
                                     "Comparable (+-5%)", 
                                     "FDF worse by at least 5%"))) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"),
                      labels = c("RRI (infections)", "RRD (deaths)", "RRE (economic harm)"))) %>%
  mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 50%")) %>%
  mutate(value = round(value, 2)) %>%
  mutate(speedup = factor(d1)) %>%
  mutate(e = factor(e)) %>%
  ggplot(aes(x = speedup, y = e, fill = fdf_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "Reduction in infections RI:") +
  theme(legend.position = "bottom") +
  facet_grid(key~model) + ylab("e1 (efficacy following 1st dose)") + 
  xlab("average wait until the first dose under default policy, 1/delta1 [days]") +
  geom_text(aes(label = value), color = "white")  





df_fdf %>% filter(model == "Constant risk", e == .95) %>% select(d1, model, policy, v1) %>% spread(policy, v1)
