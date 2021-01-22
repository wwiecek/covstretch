source("cases/prep-results.R")
comp_to_display <- c("I", "D", "cumV", "cumI", "P1", "P2")

# S1: Base case (Figure 1) ------

rescale_and_bind(list(
  "Default (4 weeks)" = sr(list_modify(pars_fdf_fast, 
                                       ta = rep(10, Ngroups),
                                       delta1 = rep(1/180, Ngroups),
                                       delta2 = rep(1/18, Ngroups)), f = "2d"),
  "FDF (12 weeks)" = sr(list_modify(pars_fdf_fast, 
                                    ta = rep(10, Ngroups),
                                    delta1 = rep(1/142, Ngroups),
                                    delta2 = rep(1/74, Ngroups)), f = "2d")
), pop) %>% plot_rcs(comp_to_display, long_names = ln, ncol = 3) + 
  ylab("Proportion of population") + theme(legend.position = "top") + 
  scale_color_discrete(name = "2nd dose delay policy:")



# FINDING VALUES FOR DELTA2 that fit DELTA1 -----
d1_default <- c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "49%" = 45, "100%" = 1)
d2_default <- sapply(d1_default, function(d){
  bsl <- sr(list_modify(pars_fdf_slow, 
                        ta = rep(10, Ngroups),
                        delta1 = rep(1/d, Ngroups),
                        delta2 = rep(1/18, Ngroups)), f = "2d") %>% 
    rescale_rcs(pop, merge = T)
  bsl_vaccinated <- bsl[c(30, 60, 90, 180, 360),"cumV",1]
  sse <- sapply(1:d, function(d2){
    fdf <- sr(list_modify(pars_fdf_slow, 
                          ta = rep(10, Ngroups),
                          delta1 = rep(1/d2, Ngroups),
                          delta2 = rep(1/74, Ngroups)), f = "2d") %>% 
      rescale_rcs(pop, merge = T)
    fdf_v <- fdf[c(30, 60, 90, 180, 360),"cumV",1]
    sum((fdf_v - bsl_vaccinated)^2)
  })
  which.min(sse)  
})
# > cbind(d1, d2)
# d1   d2
# 2%  1460 1302
# 4%   730  643
# 8%   360  308
# 15%  180  142
# 28%   90   61
# 49%   45   26


# c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "40%" = 60, "49%" = 45)
model_fdf <- function(model, d1, e, policy, comp = c("cumI", "D")) {
  if(model == "pars_le_cr")   pars <- pars_fdf_cr
  if(model == "pars_le_slow")   pars <- pars_fdf_slow
  if(model == "pars_le_fast") pars <- pars_fdf_fast
  d2 <- as.numeric(d2_default[d1_default == d1])
  # if(d2 == 1) browser()
  if(policy == "default")
    res <- sr(f = "2d", list_modify(pars, 
                   e1 = e,
                   e2 = .95, 
                   ta = rep(10, Ngroups),
                   delta1 = rep(1/d1, Ngroups),
                   delta2 = rep(1/18, Ngroups))) %>% b_any(pop, comp)
  if(policy == "fdf")
    res <- sr(f = "2d", list_modify(pars, 
                   e1 = e,
                   e2 = .95,
                   ta = rep(10, Ngroups),
                   delta1 = rep(1/d2, Ngroups),
                   delta2 = rep(1/74, Ngroups))) %>% b_any(pop, comp)
  res
}

# model_fdf("pars_le_slow", 1, .8, "fdf")
# model_fdf("pars_le_slow", 1, .8, "default")

df_fdf <- expand.grid(d1 = d1_default, 
            model = c("pars_le_cr", "pars_le_slow", "pars_le_fast"), 
            e = c(.4, .5, .6, .7, .8), 
            policy = c("default", "fdf")) %>%
  mutate(data = pmap(list(model, d1, e, policy), 
                     function(x,y,z,a) data.frame(value = model_fdf(x,y,z,a), 
                                                  var = c("i", "d")))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                        labels = c("Constant risk", "Slow growth", "Fast growth")))
  
df_fdf_ratio <- df_fdf %>%
  pivot_wider(values_from = c("d", "i"), names_from = "policy") %>%
  mutate(rd = 1 - (d_fdf/d_default)) %>%
  mutate(ri = 1 - (i_fdf/i_default)) %>%
  group_by(d1, model, e)


filter(df_fdf_ratio, e %in% c(.5, .8)) %>%
  filter(d1 %in% c(1, 45, 90, 180, 360, 730, 1460, Inf)) %>% 
  select(d1, model, ri, i_fdf, i_default) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(e, variable, model) %>% 
  ungroup() %>% 
  select(-e, -variable)

df_fdf_ratio  %>%
  filter(e %in% c(.5, .8)) %>%
  filter(d1 > 45) %>%
  select(ri, rd) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd"), 
                      labels = c("RI (infections)", "RD (deaths)"))) %>%
  mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 50%")) %>%
  ggplot(aes(x = d1, y = value, group = interaction(model, efficacy), color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(efficacy~key, scales = "free") +
  scale_x_continuous(breaks = c(90, 180, 360, 730, 1460)) +
  scale_color_discrete(name = "scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 14), legend.position = "top") +
  xlab("average time to first dose, 1/ delta1 [days]") + 
  ylab("relative benefit R")




# Figure FDF4: optimal choice FDF vs default depending on speed and efficacy ------
df_fdf_ratio %>%
  filter(d1 > 1) %>%
  select(d1, e, model, rd) %>%
  mutate(fdf_better = cut(rd, c(-Inf, -.05, .05, Inf), 
                         labels = c("FDF worse by 5% or more", 
                                    "Comparable (+-5%)", 
                                    "FDF better by at least 5%"))) %>%
  mutate(value = round(rd, 2)) %>%
  mutate(speedup = factor(d1)) %>%
  mutate(e = factor(e)) %>%
  ggplot(aes(x = speedup, y = e, fill = fdf_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "Reduction in infections RI:") +
  theme(legend.position = "bottom") +
  facet_grid(~model) + ylab("e1 (efficacy following 1st dose)") + 
  xlab("average wait until the first dose under default policy, 1/delta1 [days]") +
  geom_text(aes(label = value), color = "white")  



# Loss of immunity ------
