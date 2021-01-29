# Going from 100 million doses to 200 million in 3 months


us_fdf_slow <- lst(
  Nc = 13, 
  Ngroups, 
  Ndays = 365,
  y0 = y0_gen(13, 9, pre_immunity+.05, 1e-02),
  q = rep(1.5/(5*ev), Ngroups),
  contacts = default_cm,
  gamma1 = rep(.2, Ngroups),
  gamma2 = rep(.2, Ngroups), #duration of infectious period
  delta1 = rep(1/300, Ngroups),
  delta2 = rep(0, Ngroups),
  kappa1 = rep(0, Ngroups),
  kappa2 = rep(0, Ngroups),
  phi = rep(0, Ngroups), 
  ta = rep(0, Ngroups),
  e1 = .8,
  e2 = .95,
  pdeath = c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100,
  # vrf = 1,
  constantrisk = 0
)

us_fdf_fast <- list_modify(us_fdf_slow,
                             y0 = y0_gen(13, 9, pre_immunity+.05, 1e-03),
                             q = rep(2.5/(5*ev), Ngroups))
us_fdf_cr <- list_modify(us_fdf_slow,
                           y0 = y0_gen(13, 9, pre_immunity+.05, 0),
                           q = rep(0, Ngroups),
                           constantrisk = .01/30.5)

# S1: Base case (Figure 1) ------
models <- list(
  "no vaccination"      = sr(apap(us_fdf_cr, Inf, delay_default), f = "2d_v2"),
  "current speed"       = sr(apap(us_fdf_cr, 540, delay_default), f = "2d_v2"),
  "2x faster"           = sr(apap(us_fdf_cr, 270, delay_default), f = "2d_v2"),
  "2x faster, with FDF" = sr(apap(us_fdf_cr, 213, delay_fdf), f = "2d_v2")
)
lapply(models, main_metrics, pop, 100) %>% bind_rows(.id = "scenario") %>% 
  select(scenario, i, d, harm_vr, v1) %>%
  gather(key, value, -scenario) %>% 
  spread(scenario, value) %>%
  mutate(key = c("deaths per 100", "economic harm", "infections per 100", "doses used")) %>%
  gather(var, value, -key) %>% spread(key, value)


lapply(models, main_metrics, pop, 100) %>% bind_rows(.id = "scenario") %>% 
  select(scenario, i, d, harm_vr, v1) %>%
  gather(key, value, -scenario) %>% 
  spread(scenario, value) %>%
  mutate(key = c("deaths per 100", "economic harm", "infections per 100", "doses used")) %>%
  mutate(`2x faster` = `2x faster`/`no vaccination`) %>%
  mutate(`2x faster, with FDF` = `2x faster, with FDF`/`no vaccination`) %>%
  mutate(`current speed` = `current speed`/`no vaccination`) %>%
  select(-`no vaccination`) %>% filter(c(T,T,T,F))

rab <- rescale_and_bind(models, pop) 

rab[100, ,]

rab %>% 
  plot_rcs(comp_to_display, long_names = ln, ncol = 3) + 
  ylab("Proportion of population") + theme(legend.position = "top") + 
  scale_color_discrete(name = "scenario")




