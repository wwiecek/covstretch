df1 <- lapply(scenario_list_2v, function(pars) {
  ll <- list(
    "Vaccinate 0.3% per day" = apap_2v(pars, 100/.3),
    "Vaccinate 0.5% per day" = apap_2v(pars, 100/.5),
    "Vaccinate 1% per day"   = apap_2v(pars, 100/1),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(list_modify, kappa1 = rep(0, Ngroups)) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind()
  as.data.frame(ll[,"I",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario")

df2 <- lapply(scenario_list_2v, function(pars) {
  ll <- list(
    "Vaccinate 0.3% per day" = apap_2v(pars, 100/.3),
    "Vaccinate 0.5% per day" = apap_2v(pars, 100/.5),
    "Vaccinate 1% per day"   = apap_2v(pars, 100/1),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(list_modify, kappa1 = rep(1/360, Ngroups)) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind()
  as.data.frame(ll[,"I",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario")

fig_kappa <- rbind(mutate(df1, kappa = "No immunity loss"),
      mutate(df2, kappa = "1 year immunity")) %>%
  mutate(scenario = factor(scenario, 
                           levels = c("Slow decrease", "Slow growth", "Fast growth"))) %>% 
  gather(var, value, -time, -scenario, -kappa) %>%
  # filter(scenario != "Fast growth") %>%
  ggplot(aes(x = time, y = value*100, color = var, linetype = kappa)) + geom_line() + facet_wrap(.~scenario, scales = "free") +
  xlab("Time (days)") + scale_x_continuous(breaks = seq(0, 360, 120)) + ylab("% currently infected") +
  lightness(scale_color_brewer(palette = "YlOrRd",direction = 1,name = ""),scalefac(0.9)) +
  scale_linetype_discrete(name = "") +
  theme(legend.position = "top", 
        legend.direction = "vertical")

# Compare two kappa levels ------

model_kappa <- function(model, d1, e, kappa, rm = FALSE) {
  pars <- grab_2v_parms(model)
  pars <- apap_2v(pars, d1)
  # pars <- list_modify(pars, delta1 = rep(1/d1, Ngroups))
  y <- sr(list_modify(pars, e1 = e, kappa1 = rep(kappa, Ngroups)), "2v_v2")
  if(rm) return(y)
  main_metrics(y, pop)
}
  
df_kappa <- expand_grid(d1 = c(d1_general, Inf),
                  e = .95,
                  kappa = c(0, 1/360),
                  model = scenario_par_nms_2v) %>%
  mutate(data = pmap(list(model, d1, e, kappa), function(x,y,z,k) data.frame(value = model_kappa(x,y,z,k), 
                                                                             var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e, kappa) 

# This needs a table for appendix?
df_kappa %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) %>%
  filter(d1 %in% c(d1_general, Inf)) %>%
  mutate(ref_e = harm[d1 > 1460]) %>%
  mutate(ref_i = i[d1 > 1460]) %>%
  mutate(ref_d = d[d1 > 1460]) %>%
  ungroup() %>%
  mutate(re = 1-(harm/ref_e)) %>%
  mutate(ri = 1-(i/ref_i)) %>%
  mutate(diffi = (i-ref_i)) %>%
  mutate(rd = 1-(d/ref_d)) %>%
  mutate(diffd = (d-ref_d)) %>%
  group_by(d1, model, e, kappa) %>%
  filter(e == .95) %>%
  mutate(d1 = as.percent(1/d1)) %>%
  select(d1, model, i, ri, d, rd, diffi, diffd) %>%
  mutate(d = round(d*1e05)) %>%
  mutate(diffd = round(diffd*1e05)) %>%
  gather(variable, value, -d1, -model, -e, -kappa) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "harm", "ri", "rd", "re", "diffd", "diffi"),
                           labels = c("Infections", "Deaths per 100,000", "Economic harm",
                                      "RI", "RD", "RE", "Difference in deaths", "Difference in infections"))) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e) 
