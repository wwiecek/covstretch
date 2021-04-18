# Relative benefit of not vaccinating recovered population 

model_vr <- function(model, d1, e, vrf = 1, rm = FALSE) {
  pars <- grab_2v_parms(model) %>% 
    apap_2v(d1)
  if(vrf == 1)
    y <- sr(list_modify(pars, e1 = e, vrf = 1), "2v_v2")
  else if(vrf == 0)
    y <- sr(list_modify(pars, e1 = e, vrf = 0), "2v_v2")
  else{
    pars <- apap_2v(pars, d1*(1-vrf))
    y <- sr(list_modify(pars, e1 = e, vrf = 1), "2v_v2")
  }
  if(rm) return(y)
  main_metrics(y, pop)
}


df <- expand_grid(d1 = default_speeds,
                                     e = .95,
                                     vrf = c(0,1),
                                     model = scenario_par_nms_2v) %>%
  mutate(data = pmap(list(model, d1, e, vrf), function(x,y,z,v) data.frame(value = model_vr(x,y, z,v), 
                                                                     var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e, vrf)  %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 

# df2 <- 
#   df %>%
#   filter(d1 %in% c(d1_general, Inf)) %>%
#   select(-harm, -tt50, -v1) %>%
#   pivot_wider(names_from = c("vrf"), id_cols = c("d1", "e", "model"), values_from = c("d", "i")) %>%
#   mutate(bd = d_0/d_1, bi = i_0/i_1) %>%
#   group_by(d1, model, e)
# 
# df2 %>% 
#   filter(d1 %in% c(d1_general, Inf)) %>% 
#   mutate(d1 = as.percent(1/d1)) %>%
#   select(-d_0, -d_1, -i_0, -i_1) %>%
#   gather(key, value, -d1, -e, -model) %>%
#   spread(d1, value)
# 
# df2 %>% 
#   filter(e == .95) %>%
#   filter(d1 %in% c(d1_general, Inf)) %>% 
#   mutate(d1 = as.percent(1/d1)) %>%
#   select(d1, model, i, ri, d, rd, vrf) %>%
#   mutate(d = round(d*1e05)) %>%
#   gather(variable, value, -d1, -model, -e, -vrf) %>%
#   mutate(variable = factor(variable, 
#                            levels = c("i", "d", "harm", "ri", "rd", "re", "diffd", "diffi"),
#                            labels = c("Infections", "Deaths per 100,000", "Economic harm",
#                                       "RI", "RD", "RE", "Difference in deaths", "Difference in infections"))) %>%
#   mutate(value = round(value, 3)) %>%
#   spread(d1, value) %>% 
#   arrange(variable, model) %>% ungroup() %>%
#   select(-e)


# Prep figure

tab_screening <- df %>% 
  filter(d1 %in% c(d1_general, Inf)) %>% 
  # mutate(d1 = as.percent(1/d1)) %>%
  mutate(d1 = (1/d1)) %>%
  group_by(model) %>% mutate(d = 1 - i/max(i)) %>%
  select(d1, model, d, vrf) %>%
  spread(d1, d)
 


fig_screening <- df %>% 
  filter(d1 %in% c(d1_general, Inf)) %>% 
  # mutate(d1 = as.percent(1/d1)) %>%
  mutate(d1 = (1/d1)) %>%
  group_by(model) %>% mutate(d = 1 - i/max(i)) %>%
  mutate(vrf = factor(vrf, labels = c("yes", "no"), levels = c(0,1))) %>%
  ggplot(aes(x = d1, y = d, linetype = factor(vrf), color = model)) + geom_line() +
  xlab("Vaccinated per day (%)") + ylab("Fraction of infections averted") + 
  scale_linetype_discrete(name = "screening") +
  scale_color_discrete(name = "epidemic scenario") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general))
