# Relative benefit of 
model_i <- function(model, d1, e, vrf = 1, rm = FALSE) {
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
                                     vrf = c(0,.05,1),
                                     model = scenario_par_nms_2v) %>%
  mutate(data = pmap(list(model, d1, e, vrf), function(x,y,z,v) data.frame(value = model_i(x,y, z,v), 
                                                                     var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e, vrf)  %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 

df2 <- 
  df %>%
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
  group_by(d1, model, e)

df2 %>% 
  filter(e == .95) %>%
  filter(d1 %in% c(d1_general, Inf)) %>% 
  mutate(d1 = as.percent(1/d1)) %>%
  select(d1, model, i, ri, d, rd, vrf) %>%
  mutate(d = round(d*1e05)) %>%
  gather(variable, value, -d1, -model, -e, -vrf) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "harm", "ri", "rd", "re", "diffd", "diffi"),
                           labels = c("Infections", "Deaths per 100,000", "Economic harm",
                                      "RI", "RD", "RE", "Difference in deaths", "Difference in infections"))) %>%
  mutate(value = round(value, 3)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e) %>%
  View
