model_i <- function(model, d1, e, rm = FALSE) {
  pars <- grab_2v_parms(model)
  pars <- apap_2v(pars, d1)
  # pars <- list_modify(pars, delta1 = rep(1/d1, Ngroups))
  y <- sr(list_modify(pars, e1 = e), "2v_v2")
  if(rm) return(y)
  main_metrics(y, pop)
}

# c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "40%" = 60, "49%" = 45)

df_efficacy_delta_raw <- expand_grid(d1 = default_speeds,
                                     e = seq(.5, .95, .05),
                                     model = scenario_par_nms_2v) %>%
  mutate(data = pmap(list(model, d1, e), function(x,y, z) data.frame(value = model_i(x,y, z), 
                                                                     var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e)  %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 


saveRDS(df_efficacy_delta_raw, file = "results/df_efficacy_delta_raw.rds")
