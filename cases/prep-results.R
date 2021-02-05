source("setup.R")

model_i <- function(model, d1, e, rm = FALSE) {
  if(model == "pars_le_cr")   pars <- pars_le_cr
  if(model == "pars_le_slow") pars <- pars_le_slow
  if(model == "pars_le_fast") pars <- pars_le_fast
  
  pars <- apap_2v(pars, d1)
  # pars <- list_modify(pars, delta1 = rep(1/d1, Ngroups))
  y <- sr(list_modify(pars, e1 = e), "2v_v2")
  if(rm) return(y)
  main_metrics(y, pop)
}

# c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "40%" = 60, "49%" = 45)

df_efficacy_delta_raw <- expand_grid(d1 = c(seq(60, 360, 10), 450, 540, 630, 730, 1460, Inf),
                                     e = seq(.5, .95, .05),
                                     model = c("pars_le_cr", "pars_le_slow", "pars_le_fast")) %>%
  mutate(data = pmap(list(model, d1, e), function(x,y, z) data.frame(value = model_i(x,y, z), 
                                                                     var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e)  %>%
  mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                        labels = c("Constant risk", "Slow growth", "Fast growth"))) 


saveRDS(df_efficacy_delta_raw, file = "results/df_efficacy_delta_raw.rds")
