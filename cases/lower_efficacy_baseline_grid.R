library(stats)
library(MASS)

select <- dplyr::select

nab_factor <- c(0.2,0.4,0.8,1)
speedups <- c(1,2/3,1/2,1/3,1/4)

curve <- read.csv('data/curve.csv',header = FALSE)
colnames(curve) <- c("x","y")

#Interpolating curve
model.eff <- approxfun(curve$x,curve$y)
model.neutr <- approxfun(curve$y,curve$x)

#Iterating over baseline efficacies
le.base_eff_grid <- data.frame()
for (base_eff in c(0.7,0.95)){
  base_neutr <- model.neutr(base_eff*100)
  effs <- c()
  nab_factor_mod <- nab_factor
  for (factor in nab_factor){
    if (is.na(model.eff(factor*base_neutr))){
      effs <- c(effs, min(curve$y))
      nab_factor_mod[nab_factor_mod==factor] <- round(curve[curve$y==min(curve$y),]$x/base_neutr,2)
    } else {
      effs <- c(effs, model.eff(factor*base_neutr))
    }
  }
  effs <- effs/100
  
  df_efficacy_delta_raw.base_eff_grid <- expand_grid(d1 = speedups/default_delta_value,
                                       e = effs,
                                       model = scenario_par_nms_2v) %>%
    mutate(data = pmap(list(model, d1, e), function(x,y, z) data.frame(value = model_i(x,y, z), 
                                                                       var = metric_nms))) %>%
    unnest(data) %>%
    spread(var, value) %>%
    mutate(nab_factor = factor(round(e,3),levels=round(effs,3),labels=nab_factor_mod)) %>% 
    group_by(model, e)  %>%
    mutate(model = factor(model, levels = scenario_par_nms_2v,
                          labels = scenario_nms_2v)) 
  
  df_efficacy_delta.base_eff_grid <- df_efficacy_delta_raw.base_eff_grid %>%
    mutate(delta1 = 1/d1) %>% 
    select(delta1, e, model, i,d,harm, nab_factor) %>%
    gather(var, value, -delta1, -e, -model, -nab_factor) %>%
    group_by(model,var) %>%
    mutate(e_ref=base_eff) %>% 
    mutate(ref = value[e == base_eff & delta1 == default_delta_value]) %>%
    ungroup() %>%
    mutate(r = (value/ref)) %>% 
    mutate(le_better = cut(r, c(-Inf, 1, Inf), 
                           labels = c("Fractional dose better",
                                      "Full dose better"
                           ))) %>%
    mutate(var = factor(var, levels = c("i", "d"),
                        labels = c("Infections", "Deaths"))) %>%
    mutate(value = round(r, 2)) %>%
    mutate(speedup = factor(round(delta1/default_delta_value, 1))) %>%
    mutate(e = factor(round(e,2))) %>%
    mutate(nab_e = factor(paste0(nab_factor,'-',e))) %>%
    mutate(speedup = factor(round(delta1/default_delta_value, 4),
                            levels = round(1/speedups,4),
                            labels = as.character(fractions(speedups)))) %>% 
    filter(var != "Economic harm")
  
  le.base_eff_grid <- rbind(le.base_eff_grid,df_efficacy_delta.base_eff_grid)
  
  df_efficacy_delta.base_eff_grid.plot <- df_efficacy_delta.base_eff_grid %>%
    ggplot(aes(x = speedup, y = nab_e, fill = le_better)) + geom_tile() +
      scale_fill_manual(values = c("grey60", "grey40", "grey20"),
                        name = "") +
      theme(legend.position = "bottom", axis.text.x = element_text(hjust = 1,angle = 45)) +
      facet_grid(var~model) +
      ylab("NAb variation of fractional dose") +
      xlab("Dose fraction") +
      geom_text(aes(label = format(value,3)), color = "white", size = 2.5)
  ggsave(paste0(fig_folder, paste0("/le_optimal_baseline_",str_replace(as.character(base_eff), "0.", ""),".pdf")),
        df_efficacy_delta.base_eff_grid.plot, width = width, height=3.4/5.55*width)
}
