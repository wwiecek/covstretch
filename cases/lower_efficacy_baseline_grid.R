library(stats)
library(MASS)

select <- dplyr::select

nab_factor <- c(0.2,0.4,0.8,1)
speedups <- c(1,2/3,1/2,1/3,1/4)
baseline_effs <- c(0.7,0.95)

curve <- read.csv('data/curve.csv',header = FALSE)
colnames(curve) <- c("x","y")

##Interpolating curve----
model.eff <- approxfun(curve$x,curve$y)
model.neutr <- approxfun(curve$y,curve$x)

##Iterating over baseline efficacies----
list_plots <- vector('list', length(baseline_effs))
le.base_eff_grid <- data.frame()
for (i in 1:length(baseline_effs)){
  base_eff <- baseline_effs[i]
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
  
  df_efficacy_delta_raw.base_eff_grid <- expand_grid(d1 = c(speedups/default_delta_value),
                                       e = unique(effs),
                                       model = scenario_par_nms_2v) %>%
    mutate(data = pmap(list(model, d1, e), function(x,y, z) data.frame(value = model_i(x,y, z), 
                                                                       var = metric_nms))) %>%
    unnest(data) %>%
    spread(var, value) %>%
    mutate(nab_factor = factor(round(e,3),levels=round(effs,3),labels=nab_factor_mod)) %>% 
    group_by(model, e)  %>%
    mutate(d_rel=1-d/max(d)) %>%
    mutate(model = factor(model, levels = scenario_par_nms_2v,
                          labels = scenario_nms_2v)) 
  
  df_efficacy_delta.base_eff_grid <- df_efficacy_delta_raw.base_eff_grid %>%
    mutate(delta1 = 1/d1) %>% 
    select(delta1, e, model, i,d,d_rel,harm, nab_factor) %>%
    gather(var, value, -delta1, -e, -model, -nab_factor) %>%
    mutate(raw_value = value) %>%
    group_by(model,var) %>%
    mutate(e_ref=base_eff) %>% 
    mutate(ref = value[e == base_eff & delta1 == default_delta_value]) %>%
    ungroup() %>%
    mutate(r = (value/ref)) %>% 
    mutate(nab_factor_label=round(as.numeric(as.character(nab_factor)),1)) %>% 
    mutate(le_better = cut(r, c(-Inf, 1, Inf), 
                           labels = c("Fractional dose better",
                                      "Full dose better"
                           ))) %>%
    mutate(var = factor(var, levels = c("i", "d"),#,"d_rel"
                        labels = c("Infections", "Deaths"))) %>%#,"Relative Deaths"
    mutate(value = round(r, 2)) %>%
    mutate(speedup = factor(round(delta1/default_delta_value, 1))) %>%
    mutate(nab_e = factor(paste0(nab_factor_label,'(',round(100*e),'%)'))) %>%
    mutate(e = factor(round(e,2))) %>%
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
      ylab("NAb ratio (associated efficacy)") +
      xlab("Dose fraction") +
      geom_text(aes(label = format(value,3)), color = "white", size = 2.5)
  list_plots[[i]] <- local(print(df_efficacy_delta.base_eff_grid.plot))
}

#Matrix of results
le_baselines<-ggarrange(list_plots[[2]]+ggtitle("Baseline 95%"),list_plots[[1]]+ggtitle("Baseline 70%"), 
              common.legend = TRUE, ncol = 1, legend = "bottom")

#2d plots
le_2d_baselines <- le.base_eff_grid %>% filter((nab_factor_label==1&speedup=="1")|(nab_factor_label==0.8&speedup=="1/2")|
                              (nab_factor_label==0.4&speedup=="1/3")|(nab_factor_label==0.2&speedup=="1/10")) %>% 
  filter(var=='Relative Deaths') %>% 
  mutate(nab_factor_label=factor(nab_factor_label,levels=c(1,0.8,0.4,0.2),
                                 labels=c("Full dose (reference)","1/2 dose, 0.8 NAb ratio",
                                          "1/3 dose, 0.4 NAb ratio","1/10 dose, 0.2 NAb ratio"))) %>% 
  ggplot(aes(x = round(e_ref*100), y = 100*raw_value, color = nab_factor_label)) + 
  geom_line() + facet_grid(model~.) +
  xlab("Efficacy with full dose") +
  ylab("% burden averted") +
  lightness(scale_color_brewer(palette = "YlOrRd",direction = 1),scalefac(0.85))+
  theme(legend.position = "right")+labs(colour="Scenario")
le_2d_baselines

