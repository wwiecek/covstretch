fdf_palette <- c("grey20", "#E69F00", "#56B4E9")
# fdf_palette <- c("grey80", "#E69F00", "#56B4E9")
fdf_palette_text <- c("black", "white", "white")

# FDF model function (generating main metrics for all results) -----
model_fdf.all_k <- function(model, d1, e, policy, comp = c("cumI", "D")) {
  pars <- grab_2d_parms(model)
  
  dfdf <- as.numeric(fdf_deltas$d_fdf[fdf_deltas$d == d1])
  dh1 <- as.numeric(fdf_deltas$d_sse_h1[fdf_deltas$d == d1])
  dh2 <- as.numeric(fdf_deltas$d_sse_h2[fdf_deltas$d == d1])
  dh3 <- as.numeric(fdf_deltas$d_sse_h3[fdf_deltas$d == d1])
  dh4 <- as.numeric(fdf_deltas$d_sse_h4[fdf_deltas$d == d1])
  dh5 <- as.numeric(fdf_deltas$d_sse_h5[fdf_deltas$d == d1])
  dh6 <- as.numeric(fdf_deltas$d_sse_h6[fdf_deltas$d == d1])
  dh7 <- as.numeric(fdf_deltas$d_sse_h7[fdf_deltas$d == d1])
  dh8 <- as.numeric(fdf_deltas$d_sse_h8[fdf_deltas$d == d1])
  if(is.infinite(d1)){
    dfdf <- Inf; dh1 <- Inf; dh2 <- Inf; dh3 <- Inf; dh4 <- Inf; dh5 <- Inf; dh6 <- Inf; dh7 <- Inf; dh8 <- Inf
  }
  # if(d2 == 1) browser()
  if(policy == "default")
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    d1, delay_default))
  if(policy == "fdf")
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dfdf, delay_fdf))
  if(policy == "hybrid_1"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh1, delay_hybrid_k[,1]))
  }
  if(policy == "hybrid_2"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh2, delay_hybrid_k[,2]))
  }
  if(policy == "hybrid_3"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh3, delay_hybrid_k[,3]))
  }
  if(policy == "hybrid_4"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh4, delay_hybrid_k[,4]))
  }
  if(policy == "hybrid_5"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh5, delay_hybrid_k[,5]))
  }
  if(policy == "hybrid_6"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh6, delay_hybrid_k[,6]))
  }
  if(policy == "hybrid_7"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh7, delay_hybrid_k[,7]))
  }
  if(policy == "hybrid_8"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    dh8, delay_hybrid_k[,8]))
  }
  if(is.infinite(d1) || policy == "no_vaccination")
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

# KEEP AN EYE OUT FOR NUMERICAL ISSUES WITH THE ODEs HERE
df_fdf.all_k <- expand.grid(d1 = c(fdf_speeds,Inf),
                            model = scenario_par_nms_2v, 
                            e = c(.4, .5, .6, .7, .8, .9, .95), 
                            policy = c("default", "fdf", "hybrid_1", "hybrid_2", "hybrid_3", "hybrid_4",
                                       "hybrid_5", "hybrid_6", "hybrid_7", "hybrid_8")) %>%
  mutate(data = pmap(list(model, d1, e, policy), 
                     function(x,y,z,a) data.frame(value = model_fdf.all_k(x,y,z,a), 
                                                  var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) %>%
  group_by(d1, model, e, policy)


# optimal solution -----

fig2.all_k <- df_fdf.all_k %>% 
  filter(d1 > 40, d1 <= 1000) %>%
  filter(e >= .5) %>%
  select(d1, model, e, policy, d, harm, i) %>%
  gather(var, value, -d1, -model, -e, -policy) %>%
  group_by(d1, model, e, var) %>% 
  summarise(value_m = min(value[policy != "default"])/value[policy == "default"], 
            policy = policy[which.min(value)]
  ) %>%
  mutate(better = ifelse(value_m<=1&value_m>=0.95,1,NA)) %>% 
  mutate(better = factor(better, levels=c(1), labels = c("Less than 5% better"))) %>% 
  mutate(policy = factor(policy)) %>% #, levels = c("default", "fdf", "hybrid"), 
  #labels = c("Default (4 wks delay)", 
  #          "FDF (12 wks for all)", 
  #          "S-FDF (4 wks for 60+, 12 for rest)"))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  # mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 50%")) %>%
  mutate(delta1 = factor(1/d1,
                         levels = rev(1/fdf_speeds),
                         labels = as.percent(rev(1/fdf_speeds)))) %>%
  mutate(e = factor(e)) %>%
  mutate(value_m = round(value_m, 2)) %>%
  filter(var != "Economic harm") %>%
  ggplot(aes(x = delta1, y = e, fill = policy, alpha=better)) + 
  geom_tile() +
  # scale_fill_viridis_d(name="") +  
  # scale_fill_manual(values = c("grey20", "grey40", "grey60"),
  #                   name = "") +
  
  scale_fill_manual(values = fdf_palette, name = "") +
  scale_alpha_manual(values = c(0.7,1), name = "",labels=NULL,guide = 'none') +
  theme(legend.position = "bottom") +
  facet_grid(var~model) + 
  ylab("e1 (efficacy following 1st dose)") +
  theme(axis.text.x = element_text(angle = 45),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),text=element_text(size=8)) +
  xlab(paste0(def_labels, " (1st dose, default policy)"))

fig2s.all_k <- fig2.all_k + geom_text(aes(label = value_m), color = "white", alpha=1, size = 2)
# + scale_color_manual(values = fdf_palette_text, name = "") 
fig2s.all_k

fig_folder <- "figures"
ggsave(paste0(fig_folder, "/fdf_best_policy_allk.pdf"), fig2s.all_k, width = 6.5, height=0.6*6.5)
