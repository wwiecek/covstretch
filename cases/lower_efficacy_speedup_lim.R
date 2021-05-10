default_speeds <- round(100/c(2, rev(seq(.05, 1, .05)), .025, .01, 0), 5) # % per day
for (c in c(seq(1,3,0.05),seq(3.5,32,0.5))){
  default_speeds <- unique(c(default_speeds,
                             round(1/(c(0.001,0.0025,0.005,0.0075,0.01,0.02)*c),5)))
  
}

source("cases/prep-results.R")

le_ref <- data.frame()
for (ref_delta in c(0.001,0.0025,0.005,0.0075,0.01,0.02)){
  le_ref <- rbind(le_ref, df_efficacy_delta_raw %>%
                    filter(round(d1,5) %in% round(as.numeric(1/(ref_delta*c(seq(1,3,0.05),seq(3.5,32,0.5)))),5)) %>%
                    mutate(delta1 = round(1/d1,5)) %>%
                    select(delta1, e, model, i,d,harm) %>%
                    gather(var, value, -delta1, -e, -model) %>%
                    group_by(model,var) %>%
                    mutate(ref = value[e == .95 & delta1 == ref_delta]) %>%
                    mutate(ref_delta = ref_delta) %>%
                    filter(e %in% c(seq(.5, .9, .1), .95)) %>%
                    ungroup() %>%
                    mutate(r = (value/ref)) %>% 
                    mutate(var = factor(var, levels = c("i", "d", "harm"),
                                        labels = c("Infections", "Deaths", "Economic harm"))) %>%
                    mutate(speedup = round(delta1/ref_delta, 2)) %>%
                    # mutate(speedup = factor(as.percent(delta1, 2))) %>%
                    mutate(e = factor(e)) %>%
                    filter(var != "Economic harm"))# & speedup %in% c(1, 1.2, 1.4, 1.6, 2, 4, 8, 16, 32)))
}

#Grouped epidemic scenarios (max)
le_ref.fig <- le_ref %>% 
  select(var,e,model,ref_delta,speedup,r) %>% 
  group_by(var,e,ref_delta,speedup) %>% 
  mutate(max_r = max(r)) %>%
  ungroup() %>% 
  filter(max_r<=1) %>% 
  select(-r,-max_r) %>% 
  group_by(var,e,model,ref_delta) %>% 
  mutate(min_speedup = round(min(speedup),2)) %>%
  group_by(var,e,ref_delta) %>% 
  mutate(min_speedup = max(min_speedup)) %>%
  mutate(ref_delta = factor(ref_delta)) %>%
  select(-model,-speedup) %>%
  unique() %>%
  ggplot(aes(x = ref_delta, y = e,fill=min_speedup)) + geom_tile() +
  scale_fill_viridis_c(begin = 0.5,end = 0) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(var~.) + 
  labs(fill = "Minimum speedup for benefits") +
  ylab("Efficacy for the alternative dose") + 
  xlab("Reference vaccination speed") + 
  scale_x_discrete(breaks = c(0.001,0.0025,0.005,0.0075,0.01), labels = as.percent(c(0.001,0.0025,0.005,0.0075,0.01))) +
  geom_text(aes(label = min_speedup), color="white", size = 3)

le_ref.fig

ggsave(paste0(fig_folder, "/speedup_lim_grouped.pdf"),le_ref.fig + theme(text = element_text(size=9)), 
       width = width, height=8/5.55*width)

#Separate epidemic scenarios
le_ref.fig <- le_ref %>% 
  dplyr::select(var,e,model,ref_delta,speedup,r) %>% 
  group_by(var,e,model,ref_delta,speedup) %>% 
  mutate(max_r = max(r)) %>%
  ungroup() %>% 
  filter(max_r<=1) %>% 
  dplyr::select(-r,-max_r) %>% 
  group_by(var,e,model,ref_delta) %>% 
  mutate(min_speedup = round(min(speedup),2)) %>%
  mutate(ref_delta = factor(ref_delta)) %>%
  dplyr::select(-model,-speedup) %>%
  unique() %>%
  ggplot(aes(x = ref_delta, y = e,fill=min_speedup)) + geom_tile() +
  lightness(scale_fill_distiller(palette = "Greys",direction = 1, name = "Minimum speedup for benefits"),scalefac(0.6))+
  # scale_fill_viridis_c(begin = 0.5,end = 0) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(var~model) + 
  ylab("Efficacy of fractional dose") + 
  xlab("Percentage of pop. vaccinated daily (baseline)") + 
  scale_x_discrete(breaks = c(0.001,0.0025,0.005,0.0075,0.01,0.02), labels = as.percent(c(0.001,0.0025,0.005,0.0075,0.01,0.02))) +
  geom_text(aes(label = min_speedup), color="white", size = 2)

le_ref.fig

ggsave(paste0(fig_folder, "/speedup_lim_separate_speedup.pdf"),le_ref.fig + theme(text = element_text(size=10)), 
       width = width, height=3.4/5.55*width)
