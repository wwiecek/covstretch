source("project-setup.R")

theme_set(theme_minimal(base_size = 18))
pars_fdf_slow %>%
  apap_2d(400, 18) %>%
  sr(f = "2d_v2") %>% 
  plot_rcs(c("R", "S", "cumI", "cumV", "D"))

list_modify(pars_fdf_slow,
            Nc = 18,
            y0 = y0_gen(18, 9, pre_immunity, 5e-03, E = 2, I = 5, R = 8),
            pdeath = NULL,
            pd0 = ifr_hic,
            e1 = 0, e2 = 0,
            pd1 = (1-.95)*ifr_hic,
            pd2 = (1-.95)*ifr_hic) %>%
  apap_2d(400, 18) %>%
  sr(f = "2d_v3") %>% 
  plot_rcs(c("R", "RV", "cumI", "S", "cumV", "D"))

expand_2d_v3 <- function(pars, x, e1, e2) {
  pd <- pars$pdeath
  list_modify(pars,
              Nc = 18,
              y0 = y0_gen(18, 9, pre_immunity, x, E = 2, I = 5, R = 8),
              pdeath = NULL,
              pd0 = pd,
              pd1 = (1-e1)*pd,
              pd2 = (1-e2)*pd)
}

expand_2d_v3(pars_fdf_fast, infected0[["fast"]], 
             e1=.95,e2=.95) %>%
  apap_2d(400, group_seq = F) %>%
  list_modify(e1 = 0, e2 = 0) %>%
  sr("2d_v3") %>%
  # check0sums(14)
  # plot_rcs(c("R", "RV", "cumI", "S", "cumV", "D"))
  # plot_rcs(c("R", "RV", "I0", "I1", "I2", "S", "cumV", "D"))
  plot_rcs(c("R", "RV", "I0", "I1", "I2", "P1", "P2", "N1", "N2"))
  # b_any(pop, "cumI")





model_i <- function(model, d1, e_d, e_t, rm = FALSE) {
  if(model == "pars_linear")    
    pars <- expand_2d_v3(pars_fdf_linear, infected0[["decreasing"]], 
                         e1=e_d,e2=e_d)
  if(model == "pars_le_slow")   
    pars <- expand_2d_v3(pars_fdf_slow, infected0[["slow"]], 
                         e1=e_d,e2=e_d)
  if(model == "pars_le_fast")   
    pars <- expand_2d_v3(pars_fdf_fast, infected0[["fast"]], 
                         e1=e_d,e2=e_d)
  
  pars <- apap_2d(pars, d1, group_seq = T)
  y <- sr(list_modify(pars, e1 = e_t, e2 = e_t), "2d_v3")
  if(rm) return(y)
  
  c("i" = b_any(y, pop, "cumI"), 
    "d" = b_any(y, pop, "D"))
}

df <- expand_grid(d1 = c(d1_general, Inf),
                  e_d = .95,
                  e_t = c(0, .1, .25, .5, .75, .95),
                  model = scenario_par_nms_2v[3]) %>%
  mutate(data = pmap(list(model, d1, e_d, e_t), 
                     function(a,b,c,d) data.frame(value = model_i(a,b,c,d), 
                                                  var = c("i", "d")))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e_d, e_t)  %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 

# Burden relative to no vaccination
g1 <- df %>%
  mutate(e_t = factor(e_t, levels = unique(e_t), labels = as.percent(unique(e_t),0,T))) %>%
  ggplot(aes(x = 1/d1, y = i, color = e_t)) + 
  geom_line(size = 1.2) + 
  # facet_wrap(~model, ncol = 4) +
  scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
  xlab("Vaccination speed [%/day]") + ylab("Infections") +
  scale_color_viridis_d(option = "mako", name = "Efficacy\nagainst\ntransmission (EAT)")
g2 <- df %>%
  mutate(e_t = factor(e_t, levels = unique(e_t), labels = as.percent(unique(e_t),0,T))) %>%
  ggplot(aes(x = 1/d1, y = d, color = e_t)) + 
  geom_line(size = 1.2) + 
  # facet_wrap(~model, ncol = 4) +
  scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
  xlab("Vaccination speed [%/day]") + 
  ylab("Mortality") +
  scale_color_viridis_d(option = "mako", name = "Efficacy\nagainst\ntransmission (EAT)")

# Burden relative to no impact on transmission 
g3 <- df %>%
  ungroup() %>%
  group_by(e_d, model, d1) %>%
  mutate(d = d/max(d)) %>%
  mutate(e_t = factor(e_t, levels = unique(e_t), labels = as.percent(unique(e_t),0,T))) %>%
  ggplot(aes(x = 1/d1, y = d, color = e_t)) + 
  geom_line(size = 1.2) + 
  # facet_wrap(~model, ncol = 4) +
  scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
  xlab("Vaccination speed [%/day]") +
  ylab("Mortality reduction vs EAT 0% case") +
  scale_color_viridis_d(option = "mako", name = "Efficacy\nagainst\ntransmission (EAT)")
                     
burden_vetrans <- df %>%
  gather(var, value, -d1, -model, -e_d, -e_t) %>%
  mutate(var = factor(var, levels = c("i", "d"), 
                      labels = c("Infections", "Deaths"))) %>%
  mutate(delta1 = 1/d1) %>%
  mutate(e_t=factor(e_t,levels=unique(e_t),labels=as.percent(unique(e_t),d=0,perc=T))) %>% 
  ggplot(aes(x = delta1, y = value*100, color = e_t)) + 
  geom_line(size=1) +
  facet_wrap(~var, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = unique(x_ticks), labels = as.percent(unique(1/df$d1)[unique(1/df$d1)!=0])) +
  lightness(scale_color_brewer(palette = "YlOrRd",direction = 1),scalefac(0.95))+
  theme(axis.text.x = element_text(angle = 45, size = 10), legend.position = "bottom") +
  xlab(def_labels$speed) + 
  ylab("Burden (% pop.)")  +
  labs(color="Efficacy against transmission")

# install.packages("patchwork")
library(patchwork)
g1 + g2 + plot_layout(guides = "collect", ncol = 2)
