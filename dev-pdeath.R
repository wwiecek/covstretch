pars_fdf_slow %>%
  apap_2d(180, 18) %>%
  sr(f = "2d_v2") %>% 
  plot_rcs(c("R", "S", "cumV", "D"))

list_modify(pars_fdf_slow,
            Nc = 17,
            y0 = y0_gen(17, 9, pre_immunity, 5e-03, E = 2, I = 5, R = 8),
            pdeath = NULL,
            pd0 = 0*ifr_hic,
            pd1 = ifr_hic,
            pd2 = ifr_hic) %>%
  apap_2d(180, 18) %>%
  sr(f = "2d_v3") %>% 
  plot_rcs(c("R", "S", "cumV", "D"))


expand_2d_v3 <- function(pars, x, e1, e2) {
  pd <- pars$pdeath
  list_modify(pars,
              Nc = 17,
              y0 = y0_gen(17, 9, pre_immunity, x, E = 2, I = 5, R = 8),
              pdeath = NULL,
              pd0 = pd,
              pd1 = (1-e1)*pd,
              pd2 = (1-e2)*pd)
}

model_i <- function(model, d1, e_d, e_t, rm = FALSE) {
  if(model == "pars_linear")    
    pars <- expand_2d_v3(pars_fdf_linear, infected0[["decreasing"]], 
                         e1=.5*e_d,e2=e_d)
  if(model == "pars_le_slow")   
    pars <- expand_2d_v3(pars_fdf_slow, infected0[["slow"]], 
                         e1=.5*e_d,e2=e_d)
  if(model == "pars_le_fast")   
    pars <- expand_2d_v3(pars_fdf_fast, infected0[["fast"]], 
                         e1=.5*e_d,e2=e_d)
  
  pars <- apap_2d(pars, d1)
  y <- sr(list_modify(pars, e1 = .5*e_t, e2 = e_t), "2d_v3")
  if(rm) return(y)
  
  c("i" = b_any(y, pop, "cumI"), 
    "d" = b_any(y, pop, "D"))
}

df <- expand_grid(d1 = c(d1_general, Inf),
                  e_d = .95,
                  e_t = c(0, .25, .5, .95),
                  model = scenario_par_nms_2v) %>%
  mutate(data = pmap(list(model, d1, e_d, e_t), 
                     function(a,b,c,d) data.frame(value = model_i(a,b,c,d), 
                                                  var = c("i", "d")))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e_d, e_t)  %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 

# Reductions relative to no vaccination
df %>%
  mutate(d = 1 - d/max(d)) %>%
  mutate(e_t = factor(e_t)) %>%
  ggplot(aes(x = d1, y = d, color = e_t)) + geom_line() + facet_wrap(~model, ncol = 4)
#
# Reductions relative to no impact on transmission
df %>%
  ungroup() %>%
  group_by(e_d, model, d1) %>%
  mutate(d = d/max(d)) %>%
  mutate(e_t = factor(e_t)) %>%
  ggplot(aes(x = d1, y = d, color = e_t)) + geom_line() + facet_wrap(~model, ncol = 4)
