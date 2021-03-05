
delay_palette <- c("Switch V2->V1" = "#D55E00", 
                   "No switching (V2 only)" = "#009E73", 
                   "No early vaccine (V1 only)" = "grey20")

# Fig LE3: Lower efficacy vaccine now vs higher efficacy later ----
d <- seq(0, 24*7, 7)
e <- c(.4, .5, .6, .7, .8)


model_delay <- function(model, delay, e, scenario, rm = FALSE) {
  pars <- grab_2v_parms(model)
  
  if(scenario == "v1only")
    y <- sr(f="2v_v2", list_modify(
      apap_2v(pars, len = 1/default_delta_value, delay = 10+delay), e1 = .95, e2 = 0))
  
  if(scenario == "noswitch")
    y <- sr(f="2v_v2", list_modify(
      apap_2v(pars, len = 1/default_delta_value), e1 = e, e2 = 0))
  
  if(scenario == "switch")
    y <-   sr(f="2v_v2", list_modify(
      apap_2v(pars, len = 1/default_delta_value, switch = delay), e1 = e, e2 = .95))
  
  if(rm) return(y)
  main_metrics(y, pop)
}


df_delay <- expand.grid(delay = 15.25*(0:12),
            e = seq(.5, .9, .1),
            scenario = c("switch", "noswitch", "v1only"),
            model = scenario_par_nms_2v) %>%
  mutate(data = pmap(list(model, delay, e, scenario), function(x,y,z,s) 
    data.frame(value = model_delay(x,y,z,s), var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  group_by(model, e, delay)  %>%

  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) %>%
  mutate(scenario = factor(scenario, levels = c("switch", "noswitch", "v1only"),
                           labels = c("Switch V2->V1", "No switching (V2 only)", "No early vaccine (V1 only)")))


gg_delay <- df_delay %>%
  select(-harm, -tt50, -v1) %>%
  gather(var, value, -delay, -e, -scenario, -model) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm")))

delay_burden <- 
  gg_delay %>%
  filter(e %in% c(.5, .8)) %>%
  mutate(e=factor(e)) %>%
  filter(var == "Infections") %>%
  ggplot(aes(x = delay, y = value, linetype = e, color = scenario)) +
  geom_line() +
  facet_wrap(var ~ model, scales = "free") +
  scale_x_continuous(breaks = c(0, 60, 120, 180)) +
  scale_color_manual(values = delay_palette, name = "") +
  guides(linetype=F) +
  theme(legend.position = "top") +
  ylab("Fraction infected in 1 year") +
  xlab("time until vaccine 1 available [days]")

delay_optimal <- 
gg_delay %>% 
  filter(delay %in% (30.5*(0:6))) %>%
  spread(scenario, value) %>%
  mutate(r = `No switching (V2 only)`/`No early vaccine (V1 only)`) %>%
  select(delay, e, model, var, r) %>%
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Less effective better by 5% or more",
                                    "Comparable (+-5%)", 
                                    "95% effective better by at least 5%"
                         ))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(delay = factor(delay/30.5)) %>%
  ggplot(aes(x = delay, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom") +
  facet_grid(var~model) + ylab("e2 (efficacy for the less effective vaccine)") + 
  xlab("Months until 95% effective vaccine available") +
  geom_text(aes(label = value), color = "white", size = 2.5)  

delay_switch <- 
gg_delay %>% 
  filter(delay %in% (30.5*(0:6))) %>%
  spread(scenario, value) %>%
  mutate(r = `Switch V2->V1`/`No early vaccine (V1 only)`) %>%
  select(delay, e, model, var, r) %>%
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Less effective better by 5% or more",
                                    "Comparable (+-5%)", 
                                    "95% effective better by at least 5%"
                         ))) %>%
  mutate(value = round(r, 2)) %>%
  filter(delay > 0) %>%
  mutate(delay = factor(delay/30.5)) %>%
  ggplot(aes(x = delay, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom") +
  facet_grid(var~model) + ylab("e2 (efficacy for the less effective vaccine)") + 
  xlab("Months until 95% effective vaccine available") +
  geom_text(aes(label = value), color = "white", size = 2.5)  +
  ggtitle("Optimal policy if switching allowed, 0.25% vaccinated/day")



delay_both <- ggpubr::ggarrange(delay_burden + ggtitle("Burden of infections, 0.25% vaccinated/day"), 
                                delay_optimal + ggtitle("Optimal policy, no switching allowed, 0.25% vaccinated/day")+theme(text = element_text(size=9)),
                                ncol = 1, heights = c(.7,1))
