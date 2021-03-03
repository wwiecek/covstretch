pars_delay <- list_modify(pars_le_slow,
                          y0 = y0_gen(13, 9, pre_immunity, 3e-05),
                          q = rep(3/(5*ev), Ngroups))

# Check that the peak occurs at 120 days
sr(pars_delay) %>% rescale_rcs(pop, TRUE) -> w
cat("Peak for the delay scenario is ", which.max(w[,"I",]))

model_delay <- function(d1, delay) {
  set0 <- c(rep(1,4), 0, rep(1,7), 0) #resets cumulative deaths and infections to 0
  pars <- apap_2v(pars_delay, d1)
  if(delay > 0)
    pars <- list_modify(pars, 
                        y0 = (sr(pars_delay)[as.character(delay*30), ,])*set0)
  y <- sr(list_modify(pars, e1 = .95), "2v_v2")
  main_metrics(y, pop)
}

df <- expand_grid(d1 = default_speeds, 
                  delay = 0:9) %>%
  mutate(data = pmap(list(d1, delay), 
                     function(x,y) data.frame(value = model_delay(x,y), 
                                              var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value)

# Burdens
select(df, d1, delay, i) %>%
  spread(delay, i) %>%
  filter(d1 %in% c(d1_general, Inf)) %>%
  mutate(d1 = 100/d1)
# Reductions
select(df, d1, delay, i) %>%
  group_by(delay) %>%
  mutate(i = 1 - i/max(i)) %>%
  ungroup() %>%
  spread(delay, i) %>%
  filter(d1 %in% c(d1_general, Inf)) %>%
  mutate(d1 = 100/d1)


inf_rates <- data.frame(delay = (1:nrow(w))/30, imax = w[,"I",]/max(w[,"I",]), imin = 0)

fig_delay <- select(df, d1, delay, i) %>%
  group_by(delay) %>%
  mutate(i = 1 - i/max(i)) %>%
  ungroup() %>%
  filter(d1 %in% c(d1_general[-1])) %>%
  mutate(d1 = factor(d1, levels = c(d1_general), labels = c(as.percent(1/d1_general)))) %>%
  ggplot() + 
  geom_ribbon(data = inf_rates, aes(x=delay, ymin=imin, ymax=imax), alpha = .10) +
  geom_line(aes(x=delay, y=i, group=d1, color=d1)) +
  theme(legend.position = "top") +
  scale_color_discrete(name = "Vaccinated per day") +
  xlab("Vaccination start date (months until infection peak)") + 
  ylab("Fraction of infections averted") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6, 7, 8), labels = -(0:8) + 4, limits = c(0, 9))
