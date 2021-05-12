##Fast growth----

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

df_delay_vac <- expand_grid(d1 = default_speeds, 
                  delay = 0:8) %>%
  mutate(data = pmap(list(d1, delay), 
                     function(x,y) data.frame(value = model_delay(x,y), 
                                              var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value)

# Burdens
select(df_delay_vac, d1, delay, i) %>%
  spread(delay, i) %>%
  filter(d1 %in% c(d1_general, Inf)) %>%
  mutate(d1 = 100/d1)
# Reductions
select(df_delay_vac, d1, delay, i) %>%
  group_by(delay) %>%
  mutate(i = 1 - i/max(i)) %>%
  ungroup() %>%
  spread(delay, i) %>%
  filter(d1 %in% c(d1_general, Inf)) %>%
  mutate(d1 = 100/d1)


inf_rates <- data.frame(delay = (1:nrow(w))/30, imax = w[,"I",]*100/max(w[,"I",]), imin = 0)

fig_delay <- select(df_delay_vac, d1, delay, i) %>%
  group_by(delay) %>%
  mutate(i = 1 - i/max(i)) %>%
  ungroup() %>%
  filter(d1 %in% c(d1_general)) %>%
  mutate(d1 = factor(d1, levels = c(rev(d1_general)), labels = c(as.percent(rev(1/d1_general),perc=TRUE)))) %>%
  ggplot() + 
  geom_ribbon(data = inf_rates, aes(x=delay, ymin=imin, ymax=imax), alpha = .10) +
  geom_line(aes(x=delay, y=i*100, group=d1, color=d1),size=1) +
  theme(legend.position = "top") +
  lightness(scale_color_brewer(palette = "YlOrRd",direction = 1,),scalefac(0.8)) +
  # xlab("Vaccination start (months after infection peak)") + 
  ylab("% infections averted") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6, 7), labels = (0:7) - 4, limits = c(0, 7)) +
  labs(color="Vaccinated per day")
fig_delay

daily_deaths <- w[1:dim(w)[1],"D",]-c(0,w[1:dim(w)[1]-1,"D",])
death_rates <- data.frame(delay = (1:nrow(w))/30, dmax = daily_deaths*100/max(daily_deaths), dmin = 0)
fig_delay.deaths <- select(df_delay_vac, d1, delay, d) %>%
  group_by(delay) %>%
  mutate(d = 1 - d/max(d)) %>%
  ungroup() %>%
  filter(d1 %in% c(d1_general)) %>%
  mutate(d1 = factor(d1, levels = c(rev(d1_general)), labels = c(as.percent(rev(1/d1_general),perc=TRUE)))) %>%
  ggplot() + 
  geom_ribbon(data = inf_rates, aes(x=delay, ymin=imin, ymax=imax), alpha = .10) +
  geom_line(aes(x=delay, y=d*100, group=d1, color=d1),size=1) +
  theme(legend.position = "top") +
  lightness(scale_color_brewer(palette = "YlOrRd",direction = 1,),scalefac(0.8)) +
  # xlab("Vaccination start (months after infection peak)") + 
  ylab("% deaths averted") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6, 7), labels = (0:7) - 4, limits = c(0, 7)) +
  labs(color="Vaccinated per day")

fig_delay.both <- ggpubr::ggarrange(
  fig_delay + ggtitle("Infections")+ rremove("xlab") + theme(text = element_text(size=10),plot.title = element_text(size = 10)),
  fig_delay.deaths + ggtitle("Deaths")+ rremove("xlab")+ theme(text = element_text(size=10),plot.title = element_text(size = 10)),
  labels = NULL,
  ncol = 2,
  common.legend = TRUE)
fig_delay.both <- annotate_figure(fig_delay.both, bottom = textGrob("Vaccination start (months after infection peak)", gp = gpar(cex = 0.8)))
##Slow growth----
pars_delay <- list_modify(pars_le_slow,
                          y0 = y0_gen(13, 9, pre_immunity, 4.75e-03))

# Check that the peak occurs at 120 days
sr(pars_delay) %>% rescale_rcs(pop, TRUE) -> w
cat("Peak for the delay scenario is ", which.max(w[,"I",]))

df_delay_vac.slow <- expand_grid(d1 = default_speeds, 
                            delay = 0:12) %>%
  mutate(data = pmap(list(d1, delay), 
                     function(x,y) data.frame(value = model_delay(x,y), 
                                              var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value)

inf_rates.slow <- data.frame(delay = (1:nrow(w))/30, imax = w[,"I",]*100/max(w[,"I",]), imin = 0)

fig_delay.slow <- select(df_delay_vac.slow, d1, delay, i) %>%
  group_by(delay) %>%
  mutate(i = 1 - i/max(i)) %>%
  ungroup() %>%
  filter(d1 %in% c(d1_general)) %>%
  mutate(d1 = factor(d1, levels = c(rev(d1_general)), labels = c(as.percent(rev(1/d1_general),perc=TRUE)))) %>%
  ggplot() + 
  geom_ribbon(data = inf_rates.slow, aes(x=delay, ymin=imin, ymax=imax), alpha = .10) +
  geom_line(aes(x=delay, y=i*100, group=d1, color=d1),size=1) +
  theme(legend.position = "top") +
  lightness(scale_color_brewer(palette = "YlOrRd",direction = 1),scalefac(0.8)) +
  xlab("Vaccination start (months after infection peak)") + 
  ylab("% infections averted") +
  scale_x_continuous(breaks = 0:12, labels = (0:12) - 4, limits = c(0, 12)) +
  labs(color="Vaccinated per day")
fig_delay.slow
