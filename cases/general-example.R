
df_efficacy_delta <- 
  df_efficacy_delta_raw %>%
  # filter(d1 %in% c(90, 180, 360, 730, 1460, Inf)) %>%
  filter(round(d1,4) %in% c(round(d1_general,4), Inf)) %>%
  # mutate(benefit_vr = 1 - harm_vr) %>%
  mutate(ref_e = harm[d1 > 1460]) %>%
  mutate(ref_i = i[d1 > 1460]) %>%
  mutate(ref_d = d[d1 > 1460]) %>%
  ungroup() %>%
  mutate(re = 1-(harm/ref_e)) %>%
  mutate(ri = 1-(i/ref_i)) %>%
  mutate(diffi = (i-ref_i)) %>%
  mutate(rd = 1-(d/ref_d)) %>%
  mutate(diffd = (d-ref_d)) %>%
  group_by(d1, model, e)




# Figure G1: general ilustration of the model and benefits of vaccination -----
mlist <- list(
  "Constant risk of infection" = pars_le_cr,
  "Base case (R0 = 1.5)" = pars_le_slow,
  "Fast growth (R0 = 3)" = pars_le_fast) %>%
  setNames(scenario_names)

ln2 <- ln
ln2[["cumV1"]] <- "Courses of vaccine used"
gglist <- lapply(as.list(1:3), function(i) {
  pars <- mlist[[i]]
  list(
    "Vaccinate 3/1000 per day" = apap_2v(pars, 100/.3),
    "Vaccinate 5/1000 per day" = apap_2v(pars, 100/.5),
    "Vaccinate 10/1000 per day"  = apap_2v(pars, 100/1),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind() %>% 
    # plot_rcs(c("S", "I", "cumV1", "D"), ncol = 5, long_names = ln2,
    plot_rcs(c("S", "I", "cumV1", "D"), ncol = 5, long_names = ln2,
             start_date = NULL
             # end_date = as.Date("01-01-2021", format="%d-%m-%Y") + 300
    ) + ylab("") +
    ggtitle(names(mlist)[i]) + 
    xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120))
})

g1<-ggarrange(plotlist=gglist, common.legend = TRUE, ncol = 1, legend = "top")

gg1.df <- lapply(mlist, function(pars) {
  ll <- list(
    "Vaccinate 0.3% per day" = apap_2v(pars, 100/.3),
    "Vaccinate 0.5% per day" = apap_2v(pars, 100/.5),
    "Vaccinate 1% per day"  = apap_2v(pars, 100/1),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind()
  as.data.frame(ll[,"I",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) 

gg1 <- gg1.df %>%
  ggplot(aes(x = time, y = value, color = var)) + geom_line() + facet_wrap(.~scenario, scales = "free") +
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + ylab("fraction infected") +
  scale_color_discrete(name = "") +
  theme(legend.position = "top")

pars <- mlist[[1]]
ll <- list(
  "Vaccinate 0.3% per day" = apap_2v(pars, 100/.3),
  "Vaccinate 0.5% per day" = apap_2v(pars, 100/.5),
  "Vaccinate 1% per day"  = apap_2v(pars, 100/1),
  "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
  lapply(sr, f = "2v_v2") %>%
  lapply(rescale_rcs, pop, merge=T) %>% 
  abind::abind()
gg2.df <- as.data.frame(ll[,"cumV",]) %>%
  mutate(time = 1:nrow(.)) %>%
  gather(var, value, -time) 

gg2 <- gg2.df %>%
  ggplot(aes(x = time, y = value, color = var)) + geom_line() + 
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + 
  ylab("fraction currently vaccinated (I)") +
  scale_color_discrete(name = "") +
  theme(legend.position = "top")
gg2
g1_joint <- ggarrange(gg2 + ggtitle("Vaccinations") + theme(legend.spacing.x = unit(0.1, 'in'), text = element_text(size=7),legend.key.size = unit(0.4, "cm")), common.legend = TRUE, 
                 gg1 + ggtitle("Infections") + theme(text = element_text(size=7)), 
                 widths = c(1,2.5))


# Age-specific dynamics -----
sgg_age <- sr(f="2v_v2", apap_2v(mlist[[2]], 360)) %>% plot_rcs(c("I", "S", "cumV1"), ncol = 3) + 
  ylab("Proportion of age group")+ theme(legend.spacing.y = unit(0.1, 'in'), legend.text = element_text(size = 5), legend.key.size = unit(0.5, "cm"))



# library(scales)

# Fig G2A: absolute harm
g2a <- df_efficacy_delta  %>%
  filter(e == .95) %>%
  filter(round(d1,4) %in% round(d1_general,4)) %>%
  select(i, d) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("i", "d", "harm"), 
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(delta1 = 1/d1) %>%
  ggplot(aes(x = delta1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 2, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
  scale_color_discrete(name = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 9), legend.position = "top") +
  xlab(def_labels$speed) + 
  ylab("Total harm") 


# Fig G2B: How many infections and deaths averted with 95% efficacious vaccine -----
g2b <- df_efficacy_delta  %>%
  filter(e == .95) %>%
  filter(round(d1,4) %in% round(d1_general,4)) %>%
  select(ri, rd) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(delta1 = 1/d1) %>%
  ggplot(aes(x = delta1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 2, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
  scale_color_discrete(name = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 9), legend.position = "top") +
  xlab(def_labels$speed) + 
  ylab("Fraction of harm averted")

fig_g2<-ggarrange(g2a+ggtitle("Burden"), 
                  g2b+ggtitle("Reductions"), 
                  common.legend = TRUE, ncol = 1, legend = "top")


# Table G3 -----
df_efficacy_delta %>% 
  filter(e == .95) %>%
  filter(d1 %in% c(d1_general, Inf)) %>% 
  mutate(d1 = as.percent(1/d1)) %>%
  select(d1, model, i, ri, d, rd, diffi, diffd) %>%
  mutate(d = round(d*1e05)) %>%
  mutate(diffd = round(diffd*1e05)) %>%
  # mutate(i = i*1e05) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "harm", "ri", "rd", "re", "diffd", "diffi"),
                           labels = c("Infections", "Deaths per 100,000", "Economic harm",
                                      "RI", "RD", "RE", "Difference in deaths", "Difference in infections"))) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)

#Plot - Increasing Relative Speed
df.vax_rate <- data.frame()
for (vax_rate in c(0.001,0.0025,0.005,0.0075,0.01,0.02)){
  df.vax_rate <- rbind(df.vax_rate, df_efficacy_delta_raw %>%
                         mutate(delta1 = round(1/d1,5)) %>% 
                         select(delta1, e, model, i,d,harm) %>%
                         gather(var, value, -delta1, -e, -model) %>%
                         mutate(ref_rate = vax_rate) %>%
                         group_by(model,var) %>%
                         mutate(ref = value[e == .95 & delta1 == vax_rate]) %>%
                         filter(delta1 %in% c(vax_rate*c(1.5,2,3,4)) & e == .95) %>%
                         ungroup() %>%
                         mutate(r = (value/ref)) %>% 
                         mutate(vax_rate_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                                                      labels = c("Faster vaccination better by 5% or more",
                                                                 
                                                                 "Comparable (+-5%)", 
                                                                 "Slower vaccination better by at least 5%"
                                                      ))) %>%
                         mutate(var = factor(var, levels = c("i", "d", "harm"),
                                             labels = c("Infections", "Deaths", "Economic harm"))) %>%
                         mutate(value = round(r, 2)) %>%
                         mutate(speedup = factor(round(delta1/vax_rate, 1))) %>%
                         mutate(ref_rate = factor(ref_rate)) %>%
                         filter(var != "Economic harm"))
}

fig.vax_rate <- df.vax_rate %>% 
  ggplot(aes(x = speedup, y = ref_rate)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(var~model) + ylab("baseline vaccination rate (delta1)") + 
  xlab("delta2/delta1 (speed-up factor vs base case)") + 
  scale_y_discrete(breaks = c(0.001,0.0025,0.005,0.0075,0.01,0.02),labels = as.percent(c(0.001,0.0025,0.005,0.0075,0.01,0.02))) +
  geom_text(aes(label = value), color = "white", size = 2)

