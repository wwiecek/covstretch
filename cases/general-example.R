d1_general <- c(90, 120, 180, 360, 730, 1460)

df_efficacy_delta <- 
  df_efficacy_delta_raw %>%
  # filter(d1 %in% c(90, 180, 360, 730, 1460, Inf)) %>%
  filter(d1 %in% c(d1_general, Inf)) %>%
  # mutate(benefit_vr = 1 - harm_vr) %>%
  mutate(ref_e = harm[d1 > 1460]) %>%
  mutate(ref_i = i[d1 > 1460]) %>%
  mutate(ref_d = d[d1 > 1460]) %>%
  ungroup() %>%
  mutate(re = 1 - (harm/ref_e)) %>%
  mutate(ri = 1 - (i/ref_i)) %>%
  mutate(rd = 1 - (d/ref_d)) %>%
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
    "Mass vaccination: 360 days" = apap_2v(pars, 360),
    # "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
    "Mass vaccination: 180 days" = apap_2v(pars, 180),
    "Mass vaccination: 90 days" = apap_2v(pars, 90),
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
# ggsave("figures/g1.pdf",g1, width = 8, height=7)




gg1 <- lapply(mlist, function(pars) {
  ll <- list(
    "Mass vaccination: 360 days" = apap_2v(pars, 360),
    # "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
    "Mass vaccination: 180 days" = apap_2v(pars, 180),
    "Mass vaccination: 90 days" = apap_2v(pars, 90),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind()
  as.data.frame(ll[,"I",]) %>%
    mutate(time = 1:nrow(.))
}) %>%
  bind_rows(.id = "scenario") %>%
  gather(var, value, -time, -scenario) %>%
  ggplot(aes(x = time, y = value, color = var)) + geom_line() + facet_wrap(.~scenario, scales = "free") +
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + ylab("fraction infected") +
  scale_color_discrete(name = "") +
  theme(legend.position = "top")
# ggsave("figures/g1_mini.pdf", width = 7, height=2.5)


ll <- list(
  "Mass vaccination: 360 days" = apap_2v(mlist[[1]], 360),
  # "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
  "Mass vaccination: 180 days" = apap_2v(mlist[[1]], 180),
  "Mass vaccination: 90 days" = apap_2v(mlist[[1]], 90),
  "No vaccination" = list_modify(mlist[[1]], delta1 = rep(0, Ngroups))) %>%
  lapply(sr, f = "2v_v2") %>%
  lapply(rescale_rcs, pop, merge=T) %>% 
  abind::abind()
gg2 <- as.data.frame(ll[,"cumV",]) %>%
  mutate(time = 1:nrow(.)) %>%
  gather(var, value, -time) %>%
  ggplot(aes(x = time, y = value, color = var)) + geom_line() + 
  # facet_wrap(.~scenario, scales = "free") +
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + 
  ylab("fraction vaccinated") +
  scale_color_discrete(name = "") +
  theme(legend.position = "top")

gga <- ggarrange(gg2 + ggtitle("Vaccinations"), common.legend = TRUE, 
                 gg1 + ggtitle("Infections"), 
                 widths = c(1,2.5))
ggsave("figures/g1_joint.pdf",gga, width = 7, height=2.5)

# Fig G2A: absolute harm
g2a <- df_efficacy_delta  %>%
  filter(e == .95) %>%
  # filter(d1 < 340) %>%
  filter(d1 %in% d1_general) %>%
  select(i, d) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("i", "d", "harm"), 
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  ggplot(aes(x = d1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 2, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = d1_general[-2]) +
  scale_color_discrete(name = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 9), legend.position = "top") +
  xlab("Vaccination speed, 1/delta") + 
  ylab("Total harm")

# Fig G2B: How many infections and deaths averted with 95% efficacious vaccine -----
g2b <- df_efficacy_delta  %>%
  filter(e == .95) %>%
  # filter(d1 < 340) %>%
  filter(d1 %in% d1_general) %>%
  select(ri, rd) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  ggplot(aes(x = d1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 2, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = d1_general[-2]) +
  scale_color_discrete(name = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 9), legend.position = "top") +
  xlab("Vaccination speed, 1/delta") + 
  ylab("Fraction of harm averted")

g2<-ggarrange(g2a+ggtitle("Burden"), 
          g2b+ggtitle("Reductions"), 
          common.legend = TRUE, ncol = 1, legend = "top")

ggsave("figures/g2.pdf",g2, width = 6, height=6)

# Table G3 -----
df_efficacy_delta %>% 
  filter(e == .95) %>%
  filter(d1 %in% c(d1_general, Inf)) %>% 
  select(d1, model, i, ri, d, rd, harm, re) %>%
  mutate(d = round(d*1e05)) %>%
  # mutate(i = i*1e05) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "harm", "ri", "rd", "re"),
                           labels = c("Infections", "Deaths per 100,000", "Economic harm",
                                      "RI", "RD", "RE"))) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)


