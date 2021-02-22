# Fig LE1: general impact of rate and efficacy on infections ------
le1 <- df_efficacy_delta_raw %>%
  filter(e %in% c(.5, .75, .95)) %>%
  mutate(delta1 = 1/d1) %>%
  # filter(d1 > 50, d1 < 370) %>%
  filter(delta1 >= .0005, delta1 <= .02) %>%
  select(delta1, d, i, harm, model, e) %>%
  gather(var, value, -model, -e, -delta1) %>%
  mutate(e = as.numeric(e)) %>%
  mutate(lab_e = paste("e =", e)) %>%
  group_by(e, model, var) %>% #
  mutate(lab_y = head(value, 1)) %>%
  filter(var == "i") %>%
  ggplot(aes(x = delta1, y = value, group = lab_e)) + 
  geom_line() +
  geom_text(aes(x = max(1/d1_general), y = lab_y, label = lab_e), hjust = 0, size = 2) +
  ylab("Fraction infected in 1 year") + 
  xlab(def_labels$speed) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  facet_wrap(~model, ncol = 3, scales = "free") +
  scale_x_continuous(breaks = 1/d1_general[-5], 
                     labels = as.percent(1/d1_general[-5]), 
                     limits = c(0,1.2*max(1/d1_general[-5])))
  # scale_x_continuous(breaks = seq(0, 360, 120), limits = c(60,500))


# Fig LE2 -----

le2 <- df_efficacy_delta_raw %>%
  filter(d1 %in% le_speeds) %>%
  mutate(delta1 = 1/d1) %>% 
  select(delta1, e, model, i,d,harm) %>%
  gather(var, value, -delta1, -e, -model) %>%
  group_by(model,var) %>%
  mutate(ref = value[e == .95 & delta1 == default_delta_value]) %>%
  filter(e %in% seq(.5, .9, .1)) %>%
  ungroup() %>%
  mutate(r = (value/ref)) %>% 
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Less effective better by 5% or more",
                                    
                                    "Comparable (+-5%)", 
                                    "95% effective better by at least 5%"
                                    ))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(speedup = factor(round(delta1/default_delta_value, 1))) %>%
  # mutate(speedup = factor(as.percent(delta1, 2))) %>%
  mutate(e = factor(e)) %>%
  filter(var != "Economic harm") %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(var~model) + ylab("e2 (efficacy for the less effective vaccine)") + 
  xlab("delta2/delta1 (speed-up factor vs 0.25% base case)") + 
  # scale_x_continuous(breaks = 1/d1_general,
                     # labels = as.percent(1/d1_general))
  geom_text(aes(label = value), color = "white", size = 2.5)


le_both <- ggpubr::ggarrange(le1 + ggtitle("Burden of infections"), 
                  le2 + ggtitle("Optimal policy (relative burden)"),# + theme(legend.direction = "vertical"), 
                  ncol = 1, heights = c(.7,1))


