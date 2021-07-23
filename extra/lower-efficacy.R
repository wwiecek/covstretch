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
  filter(round(d1,5) %in% le_speeds) %>%
  mutate(delta1 = 1/d1) %>% 
  select(delta1, e, model, i,d,harm) %>%
  gather(var, value, -delta1, -e, -model) %>%
  group_by(model,var) %>%
  mutate(ref = value[e == .95 & delta1 == default_delta_value]) %>%
  filter(e %in% c(seq(.5, .9, .1), .95)) %>%
  ungroup() %>%
  mutate(r = (value/ref)) %>% 
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Fractional dose better by >5%",
                                    
                                    "Comparable (within 5%)", 
                                    "Full dose better by >5%"
                                    ))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(speedup = factor(round(delta1/default_delta_value, 1))) %>%
  # mutate(speedup = factor(as.percent(delta1, 2))) %>%
  mutate(e = factor(e)) %>%
  mutate(speedup = factor(round(delta1/default_delta_value, 4),
                          levels = round(1/c(1,3/4,1/2,1/3,1/4,1/8),4),
                          labels = c("1","3/4","1/2","1/3","1/4","1/8"))) %>%
  filter(var != "Economic harm") %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(hjust = 1,angle = 45)) +
  facet_grid(var~model) + 
  ylab("Efficacy of fractional dose") +
  xlab("Dose fraction") +
  # scale_x_continuous(breaks = 1/d1_general,
                     # labels = as.percent(1/d1_general))
  geom_text(aes(label = format(value,3)), color = "white", size = 2.5)
le2

le_both <- ggpubr::ggarrange(le1 + ggtitle("Burden of infections")+theme(text = element_text(size=9)), 
                  le2 + ggtitle("Optimal policy (relative burden)")+theme(text = element_text(size=9)),
                  ncol = 1, heights = c(.7,1))

le2.plot <- le2


