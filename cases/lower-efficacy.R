source("setup.R")

# Fig LE1: general impact of rate and efficacy on infections ------
le1 <- df_efficacy_delta_raw %>%
  filter(e %in% c(.5, .75, .95)) %>%
  filter(d1 > 50, d1 < 370) %>%
  select(d1, d, i, harm, model, e) %>%
  gather(var, value, -model, -e, -d1) %>%
  mutate(e = as.numeric(e)) %>%
  mutate(lab_e = paste("e =", e)) %>%
  # mutate(t = as.numeric(t)) %>%
  group_by(e, model, var) %>% #
  mutate(lab_y = tail(value, 1)) %>%
  filter(var == "i") %>%
  ggplot(aes(x = d1, y = value, group = lab_e)) + 
  # geom_point() +
  geom_line() +
  geom_text(aes(x = 370, y = lab_y, label = lab_e), hjust = 0, size = 2) +
  ylab("Fraction infected in 1 year") + 
  xlab("Vaccination speed, 1/delta") +
  # geom_hline(yintercept = 0, lty = "dashed") +
  facet_wrap(~model, ncol = 3, scales = "free") +
  # xlim(60, 500) +
  # scale_color_discrete(guide = NULL) + 
  scale_x_continuous(breaks = seq(0, 360, 120), limits = c(60,500))
# ggsave("figures/le1.pdf", le1, width = 19, height=12)


# Fig LE2 -----

le2 <- df_efficacy_delta_raw %>%
  filter(d1 %in% c(60, 90, 180, 240, 300, 360)) %>%
  # filter(d1 %in% seq(30, 360, 30)) %>%
  # filter(d1 <= 360) %>%
  select(d1, e, model, i,d,harm) %>%
  gather(var, value, -d1, -e, -model) %>%
  group_by(model,var) %>%
  mutate(ref = value[e == .95 & d1 == 360]) %>%
  ungroup() %>%
  # filter(d1 < 360) %>%
  filter(e %in% seq(.5, .9, .1)) %>%
  mutate(r = (value/ref)) %>% 
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Less effective better by 5% or more",
                                    
                                    "Comparable (+-5%)", 
                                    "95% effective better by at least 5%"
                                    ))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(speedup = factor(360/d1)) %>%
  mutate(e = factor(e)) %>%
  filter(var != "Economic harm") %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom") +
  facet_grid(var~model) + ylab("e2 (efficacy for the less effective vaccine)") + 
  xlab("Speed-up factor (delta2/delta1 = 360*delta2)") +
  geom_text(aes(label = value), color = "white", size = 2.5)  

# ggsave("figures/le2.pdf", le2,width = 6, height=4)


le12 <- ggpubr::ggarrange(le1 + ggtitle("Burden"), 
                  le2 + ggtitle("Optimal policy (relative burden)"),# + theme(legend.direction = "vertical"), 
                  ncol = 1, heights = c(.7,1))
ggsave("figures/le_both.pdf", le12, width = 6, height = 7)

