
# Table with all results -----
df_fdf_ratio <- df_fdf %>% select(d, i, harm) %>%
  pivot_wider(values_from = c("d", "i", "harm"), names_from = "policy") %>%
  # mutate(re = 1 - (harm_vr_fdf/harm_vr_default)) %>%
  # mutate(re = harm_vr_fdf - harm_vr_default) %>%
  # mutate(rd = d_fdf - d_default) %>%
  # mutate(ri = i_fdf - i_default) %>%
  mutate(re = (harm_fdf / harm_default)) %>%
  mutate(rd = (   d_fdf / d_default)) %>%
  mutate(ri = (   i_fdf / i_default)) %>%
  # mutate(ri_fn = 1 - (i_fdf/i_no_vaccination)) %>%
  # mutate(ri_dn = 1 - (i_default/i_no_vaccination)) %>%
  group_by(d1, model, e)


filter(df_fdf_ratio, e %in% c(.6, .8)) %>%
  filter(d1 %in% c(d1_default, Inf)) %>% 
  # select(d1, model, ri_fd, ri_fn, ri_dn, i_fdf, i_default) %>%
  # select(d1, model, i_fdf, i_default) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(e, variable, model) %>% 
  ungroup() 
# select(-e, -variable)

# Dependence between speed of vaccination and magniutude of benefits -----
df_fdf_ratio  %>%
  filter(e %in% c(.6, .8)) %>%
  filter(d1 > 60, d1 < 1000) %>%
  select(ri, rd, re) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("RRI (infections)", "RRD (deaths)", "RRE (economic harm)"))) %>%
  mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 60%")) %>%
  ggplot(aes(x = d1, y = value, group = interaction(model, efficacy), lty = efficacy, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap( ~ key, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = c(90, 180, 360, 730, 1460)) +
  scale_color_discrete(name = "scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 14), legend.position = "top") +
  xlab("average time to first dose, 1/ delta1 [days]") + 
  ylab("relative risk FDF vs default policy")




# Figure FDF4: optimal choice FDF vs default depending on speed and efficacy ------
FDF4<-df_fdf_ratio %>%
  filter(d1 > 1, d1 < Inf) %>%
  # mutate(rd = re) %>%
  select(d1, e, model, ri, re, rd) %>%
  gather(key, value, -d1, -e, -model) %>%
  mutate(rd = value) %>%
  mutate(fdf_better = cut(rd, c(-Inf, .95, 1.05, Inf), 
                          labels = c("FDF better by 5% or more", 
                                     "Comparable (+-5%)", 
                                     "FDF worse by at least 5%"))) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"),
                      labels = c("RRI (infections)", "RRD (deaths)", "RRE (economic harm)"))) %>%
  mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 50%")) %>%
  mutate(value = round(value, 2)) %>%
  mutate(speedup = factor(d1)) %>%
  mutate(e = factor(e)) %>%
  ggplot(aes(x = speedup, y = e, fill = fdf_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "Reduction in infections RI:") +
  theme(legend.position = "bottom") +
  facet_grid(key~model) + ylab("e1 (efficacy following 1st dose)") + 
  xlab("average wait until the first dose under default policy, 1/delta1 [days]") +
  geom_text(aes(label = value), color = "white")  

ggsave("figures/fdf4.pdf",FDF4, width = 19, height=12)




df_fdf %>% filter(model == "Constant risk", e == .95) %>% select(d1, model, policy, v1) %>% spread(policy, v1)


