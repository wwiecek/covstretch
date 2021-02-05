
par(mfrow = c(1,3))
e <- c(.6, .8, .95)
x <- seq(10, 360, 10)


y0 <- sapply(e, function(e) sapply(x, function(d1) sr(list_modify(pars_le_cr, e1 = e,
                                                                  delta1 = rep(1/d1, Ngroups))) %>% 
                                     benefit(pop, r = TRUE)))
y1 <- sapply(e, function(e) sapply(x, function(d1) sr(list_modify(pars_le_slow, e1 = e,
                                                                  delta1 = rep(1/d1, Ngroups))) %>% 
                                     benefit(pop, r = TRUE)))
y2 <- sapply(e, function(e) sapply(x, function(d1) sr(list_modify(pars_le_fast, e1 = e,
                                                                  delta1 = rep(1/d1, Ngroups))) %>% 
                                     benefit(pop, r = TRUE)))

colnames(y0) <- colnames(y1) <- colnames(y2) <- e
rownames(y0) <- rownames(y1) <- rownames(y2) <- x

fdf2<-list(y0, y1, y2) %>% lapply(as.data.frame) %>% lapply(rownames_to_column, "t") %>%
  lapply(function(x) gather(x, e, value, -t)) %>%
  setNames(scenario_names) %>%
  bind_rows(.id = "model") %>%
  mutate(model = factor(model, levels = scenario_names)) %>%
  mutate(e = as.numeric(e)) %>%
  mutate(lab_e = paste("Efficacy =", e)) %>%
  mutate(t = as.numeric(t)) %>%
  group_by(e, model) %>% #
  mutate(lab_y = tail(value, 1)) %>%
  mutate(value = 1 - value/365) %>%
  ggplot(aes(x = t, y = value, group = interaction(e, model), color = lab_e, lty = model)) + 
  geom_line() +
  # geom_text(aes(x = 400, y = lab_y, label = lab_e)) +
  ylab("BE (yearly economic harm)") + 
  xlab("Vaccination speed, 1/delta")  +
  theme(legend.position = "top", legend.direction = "vertical") +
  # geom_hline(yintercept = 0, lty = "dashed") +
  # facet_wrap(~model, ncol = 3) + xlim(0, 450) +
  scale_color_discrete(name = "") + 
  scale_linetype_discrete(name = "") 

ggsave("figures/fdf2.pdf", fdf2,width = 19, height=12)
