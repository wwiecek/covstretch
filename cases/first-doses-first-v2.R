fdf_palette <- c("grey20", "#E69F00", "#56B4E9")
# FDF model function (generating main metrics for all results) -----
model_fdf <- function(model, d1, e, policy, comp = c("cumI", "D")) {
  if(model == "pars_le_cr")   pars <- pars_fdf_cr
  if(model == "pars_le_slow")   pars <- pars_fdf_slow
  if(model == "pars_le_fast") pars <- pars_fdf_fast
  d2 <- as.numeric(fdf_deltas[2, fdf_speeds == d1])
  d3 <- as.numeric(fdf_deltas[3, fdf_speeds == d1])
  if(is.infinite(d1)){
    d2 <- Inf; d3 <- Inf
  }
  # if(d2 == 1) browser()
  if(policy == "default")
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    d1, delay_default))
  if(policy == "fdf")
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    d2, delay_fdf))
  if(policy == "hybrid"){
    res <- sr(f = "2d_v2",  apap_2d(list_modify(pars, e1 = e, e2 = .95),
                                    d3, delay_hybrid))
  }
  if(is.infinite(d1) || policy == "no_vaccination")
    res <- sr(f = "2d_v2", 
              list_modify(pars, 
                          e1 = 0,
                          e2 = 0,
                          ta = rep(0, Ngroups),
                          delta1 = rep(0, Ngroups),
                          delta2 = rep(0, Ngroups)))
  main_metrics(res, pop, vat = 91)
}



# Generate a data frame with all values -----

# KEEP AN EYE OUT FOR NUMERICAL ISSUES WITH THE ODEs HERE
df_fdf <- expand.grid(d1 = c(fdf_speeds, Inf),
                      model = c("pars_le_cr", "pars_le_slow", "pars_le_fast"), 
                      e = c(.4, .5, .6, .7, .8, .9, .95), 
                      policy = c("default", "fdf", "hybrid")) %>%
  mutate(data = pmap(list(model, d1, e, policy), 
                     function(x,y,z,a) data.frame(value = model_fdf(x,y,z,a), 
                                                  var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                        labels = c("Constant risk", "Slow growth", "Fast growth"))) %>%
  group_by(d1, model, e, policy)




# Total doses used: illustration -----
lapply(as.list(1:length(fdf_speeds)), function(i) {
  x <- fdf_deltas[,i]
  w <- rescale_and_bind(list(
    "Default (4 weeks)"      = sr(apap_2d(pars_fdf_cr, x[1], delay_default) %>% list_modify(e1 = .8), f = "2d_v2"),
    "FDF (12 weeks)"         = sr(apap_2d(pars_fdf_cr, x[2], delay_fdf) %>% list_modify(e1 = .8), f = "2d_v2"),
    "FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_cr, c(rep(x[3], 6), rep(x[1],3)), delay_hybrid) %>% list_modify(e1 = .8), f = "2d_v2")
    # "FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_cr, 156, delay_hybrid) %>% list_modify(e1 = .8), f = "2d_v2")
  ), pop)
  as.data.frame(w[,"cumV",]) %>% rownames_to_column("time")
}) %>% 
  setNames(fdf_deltas[1,]) %>%
  bind_rows(.id = "d") %>%
  mutate(time = as.numeric(time)) %>%
  mutate(d = as.numeric(d)) %>%
  gather(policy, value, -d, -time) %>%
  filter(d > 1) %>%
  ggplot(aes(x=time, y=value, color=policy)) + geom_line() + facet_wrap(~d) + ylim(0, .75)



# Fig 1: Base case example ------
f1<-rescale_and_bind(list(
  "Default (4 weeks)"        = sr(apap_2d(pars_fdf_slow, fdf_deltas[1,6], delay_default) %>% 
                                    list_modify(e1 = .8), f = "2d_v2"),
  "FDF (12 weeks)"           = sr(apap_2d(pars_fdf_slow, fdf_deltas[2,6], delay_fdf) %>% 
                                    list_modify(e1 = .8), f = "2d_v2"),
  "S-FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_slow, c(rep(fdf_deltas[3,6], 6), rep(fdf_deltas[1,6],3)), delay_hybrid) %>% 
                                    list_modify(e1 = .8), f = "2d_v2")
  # "FDF (12 weeks, hybrid)" = sr(apap_2d(pars_fdf_cr, 156, delay_hybrid) %>% list_modify(e1 = .8), f = "2d_v2")
), pop) %>% 
  plot_rcs(c("P1", "P2", "cumV"), long_names = ln, ncol = 3, start_date = NULL) + 
  ylab("Proportion of population") + theme(legend.position = "top") + 
  scale_color_manual(values = fdf_palette, name = "2nd dose delay policy:") +
  scale_x_continuous(name = "time [days]", breaks = seq(0, 360, 120))

# sr(apap_2d(pars_fdf_fast, 180, delay_default) %>% list_modify(e1 = .8), f = "2d_v2") %>% main_metrics(pop)
# sr(apap_2d(pars_fdf_fast, 145, delay_fdf) %>% list_modify(e1 = .8), f = "2d_v2") %>% main_metrics(pop)
ggsave("figures/fdf1.pdf",f1, width = 6, height=2.5)



df_gg <- df_fdf %>% 
  filter(e %in% c(.5, .95)) %>% 
  filter(d1 > 60, d1 < 540) %>% 
  # select(d1, model, e, policy, d, harm, i) %>%
  select(d1, model, e, policy, d, i) %>%
  gather(var, value, -d1, -model, -e, -policy) %>%
  # mutate(d1 = factor(d1)) %>%
  mutate(policy = factor(policy, levels = c("default", "fdf", "hybrid"), 
                         labels = c("Default (4 wks delay)", 
                                    "FDF (12 wks for all)", 
                                    "S-FDF (4 wks for 60+, 12 for rest)"))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(e = factor(e))

df_gg %>%
  mutate(delta1 = 1/d1) %>%
  ggplot(aes(x = delta1, y = value, lty = e, color = policy)) + 
  geom_line(size=1) +
  # geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(var~model, scales = "free") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45)) +
  scale_color_manual(values = fdf_palette) +
  scale_x_continuous(breaks = 1/fdf_speeds, labels = as.percent(1/fdf_speeds, 1)) +
  scale_linetype_manual(name = "efficacy after 1st dose", values = c("solid", "dashed")) +
  xlab(paste0(def_labels, " (1st dose, default policy)")) + ylab("Burden")

ggsave("figures/sfdf2.pdf", width = 10, height=6)

fdf_burden <- df_gg %>%
  filter(model == "Slow growth", var %in% c("Infections", "Deaths")) %>%
  mutate(delta1 = 1/d1) %>%
  ggplot(aes(x = delta1, y = value, lty = e, color = policy)) + 
  geom_line(size=1) +
  # geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(var~model, scales = "free") +
  theme(legend.position = "top", legend.direction = "vertical", axis.text.x = element_text(angle = 45)) +
  scale_color_manual(values = fdf_palette) +
  scale_x_continuous(breaks = 1/fdf_speeds, labels = as.percent(1/fdf_speeds, 1)) +
  scale_linetype_manual(name = "efficacy after 1st dose", values = c("solid", "dashed")) +
  xlab(paste0(def_labels, " (1st dose, default policy)")) + ylab("") +
  guides(linetype = FALSE, color = FALSE)

# ggsave("figures/fdf_burden.pdf", width = 4, height=2)



# Reductions as a function of vaccination speed -----
fdf_reductions <- df_fdf %>% 
  filter(e %in% c(.5, .95)) %>% 
  # select(d1, model, e, policy, d, harm, i) %>%
  select(d1, model, e, policy, d, i) %>%
  group_by(model, e, policy) %>%
  mutate(d = 1- d/d[is.infinite(d1)], i = 1 - i/i[is.infinite(d1)]) %>%
  ungroup() %>%
  filter(d1 > 60, d1 <= 1000) %>%
  gather(var, value, -d1, -model, -e, -policy) %>%
  # mutate(d1 = factor(d1)) %>%
  mutate(policy = factor(policy, levels = c("default", "fdf", "hybrid"), 
                         labels = c("Default (4 wks delay)", 
                                    "FDF (12 wks for all)", 
                                    "S-FDF (4 wks for 60+, 12 for rest)"))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(e = factor(e)) %>%
  filter(model == "Slow growth", var %in% c("Infections", "Deaths")) %>%
  mutate(delta1 = 1/d1) %>%
  ggplot(aes(x = delta1, y = value, lty = e, color = policy)) + 
  geom_line(size=1) +
  # geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(var~model, scales = "free") +
  theme(legend.position = "top", 
        legend.box = "vertical",
        axis.text.x = element_text(angle = 45)) +
  scale_color_manual(values = fdf_palette) +
  scale_x_continuous(breaks = 1/fdf_speeds, labels = as.percent(1/fdf_speeds, 1)) +
  scale_linetype_manual(name = "efficacy after 1st dose", values = c("solid", "dashed")) +
  xlab(paste0(def_labels, " (1st dose, default policy)")) + ylab("Fraction of harm averted")+
  guides(linetype=FALSE)

# ggsave("figures/fdf_reductions.pdf", width = 6, height=3)



# optimal solution -----

fig2 <- df_fdf %>% 
  filter(d1 > 60, d1 <= 1000) %>%
  filter(e >= .5) %>%
  select(d1, model, e, policy, d, harm, i) %>%
  gather(var, value, -d1, -model, -e, -policy) %>%
  group_by(d1, model, e, var) %>% 
  summarise(value_m = min(value[policy != "default"])/value[policy == "default"], 
            policy = policy[which.min(value)]
  ) %>%
  mutate(policy = factor(policy, levels = c("default", "fdf", "hybrid"), 
                         labels = c("Default (4 wks delay)", 
                                    "FDF (12 wks for all)", 
                                    "S-FDF (4 wks for 60+, 12 for rest)"))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  # mutate(efficacy = ifelse(e == .8, "Efficacy after 1 dose = 80%", "Efficacy after 1 dose = 50%")) %>%
  mutate(delta1 = factor(1/d1,
                          labels = as.percent(unique(1/d1, 1)))) %>%
  mutate(e = factor(e)) %>%
  mutate(value_m = round(value_m, 2)) %>%
  filter(var != "Economic harm") %>%
  ggplot(aes(x = delta1, y = e, fill = policy)) + 
  geom_tile() +
  # scale_fill_viridis_d(name="") +  
  # scale_fill_manual(values = c("grey20", "grey40", "grey60"),
  #                   name = "") +
  scale_fill_manual(values = fdf_palette, name = "") +
  theme(legend.position = "bottom") +
  facet_grid(var~model) + 
  ylab("e1 (efficacy following 1st dose)") +
  theme(axis.text.x = element_text(angle = 45)) +
  xlab(paste0(def_labels, " (1st dose, default policy)"))

fig2
# ggsave("figures/fdf4.pdf", fig4, width = 6, height=6)

fig2s <- fig2 + geom_text(aes(label = value_m), size = 2, color = "white")
fig2s
# ggsave("figures/sfdf4.pdf", fig4s, width = 9, height=6)

# One big figure for FDF section -----

w <- ggpubr::ggarrange(fdf_reductions + ggtitle("Reductions"), 
                       fig2s + ggtitle("Optimal policy (relative burden)"), ncol = 1, heights = c(2, 2),
                       common.legend = TRUE)
ggsave("figures/fdf2.pdf",w , width = 6.5, height=8)




