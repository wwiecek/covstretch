
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
    "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
    # "Mass vaccination (no prioritsation), delta = 1/360" = list_modify(pars, delta1 = rep(1/360, Ngroups)),
    "Mass vaccination (no prioritsation), delta = 1/60" = list_modify(pars, delta1 = rep(1/60, Ngroups)),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(sr) %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind() %>% 
    plot_rcs(c("S", "I", "cumV1", "D"), ncol = 5, long_names = ln2,
             end_date = as.Date("01-01-2021", format="%d-%m-%Y") + 300) + ylab("") +
    ggtitle(names(mlist)[i])
})

ggarrange(plotlist=gglist, common.legend = TRUE, ncol = 1, legend = "top")

# sr(list_modify(pars_le_cr, delta1 = rep(0, Ngroups))) %>% b_any(pop, "D")
# sr(list_modify(pars_le_slow, delta1 = rep(0, Ngroups))) %>% b_any(pop, "D")
# sr(list_modify(pars_le_fast, delta1 = rep(0, Ngroups))) %>% b_any(pop, "D")


# Fig G2: How many infections and deaths averted with 95% efficacious vaccine -----
d1_general <- c(60, 90, 180, 360, 730, 1460)
df_efficacy_delta  %>%
  filter(e == .95) %>%
  # filter(d1 < 340) %>%
  filter(d1 %in% d1_general) %>%
  select(ri, rd, re) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("RI (infections)", "RD (deaths)", "RE (economic harm)"))) %>%
  ggplot(aes(x = d1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 3, fill = "white") +
  facet_wrap(~key, scales = "free", ncol = 2) +
  scale_x_continuous(breaks = d1_general[-1]) +
  scale_color_discrete(name = "scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 12), legend.position = "top") +
  xlab("average time to vaccination, 1/ delta1 [days]") + ylab("proportion of harm averted")



# Table G3 -----
df_efficacy_delta %>% 
  filter(e == .95) %>%
  filter(d1 %in% c(d1_general, Inf)) %>% select(d1, model, i, ri, benefit_r, v1, tthi) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)


