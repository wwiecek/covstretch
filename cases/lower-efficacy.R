source("cases/prep-results.R")

# Fig LE1: general impact of rate and efficacy on infections ------

par(mfrow = c(1,3))
e <- c(.6, .8, .95)
x <- seq(10, 360, 10)


y0 <- sapply(e, function(e) sapply(x, function(d1) sr(list_modify(pars_le_cr, e1 = e,
                                                                  delta1 = rep(1/d1, Ngroups))) %>% 
                                     b_any(pop, "cumI")))
y1 <- sapply(e, function(e) sapply(x, function(d1) sr(list_modify(pars_le_slow, e1 = e,
                                                                  delta1 = rep(1/d1, Ngroups))) %>% 
                                     b_any(pop, "cumI")))
y2 <- sapply(e, function(e) sapply(x, function(d1) sr(list_modify(pars_le_fast, e1 = e,
                                                                  delta1 = rep(1/d1, Ngroups))) %>% 
                                     b_any(pop, "cumI")))

colnames(y0) <- colnames(y1) <- colnames(y2) <- e
rownames(y0) <- rownames(y1) <- rownames(y2) <- x

list(y0, y1, y2) %>% lapply(as.data.frame) %>% lapply(rownames_to_column, "t") %>%
  lapply(function(x) gather(x, e, value, -t)) %>%
  setNames(scenario_names) %>%
  bind_rows(.id = "model") %>%
  mutate(model = factor(model, levels = scenario_names)) %>%
  mutate(e = as.numeric(e)) %>%
  mutate(lab_e = paste("e =", e)) %>%
  mutate(t = as.numeric(t)) %>%
  group_by(e, model) %>% #
  mutate(lab_y = tail(value, 1)) %>%
  ggplot(aes(x = t, y = value, group = e, color = lab_e)) + 
  geom_line() +
  geom_text(aes(x = 400, y = lab_y, label = lab_e)) +
  ylab("BI (fraction infected in 1 year)") + 
  xlab("average time to vaccination 1/delta1 [days]") +
  # geom_hline(yintercept = 0, lty = "dashed") +
  facet_wrap(~model, ncol = 3, scales = "free_y") + xlim(0, 450) +
  scale_color_discrete(guide = NULL)





# Fig LE2 -----

df_efficacy_delta %>%
  filter(d1 %in% c(60, 90, 180, 240, 360)) %>%
  # filter(d1 %in% seq(30, 360, 30)) %>%
  # filter(d1 <= 360) %>%
  select(d1, e, model, i) %>%
  group_by(model) %>%
  mutate(ref = i[e == .95 & d1 == 360]) %>%
  ungroup() %>%
  filter(e != .95) %>%
  filter(e %in% seq(.5, .9, .1)) %>%
  mutate(ri = 1 - (i/ref)) %>% 
  mutate(le_better = cut(ri, c(-Inf, -.05, .05, Inf), 
                         labels = c("Better by 5% or more", 
                                    "Comparable (+-5%)", 
                                    "Worse by at least 5%"))) %>%
  mutate(value = round(ri, 2)) %>%
  mutate(speedup = factor(360/d1)) %>%
  mutate(e = factor(e)) %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "RI for less effective vs 95% effective vaccine:") +
  theme(legend.position = "bottom") +
  facet_grid(~model) + ylab("e2 (efficacy for the worse vaccine)") + 
  xlab("Speed-up factor (delta2/delta1 = 360*delta2)") +
  geom_text(aes(label = value), color = "white")  
  # ggtitle("When to choose less efficacious vaccine?")






# Fig LE3: Lower efficacy vaccine now vs higher efficacy later ----
d <- seq(0, 24*7, 7)
e <- c(0, .25, .4, .5, .6, .7)

par(mfrow = c(1,3), oma = rep(1, 4))

plot_low_eff <- function(pars, label, switch = FALSE){
  w <- sapply(d, function(d1) {
    # No vaccine 2 available earlier:
    case1 <- sr(list_modify(pars, 
                            e1 = 0,
                            e2 = .95,
                            delta1 = 0*pars_le_fast$delta1,
                            delta2 = pars_le_fast$delta1,
                            ts1 = rep(d1, 9),
                            ta1 = rep( 0, 9),
                            ta2 = rep(d1, 9))) %>% b_any(pop, "cumI")
    # Vaccine 2 available from day 1 until d:
    if(!switch)
      case2 <- sapply(e, function(e1) sr(list_modify(pars, 
                                                     e1 = 0,
                                                     e2 = e1,
                                                     delta1 = rep(0, Ngroups),
                                                     delta2 = pars_le_fast$delta1,
                                                     ts1 = rep(d1, 9),
                                                     ta1 = rep( 0, 9),
                                                     ta2 = rep( 0, 9))) %>% b_any(pop, "cumI"))
    else
      case2 <- sapply(e, function(e1) sr(list_modify(pars, 
                                                     e1 = e1,
                                                     e2 = .95,
                                                     delta1 = rep(1/360, Ngroups),
                                                     delta2 = rep(1/360, Ngroups),
                                                     ts1 = rep(d1, 9),
                                                     ta1 = rep( 0, 9),
                                                     ta2 = rep(d1, 9))) %>% b_any(pop, "cumI"))
    c(d1, case1, case2)
  })
  
  d <- d/30.5
  plot(w[2, ] ~ d, type = "l", lty = "solid", xlim = c(0, 200/30.5), 
       ylab = "BI: fraction infected in 1 year", 
       xlab = "time until better vaccine available [months]", ylim = c(min(w[-1,]), max(w[-1,])),
       main = label)
  for(i in 3:nrow(w)) {
    lines(w[i, ] ~ d, lty = "dotted")
    text(x = 180/30.5, y = max(w[i, ]), paste0("e = ", e[i-2]))
  }
  
}

plot_low_eff(pars_le_cr, "Constant risk")
plot_low_eff(pars_le_slow, "Slow growth")
plot_low_eff(pars_le_fast, "Fast growth")
title(outer = TRUE, "SWITCH to more effective once available = No")

plot_low_eff(pars_le_cr, "Constant risk", TRUE)
plot_low_eff(pars_le_slow, "Slow growth", TRUE)
plot_low_eff(pars_le_fast, "Fast growth", TRUE)
title(outer = TRUE, "SWITCH = Yes")





# Fig LE4: Lattice comparing vaccination rates, vaccine efficacy, length of delay -----
# switch = does a better vaccine replace the worse one once available
ratio_i_delay <- function(model, delay, delta, efficacy, switch = TRUE) {
  # Naive get()
  if(model == "pars_le_cr")   pars <- pars_le_cr
  if(model == "pars_le_slow")   pars <- pars_le_slow
  if(model == "pars_le_fast") pars <- pars_le_fast
  # No vaccine 2 available earlier:
  case1 <- sr(list_modify(pars, 
                          e1 = 0,
                          e2 = .95,
                          delta1 = rep(0, Ngroups),
                          delta2 = rep(delta, Ngroups),
                          ts1 = rep(delay, 9),
                          ta1 = rep( 0, 9),
                          ta2 = rep(delay, 9))) %>% b_any(pop, "cumI")
  if(switch)
    # Vaccine 2 available from day 1 until d:
    case2 <- sr(list_modify(pars, 
                            e1 = efficacy,
                            e2 = .95,
                            delta1 = rep(delta, Ngroups),
                            delta2 = rep(delta, Ngroups),
                            ts1 = rep(delay, 9),
                            ta1 = rep( 0, 9),
                            ta2 = rep(delay, 9))) %>% b_any(pop, "cumI")
  else
    # Vaccine 2 available from day 1 forever:
    case2 <- sr(list_modify(pars, 
                            e1 = 0,
                            e2 = efficacy,
                            delta1 = rep(0, Ngroups),
                            delta2 = rep(delta, Ngroups),
                            ts1 = rep(delay, 9),
                            ta1 = rep( 0, 9),
                            ta2 = rep( 0, 9))) %>% b_any(pop, "cumI")
  c(case1, case2, case2/case1)
}
# ratio_i_delay("pars_le_cr", 210, 1/180, .95)

ri_df <- expand.grid(model = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                     delay = (0:6), 
                     efficacy = seq(.45, .9, .05), 
                     mvr = c(180),
                     switch = c(TRUE, FALSE)) %>%
  mutate(ri = pmap_dbl(list(model, delay, mvr, efficacy, switch), 
                       function(g,x,y,z,s) {ratio_i_delay(g,30.5*x,1/y,z,s)[3]}))

ri_df %>%
  mutate(ri = 1 - ri) %>%
  mutate(az_better = cut(ri, c(-Inf, -.05, .05, Inf), 
                         labels = c("Worse by at least 5%", 
                                    "Comparable (+-5%)", 
                                    "Better by at least 5%"))) %>%
  mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                        labels = c("Constant risk", "Slow growth", "Fast growth"))) %>%
  mutate(value = round(ri, 2)) %>%
  mutate(mvr = paste("1/delta =", mvr)) %>%
  mutate(switch = ifelse(switch, "SWITCH = Yes", "SWITCH = No")) %>%
  filter(model == "Constant risk") %>%
  mutate(delay = factor(delay)) %>%
  # "Better by less than 50%", "Better by 50% or more"))) %>%
  # select(delay, efficacy, az_better) %>% 
  # spread(delay, az_better) %>%
  ggplot(aes(x = delay, y = efficacy, fill = az_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "RI for using less effective now vs waiting") +
  theme(legend.position = "bottom") +
  facet_grid(~switch) + ylab("e1 (efficacy for the worse vaccine)") + 
  xlab("time until better vaccine available [months]") +
  geom_text(aes(label = value), color = "white")  +
  ggtitle("Constant risk model")
