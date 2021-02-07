par(mfrow = c(1,3), oma = rep(1, 4))


plot_low_eff <- function(pars, label, switch = FALSE){
  w <- sapply(d, function(d1) {
    # No vaccine 2 available earlier:
    case1 <- sr(f="2v_v2", list_modify(
      apap_2v(pars, len = 360, delay = 10+d1), e1 = .95, e2 = 0)) %>% b_any(pop, "cumI")
    # Vaccine 2 available from day 1 until d:
    if(!switch)
      case2 <- sapply(e, function(e1) 
        sr(f="2v_v2", list_modify(
          apap_2v(pars, len = 360), e1 = e1, e2 = 0)) %>% b_any(pop, "cumI"))
    
    else
      case2 <- sapply(e, function(e1) 
        sr(f="2v_v2", list_modify(
          apap_2v(pars, len = 360, switch = d1), e1 = e1, e2 = .95)) %>% b_any(pop, "cumI"))
    
    c(d1, case1, case2)
  })
  
  d <- d/30.5
  plot(w[2, ] ~ d, type = "l", lty = "solid", xlim = c(0, 200/30.5), 
       ylab = "BI: fraction infected in 1 year", 
       xlab = "time until better vaccine available [months]", 
       # ylim = c(min(w[-1,]), max(w[-1,])),
       ylim = c(0, max(w[-1,])),
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
  # No vaccine 2 available earlier:
  case1 <- sr(f="2v_v2", list_modify(
    apap_2v(pars, len = delta, delay = 10+delay), e1 = .95, e2 = 0)) %>% b_any(pop, "cumI")
  # Vaccine 2 available from day 1 until d:
  if(!switch)
    case2 <- 
    sr(f="2v_v2", list_modify(
      apap_2v(pars, len = delta), e1 = efficacy, e2 = 0)) %>% b_any(pop, "cumI")
  
  else
    case2 <- 
    sr(f="2v_v2", list_modify(
      apap_2v(pars, len = delta, switch = delay), e1 = efficacy, e2 = .95)) %>% b_any(pop, "cumI")
  c(case1, case2, case2/case1)
}
# ratio_i_delay("pars_le_cr", 210, 1/180, .95)

ri_df <- expand.grid(model = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                     delay = (0:6), 
                     efficacy = seq(.5, .9, .1), 
                     mvr = c(360),
                     switch = c(TRUE, FALSE)) %>%
  mutate(ri = pmap_dbl(list(model, delay, mvr, efficacy, switch), 
                       function(g,x,y,z,s) {ratio_i_delay(g,30.5*x,y,z,s)[3]}))

gg_df <- ri_df %>%
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
  # filter(model == "Constant risk") %>%
  mutate(delay = factor(delay))
# "Better by less than 50%", "Better by 50% or more"))) %>%
# select(delay, efficacy, az_better) %>% 
# spread(delay, az_better) %>%

gg_df %>% 
  ggplot(aes(x = delay, y = efficacy, fill = az_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "RI for using less effective now vs waiting") +
  theme(legend.position = "bottom") +
  facet_grid(switch ~ model) + ylab("e2 (efficacy for the worse vaccine)") + 
  xlab("time until better vaccine available [months]") +
  geom_text(aes(label = value), color = "white") 

gg_df %>% 
  ggplot(aes(x = delay, y = efficacy, fill = az_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "RI for using less effective now vs waiting") +
  theme(legend.position = "bottom") +
  facet_grid(switch ~ model) + ylab("e2 (efficacy for the worse vaccine)") + 
  xlab("time until better vaccine available [months]")

gg_df %>% 
  filter(switch == "SWITCH = No", model == "Constant risk") %>%
  ggplot(aes(x = delay, y = efficacy, fill = az_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "RI for using less effective now vs waiting") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  xlab("time until better vaccine available [months]") +
  ylab("e2 (efficacy for the worse vaccine)") + 
  geom_text(aes(label = value), color = "white")  +
  ggtitle("Constant risk, no switching")


