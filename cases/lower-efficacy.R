source("cases/prep-results.R")

# Fig LE1: general impact of rate and efficacy on infections ------
df_efficacy_delta_raw %>%
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
  ggplot(aes(x = d1, y = value, color = lab_e)) + 
  # geom_point() +
  geom_line(size = 1.1) +
  geom_text(aes(x = 370, y = lab_y, label = lab_e), hjust = 0) +
  ylab("fraction infected in 1 year") + 
  xlab("length of mass vaccination program, 1/delta [days]") +
  # geom_hline(yintercept = 0, lty = "dashed") +
  facet_wrap(~model, ncol = 3, scales = "free") +
  xlim(60, 500) +
  scale_color_discrete(guide = NULL)



# Fig LE2 -----

df_efficacy_delta_raw %>%
  filter(d1 %in% c(60, 90, 180, 240, 300, 360)) %>%
  # filter(d1 %in% seq(30, 360, 30)) %>%
  # filter(d1 <= 360) %>%
  select(d1, e, model, i) %>%
  group_by(model) %>%
  mutate(ref = i[e == .95 & d1 == 360]) %>%
  ungroup() %>%
  filter(d1 < 360) %>%
  filter(e != .95) %>%
  filter(e %in% seq(.5, .9, .1)) %>%
  mutate(ri = 1 - (i/ref)) %>% 
  mutate(le_better = cut(ri, c(-Inf, -.05, .05, Inf), 
                         labels = c("95% effective better by at least 5%",
                                    "Comparable (+-5%)", 
                                    "Less effective better by 5% or more"))) %>%
  mutate(value = round(ri, 2)) %>%
  mutate(speedup = factor(360/d1)) %>%
  mutate(e = factor(e)) %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "Infection reduction:") +
  theme(legend.position = "bottom") +
  facet_grid(~model) + ylab("e2 (efficacy for the worse vaccine)") + 
  xlab("Speed-up factor (delta2/delta1 = 360*delta2)") +
  geom_text(aes(label = value), color = "white")  
# ggtitle("When to choose less efficacious vaccine?")






# Fig LE3: Lower efficacy vaccine now vs higher efficacy later ----
d <- seq(0, 24*7, 7)
e <- c(.4, .5, .6, .7, .8)

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
  filter(switch == "SWITCH = No", model == "Constant risk") %>%
  ggplot(aes(x = delay, y = efficacy, fill = az_better)) + geom_tile() +
  scale_fill_manual(values = c("grey20", "grey40", "grey60"), 
                    name = "RI for using less effective now vs waiting") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  xlab("time until better vaccine available [months]") +
  ylab("e2 (efficacy for the worse vaccine)") + 
  geom_text(aes(label = value), color = "white")  +
  ggtitle("Constant risk, no switching")
