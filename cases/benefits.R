

# Sequencing vaccination (prevaccination scenario) ------


pars_d <- lst(
  Nc = 9, 
  Ngroups, 
  Ndays,
  y0 = matrix(c(1-1e-05, 1e-05, rep(0, 7)), 9, Ngroups),
  q = rep(3/(5*ev), Ngroups),
  contacts = default_cm,
  gamma1 = rep(.2, Ngroups),
  gamma2 = rep(.2, Ngroups),
  delta1 = rep(0, Ngroups),
  delta2 = rep(0, Ngroups),
  kappa1 = rep(0, Ngroups),
  kappa2 = rep(0, Ngroups),
  phi = rep(0, Ngroups), 
  doses_y1 = rep(0, Ndays),
  doses_y2 = rep(0, Ndays),
  doses_x  = 1:Ndays,
  e1 = 0,
  e2 = 0,
  pdeath = c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100,
  # pdeath = rep(.01, Ngroups),
  pop_size = hic_pop/sum(hic_pop)
)

vac_top_p <- function(p, pop) {
  pop <- pop/sum(pop)
  w <- rev(pop)
  ptemp <- p
  vrev <- vector(length = 9)
  for(j in 1:9){
    if(ptemp > 0)
      vrev[j] <- min(ptemp, w[j])
    ptemp <- ptemp - w[j]
  }
  rev(vrev)/pop
}

benefits_curve <- function(e = .95, pop = hic_pop/sum(hic_pop)){
  ii <- 1e-04
  sapply(
    seq(0,1,length=21),
    function(p) {
      v <- vac_top_p(p, pop)
      # y0_v <- matrix(c(1-1e-05, 1e-05, rep(0, 10)), 12, Ngroups)
      # y0_v[6,] <- v
      # y0_v[1,] <- (1-v) - (1-v)*ii
      # y0_v[2,] <- (1-v)*ii
      y0_v <- y0_gen(12, Ngroups, pre_immunity = pre_immunity*v + (1-pre_immunity)*e*v)
      sr <- sr(list_modify(pars_le_cr, 
                           y0 = y0_v))
      c(p = p,
        bd = bd(sr, pop),
        bi = b_any(sr, pop, "cumI"))
    }) %>% t() %>% as.data.frame()
}

bf <- rbind(
  benefits_curve(.8) %>% mutate(scenario = "80% efficacy"),
  benefits_curve(.6) %>% mutate(scenario = "60% efficacy"),
  benefits_curve(.4) %>% mutate(scenario = "40% efficacy"),
  benefits_curve(.95) %>% mutate(scenario = "95% efficacy")) %>%
  mutate(bi = 1-(bi/max(bi))) %>%
  mutate(bd = 1-(bd/max(bd)))

theme_set(theme_minimal(base_size = 14))
bf %>% 
  setNames(c("p", "Deaths averted", "Infections averted", "scenario")) %>%
  gather(var, value, -p, -scenario) %>%
  ggplot(aes(x = p, y = value, group = scenario, lty = scenario)) + 
  geom_line() + facet_wrap(~var, scales = "free") + xlab("Fraction vaccinated") + ylab("Fraction averted")

