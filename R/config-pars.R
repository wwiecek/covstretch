Ndays <- 360
Ngroups <- 9
pop <- hic_pop/sum(hic_pop)
pre_immunity <- c(.5, .5, rep(.2, 7))
pre_immunity_prop <- sum(pre_immunity*pop)


# Two doses model -----
pars_fdf_slow <- lst(
  Nc = 11, 
  Ngroups, Ndays,
  y0 = y0_gen(11, 9, pre_immunity, 1e-02),
  q = rep(1.5/(5*ev), Ngroups),
  contacts = default_cm,
  gamma1 = rep(.2, Ngroups),
  gamma2 = rep(.2, Ngroups), #duration of infectious period
  delta1 = rep(1/360, Ngroups),
  delta2 = rep(0, Ngroups),
  kappa1 = rep(0, Ngroups),
  kappa2 = rep(0, Ngroups),
  phi = rep(0, Ngroups), 
  doses_y = rep(2, Ndays),
  doses_x  = 1:Ndays,
  ta = rep(0, Ngroups),
  e1 = .8,
  e2 = .95,
  pdeath = c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100,
  pop_size = hic_pop/sum(hic_pop),
  v_recovered_flag = 1,
  constantrisk = 0
)

pars_fdf_fast <- list_modify(pars_fdf_slow,
                                 y0 = y0_gen(11, 9, pre_immunity, 1e-03),
                                 q = rep(3/(5*ev), Ngroups))
pars_fdf_cr <- list_modify(pars_fdf_slow,
                               y0 = y0_gen(11, 9, pre_immunity, 0),
                               q = rep(0, Ngroups),
                               constantrisk = .01/30.5)

# Two vaccines model -----
pars_le_slow <- list_modify(pars_fdf_slow,
                          e1 = 0.95, e2 = 0, 
                          Nc = 12, 
                          y0 = y0_gen(12, 9, pre_immunity, 1e-02),
                          doses_y1 = rep(1, Ndays), 
                          doses_y2 = rep(1, Ndays),
                          ta1 = rep(0, Ngroups), 
                          ta2 = rep(0, Ngroups), 
                          ts1 = rep(Ndays, Ngroups))
pars_le_slow$doses_y <- pars_le_slow$ta  <- NULL
pars_le_fast <- list_modify(pars_le_slow,
                            y0 = y0_gen(12, 9, pre_immunity, 1e-03),
                            q = rep(3/(5*ev), Ngroups))
pars_le_cr <- list_modify(pars_le_slow,
                          y0 = y0_gen(12, 9, pre_immunity, 0),
                          q = rep(0, Ngroups),
                          constantrisk = .01/30.5)
