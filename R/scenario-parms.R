# Basic inputs for deriving parameters -----
# Nc <- 9
# Ngroups <- 9
# Ndays <- 1000
r0 <- 2.5
# Use realistic age structure and infection structure
# For now in the HICs
ev <- eigen(default_cm)$values[1]



# Parameters for each scenario -----
Ndays <- 360
Ngroups <- 9
pre_immunity <- c(.5, .5, rep(.2, 7))
pre_immunity_prop <- sum(pre_immunity*pop)


scenario_names <- c("Constant risk", "Slow growth", "Fast growth")

# Two doses model -----
pars_fdf_slow <- lst(
  Nc = 13, 
  Ngroups, 
  Ndays,
  y0 = y0_gen(13, 9, pre_immunity, 5e-03),
  q = rep(1.6/(5*ev), Ngroups),
  contacts = default_cm,
  gamma1 = rep(.2, Ngroups),
  gamma2 = rep(.2, Ngroups), #duration of infectious period
  delta1 = rep(0, Ngroups),
  delta2 = rep(0, Ngroups),
  kappa1 = rep(kappa_default, Ngroups),
  kappa2 = rep(kappa_default, Ngroups),
  phi = rep(0, Ngroups), 
  ta = rep(0, Ngroups),
  e1 = .8,
  e2 = .95,
  pdeath = default_pdeath,
  vrf = 1,
  vstop = .8, #around 80% vaccinated we slow down
  constantrisk = 0
)
pars_fdf_fast <- list_modify(pars_fdf_slow,
                             y0 = y0_gen(13, 9, pre_immunity, 1e-03),
                             q = rep(3/(5*ev), Ngroups))
pars_fdf_cr <- list_modify(pars_fdf_slow,
                           y0 = y0_gen(13, 9, pre_immunity, .1/30.5),
                           q = rep(0, Ngroups),
                           constantrisk = .01/30.5)

# Two vaccines model
pars_le_slow <- list_modify(pars_fdf_slow,
                            e1 = 0.95, e2 = 0, 
                            ta1 = rep(0, Ngroups), 
                            ta2 = rep(0, Ngroups), 
                            ts1 = rep(Ndays, Ngroups))
pars_le_slow$ta  <- NULL
pars_le_fast <- list_modify(pars_le_slow,
                            y0 = y0_gen(13, 9, pre_immunity, 1e-03),
                            q = rep(3/(5*ev), Ngroups))
pars_le_cr <- list_modify(pars_le_slow,
                          y0 = y0_gen(13, 9, pre_immunity, .1/30.5),
                          q = rep(0, Ngroups),
                          constantrisk = .01/30.5)
