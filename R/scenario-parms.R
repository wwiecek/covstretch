# Basic inputs for deriving parameters -----
# Nc <- 9
# Ngroups <- 9
# Ndays <- 1000
r0 <- 2.5
# Use realistic age structure and infection structure
# For now in the HICs


# Parameters for each scenario -----
Ndays <- 360
Ngroups <- 9
pre_immunity <- c(.5, .5, rep(.2, 7))
pre_immunity_prop <- sum(pre_immunity*pop)

# R0 adjustment:
ev <- eigen(default_cm)$values[1]
# Pre-immunity adjustment:
ev_pi <- eigen(default_cm*(1-pre_immunity))$values[1]
r0 <- function(r) r/(5*ev_pi)
ev/ev_pi


# Two doses model -----
pars_fdf_slow <- lst(
  Nc = 13, 
  Ngroups, 
  Ndays,
  y0 = y0_gen(13, 9, pre_immunity, 5e-03),
  q = rep(r0(1.1), Ngroups),
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
  vstop = rep(.8, Ngroups), #around 80% vaccinated we slow down
  #this may be modified in apap() function(s)
  constantrisk = 0
)
pars_fdf_fast <- list_modify(pars_fdf_slow,
                             y0 = y0_gen(13, 9, pre_immunity, 1e-03),
                             q = rep(r0(2), Ngroups))
pars_fdf_linear <- list_modify(pars_fdf_slow,
                               y0 = y0_gen(13, 9, pre_immunity, 1e-02),
                               q = rep(r0(1), Ngroups))
pars_fdf_cr <- list_modify(pars_fdf_slow,
                           y0 = y0_gen(13, 9, pre_immunity, .1/30.5),
                           q = rep(0, Ngroups),
                           constantrisk = .01/30.5)
pars_fdf_end <- list_modify(pars_fdf_fast,
                            y0 = y0_gen(13, 9, rep(.5, Ngroups), 1e-03))
# c(rep(12,1), 0) resets cumI (compartment 13) to 0:
set0 <- c(rep(1,4), 0, rep(1,7), 0)
pars_fdf_late <- list_modify(pars_fdf_fast, 
                             y0 = (sr(pars_fdf_fast, f = "2d_v2")["120", ,])*set0)

# Two vaccines model -----
pars_le_slow <- list_modify(pars_fdf_slow,
                            e1 = 0.95, e2 = 0, 
                            ta1 = rep(0, Ngroups), 
                            ta2 = rep(0, Ngroups), 
                            tmore1 = rep(Inf, Ngroups),
                            tmore2 = rep(Inf, Ngroups),
                            ts1 = rep(Ndays, Ngroups))
pars_le_slow$ta  <- NULL
pars_le_fast <- list_modify(pars_le_slow,
                            y0 = y0_gen(13, 9, pre_immunity, 1e-03),
                            q = rep(r0(2), Ngroups))
pars_le_cr <- list_modify(pars_le_slow,
                          y0 = y0_gen(13, 9, pre_immunity, .1/30.5),
                          q = rep(0, Ngroups),
                          constantrisk = .01/30.5)
pars_le_linear <- list_modify(pars_le_slow,
                              y0 = y0_gen(13, 9, pre_immunity, 1e-02),
                              q = rep(r0(1), Ngroups))

pars_le_late <- list_modify(pars_le_fast, 
                            y0 = (sr(pars_le_fast)["120", ,])*set0)

scenario_par_nms_2v <- c("pars_linear", "pars_le_slow", "pars_le_fast")#, "pars_le_late", "pars_le_cr")
scenario_nms_2v <- c("Linear growth", "Slow growth", "Fast growth")#, "Declining risk", "Constant risk")
scenario_list_2v <- lst(
  # "Constant risk of infection" = pars_le_cr,
  "Linear growth (R0 = 1)" = pars_le_linear,
  "Slow growth (R0 = 1.1)" = pars_le_slow,
  "Fast growth (R0 = 2)" = pars_le_fast,
  # "Declining risk (R0 = 3, after peak)" = pars_le_late
) %>%
  setNames(scenario_nms_2v)

