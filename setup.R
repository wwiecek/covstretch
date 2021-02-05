# Base case settings -----
library(tidyverse)
library(ggpubr)
theme_set(theme_minimal(base_size = 10))
source("R/ode_2doses.R")
source("R/ode_2doses_v2.R")
source("R/ode_2vaccines.R")
source("R/ode_2vaccines_v2.R")
source("R/helpers.R")
source("R/output-helpers.R")
source("R/config.R")
source("R/config-pars.R")
source("R/harm_function.R")
source("R/prioritisation.R")

# Basic inputs for deriving parameters -----
# We use pbc_spread and default_cm from the default data inputs file:
load("data/default_inputs.Rdata")
# Nc <- 9
# Ngroups <- 9
# Ndays <- 1000
r0 <- 2.5
# Use realistic age structure and infection structure
# For now in the HICs
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
ev <- eigen(default_cm)$values[1]





# Names of compartments and scenarios ------
ln <- c(
  "S" = "Susceptible",
  "R" = "Natural immunity (no vaccine)",
  "E" = "Exposed (latent phase)",
  "I" = "Currently infected",
  "P" = "Vaccinated, protected",
  "N" = "Vaccinated, not protected",
  "V1" = "Had vaccine 1, protected",
  "V2" = "Had vaccine 2, protected",
  "P1" = "Had 1 dose, protected",
  "N1" = "Had 1 dose, not protected",
  "P2" = "Had 2 doses, protected",
  "N2" = "Had 2 doses, not protected",
  "cumV1" = "Courses of vaccine 1 used to date",
  "cumV2" = "Courses of vaccine 2 used to date",
  "cumV" = "Total doses used to date",
  "cumI" = "Total new infections to date",
  "D" = "Total deaths to date"
)

scenario_names <- c("Constant risk", "Slow growth", "Fast growth")





# Parameters for each scenario -----
Ndays <- 360
Ngroups <- 9
pop <- hic_pop/sum(hic_pop)
pre_immunity <- c(.5, .5, rep(.2, 7))
pre_immunity_prop <- sum(pre_immunity*pop)

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
  kappa1 = rep(0, Ngroups),
  kappa2 = rep(0, Ngroups),
  phi = rep(0, Ngroups), 
  ta = rep(0, Ngroups),
  e1 = .8,
  e2 = .95,
  pdeath = c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100,
  vrf = 1,
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




# For generic and lower efficacy cases -----
df_efficacy_delta_raw <- readRDS(file = "results/df_efficacy_delta_raw.rds")



# Settings for FDF -----
load("results/fdf-deltas.Rdata")
comp_to_display <- c("I", "D", "cumV", "cumI", "P1", "P2")
delay_default <- 18
delay_fdf <- 74
delay_hybrid <- c(rep(delay_fdf, 6), rep(delay_default, 3))

