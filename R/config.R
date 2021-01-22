library(tidyverse)
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

# pars <- lst(
#   Nc, 
#   Ngroups, Ndays,
#   y0 = matrix(c(1-1e-05, 1e-05, rep(0, Nc - 2)), Nc, Ngroups),
#   q = rep(r0/(5*ev), Ngroups),
#   contacts = default_cm,
#   gamma1 = rep(.2, Ngroups),
#   gamma2 = rep(.2, Ngroups),
#   delta1 = rep(1/365, Ngroups),
#   delta2 = rep(1/365, Ngroups),
#   kappa1 = rep(0, Ngroups),
#   kappa2 = rep(0, Ngroups),
#   phi = rep(1/180, Ngroups), #loss of immunity: 6 months
#   doses_y1 = rep(.4, Ndays),
#   doses_y2 = rep(.2, Ndays),
#   doses_x  = 1:Ndays,
#   e1 = .95,
#   e2 = .62,
#   pdeath = rep(.01, Ngroups),
#   pop_size = hic_pop/sum(hic_pop)
# )

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
