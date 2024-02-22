#-------------------------------------------------------------------------------
# This file defines general parameters regardless of scenarios/models
#-------------------------------------------------------------------------------

# Load pbc_spread and default contact matrix from the default data inputs file.
# This loads four items in the environment: 
# - pbc_spread 
# - default_cm 
# - countries
# - cases_csv_clean

library(tidyverse)
load("setup/default_input.Rdata")

# Baseline efficacy
default_e <- 0.95

# Main FDF assumptions
delay_default <- 28 - 10
delay_fdf <- 84 - 10
delay_hybrid_k <- sapply(1:8, function(k){
  c(rep(delay_fdf, k), rep(delay_default, 9-k))
})
colnames(delay_hybrid_k) <- 1:8
all_k <- TRUE
default_group_seq <- TRUE

# Demographics (for comparing HIC vs LIC)
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
lic_pop <- pbc_spread[countries["Low-income countries"],] %>% as.numeric()
ifr_hic <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100
ifr_lic <- ifr_hic*(3.2/2)^(5:(-3))
pop <- hic_pop/sum(hic_pop)
default_pdeath <- ifr_hic
# use_delta <- TRUE

# Case with losing immunity
kappa_default <- 0
# Case with lower supply (25% vs 100%)
default_supply_ceiling <- 1

def_labels <- list(
  "speed" = "Percentage of pop. vaccinated daily"
)

# default_speeds <- c(seq(60, 360, 10), 450, 540, 630, 730, 1460, Inf)
# main3speeds <- c(360, 180, 90)
default_speeds <- round(100/c(2, rev(seq(.05, 1, .05)), .025, .01, 0), 5) # % per day
default_speeds <- unique(c(default_speeds,
                           round(1/c(c(0.001,0.0025,0.005,0.0075,0.01,0.02),
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*1.2,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*1.4,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*1.6,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*2,
                                     c(0.001,0.0025,0.005,0.0075,0.01,0.02)*4),5)))
fdf_speeds <- rev(round(100/c(.1, .25, .5, .75, 1, 2), 5))
# d1_general <- c(90, 120, 180, 360, 730, 1460)
d1_general <- 100/c(2, 1, .75, .5, .25, .1) # % per day
default_delta_value <- .0025 #for LE scenario
# le_speeds <- round(100/c(.25, .3, .35, .4, .5, 1, 2), 5)
le_speeds <- round(1/0.0025*c(1,3/4,1/2,1/4,1/3,1/8), 5)
default_speeds <- unique(c(default_speeds,le_speeds))


as.percent <- function(x, d=2, perc=FALSE){
  if (perc) paste0(format(round(100*x, d), nsmall = d),'%')
  else format(round(100*x, d), nsmall = 2)
}

# Long names of compartments ------
ln <- c(
  "S" = "Susceptible",
  "R" = "Natural immunity (no vaccine)",
  "RV" = "Recovered following vaccination",
  "E" = "Exposed (latent phase)",
  "E0" = "Exposed (not vaccinated)",
  "E1" = "Exposed (had 1 dose)",
  "E2" = "Exposed (had 2 doses)",
  "I"  = "Infectious",
  "I0" = "Infectious (not vaccinated)",
  "I1" = "Infectious (had 1 dose)",
  "I2" = "Infectious (had 2 doses)",
  "P"  = "Vaccinated, protected",
  "N"  = "Vaccinated, susceptible",
  "V1" = "Had vaccine 1, protected",
  "V2" = "Had vaccine 2, protected",
  "P1" = "Had one dose, protected",
  "N1" = "Had one dose, susceptible",
  "P2" = "Had two doses, protected",
  "N2" = "Had two doses, susceptible",
  "cumV1" = "Vaccinated",
  "cumV2" = "Courses of vaccine 2 used to date",
  "cumV" = "Cumulative doses",
  "cumI" = "Total new infections to date",
  "D" = "Total deaths to date"
)

# For generic and lower efficacy cases -----
df_efficacy_delta_raw <- readRDS(file = "setup/df_efficacy_delta_raw.rds")

# Settings for FDF -----
comp_to_display <- c("I", "D", "cumV", "cumI", "P1", "P2")

y0_gen <- function(Nc, Ngroups, 
                   pre_immunity = rep(0, Ngroups),
                   # v = rep(0, Ngroups), 
                   ii = 5e-03, 
                   S = 1, E = 2, I = 3, R = 4){
  y0_default <- matrix(0, Nc, Ngroups)
  y0_default[S,] <- 1-pre_immunity
  y0_default[E,] <- (y0_default[1,]*ii)/2
  y0_default[I,] <- (y0_default[1,]*ii)/2
  # subtract initial infections
  y0_default[S,] <- y0_default[1,] - y0_default[2,]
  y0_default[R,] <- pre_immunity
  y0_default
}

