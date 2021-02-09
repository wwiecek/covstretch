# Base case settings -----
library(tidyverse)
library(ggpubr)
theme_set(theme_minimal(base_size = 10))
# source("R/ode_2doses.R")
source("R/ode_2doses_v2.R")
# source("R/ode_2vaccines.R")
source("R/ode_2vaccines_v2.R")
source("R/helpers.R")
source("R/output-helpers.R")
source("R/scenario-parms.R")
source("R/harm_function.R")
source("R/prioritisation.R")



# Names of compartments ------
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


# For generic and lower efficacy cases -----
df_efficacy_delta_raw <- readRDS(file = "results/df_efficacy_delta_raw.rds")



# Settings for FDF -----
load("results/fdf-deltas.Rdata")
comp_to_display <- c("I", "D", "cumV", "cumI", "P1", "P2")
delay_default <- 18
delay_fdf <- 74
delay_hybrid <- c(rep(delay_fdf, 6), rep(delay_default, 3))

