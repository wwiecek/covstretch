# Generate all results

# Sensitivity analyses (global parameters to modify) ------
# We use pbc_spread and default_cm from the default data inputs file:
library(tidyverse)
load("data/default_inputs.Rdata")
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
lic_pop <- pbc_spread[countries["Low-income countries"],] %>% as.numeric()
ifr_hic <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100
ifr_lic <- ifr_hic*(3.2/2)^(5:(-3))

pop <- hic_pop/sum(hic_pop)
default_pdeath <- ifr_hic

use_delta <- TRUE
kappa_default <- 0
delay_default <- 28 - 10
delay_fdf <- 84 - 10
delay_hybrid <- c(rep(delay_fdf, 6), rep(delay_default, 3))


def_labels <- list(
  "speed" = "Fraction vaccinated each day, delta"
)

# default_speeds <- c(seq(60, 360, 10), 450, 540, 630, 730, 1460, Inf)
# main3speeds <- c(360, 180, 90)
default_speeds <- round(100/c(2, rev(seq(.05, 1, .05)), .025, .01, 0), 5) # % per day
fdf_speeds <- rev(round(100/c(.1, .25, .35, .5, .75, 1), 5))
# d1_general <- c(90, 120, 180, 360, 730, 1460)
d1_general <- 100/c(1, .75, .5, .25, .1) # % per day
default_delta_value <- .0025 #for LE scenario
# le_speeds <- round(100/c(.25, .3, .4, .5, .75, 1), 5)
le_speeds <- round(100/c(.25, .3, .35, .4, .5, 1), 5)

source("setup.R")
source("cases/benefits.R")

source("cases/prep-delta-for-fdf.R")
source("cases/first-doses-first-v2.R")

source("cases/prep-results.R")
source("cases/general-example.R")
source("cases/lower-efficacy.R")
source("cases/lower-efficacy-delay.R")
