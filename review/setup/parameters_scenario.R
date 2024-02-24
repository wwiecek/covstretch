# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# This file contains everything needed to set up and adjust parameters to be 
# passed to the epidemiological model simulation. It is organized as follows:
# 
# 1. Initial compartment values
# 2. Universal parameters
# 3. Two-dose model parameters
# 4. Two-vaccine model parameters
# 5. Selection function of parameters: grab_2v_parms
# 6. Adjustment function of parameters: apap_2v
# 7. Unknown
# 
# 
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------





# ------------------------------------------------------------------------------
# 1. Initial compartment values
# ------------------------------------------------------------------------------
# First load the default input data which contains:
# - default_cm: the contact matrix C
# - pbc_spread: population by age groups for 249 countries
# - cases_csv_clean: N/A
# - countries: correspondence between country and country code
load("setup/default_input.Rdata")

# Extract high-income and low-income demographics
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
lic_pop <- pbc_spread[countries["Low-income countries"],] %>% as.numeric()
pop <- hic_pop / sum(hic_pop) # Fractions of total population by age groups

# Proportion of natrual immunity by age groups
pre_immunity <- c(.5, .5, rep(.2, 7))
pre_immunity_prop <- sum(pre_immunity * pop)

# Proportion of initial infection, same across all age groups
infected0 <- c("fast" = 1e-03,
               "slow" = 5e-03,
               "decreasing" = 1e-02)

# Function for generating initial compartment values
# - Nc: number of compartments
# - Ngroups: number of age groups
# - pre_immunity
# - ii: essentially infected0 in its actual applications
# - S/E/I/R: order of these four compartments in the output matrix
y0_gen <- function(Nc, Ngroups, 
                   pre_immunity = rep(0, Ngroups),
                   ii = 5e-03, 
                   S = 1, E = 2, I = 3, R = 4){
  y0_default <- matrix(0, Nc, Ngroups)
  y0_default[S,] <- 1 - pre_immunity
  y0_default[E,] <- (y0_default[1,]*ii)/2
  y0_default[I,] <- (y0_default[1,]*ii)/2 # Exposed == Infected initially
  y0_default[S,] <- y0_default[1,] - y0_default[2,] # Subtract initial infections
  y0_default[R,] <- pre_immunity
  y0_default
}
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# 2. Universal parameters
#-------------------------------------------------------------------------------
Ndays <- 360
Ngroups <- 9
default_e <- 0.95 # Baseline efficacy
kappa_default <- 0 # Originally marked as "losing immunity", which conflicts with phi, but it is always 0
default_supply_ceiling <- 1
default_group_seq <- TRUE # Vaccinate from older to younger

# Mortality risk p
ifr_hic <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3) / 100
ifr_lic <- ifr_hic * (3.2 / 2) ^ (5 : (-3))
default_pdeath <- ifr_hic

# Recall in the PNAS appendix, q is said to be set to match R0, the reproduction 
# number of the virus. The chunk below defines the mapping from R0 to q. 
ev <- eigen(default_cm)$values[1]
ev_pi <- eigen(default_cm*(1-pre_immunity))$values[1]
r0 <- function(r) r / (5 * ev_pi)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# 3. Two-dose model parameters
# ------------------------------------------------------------------------------
# This first one is the foundation upon which later ones simply modify
pars_fdf_slow <- lst(
  Nc = 13, 
  Ngroups, 
  Ndays,
  y0 = y0_gen(13, 9, pre_immunity, infected0[["slow"]]),
  q = rep(r0(1.1), Ngroups), # R0 == 1.1, as noted in PNAS appendix
  contacts = default_cm,
  gamma1 = rep(.2, Ngroups),
  gamma2 = rep(.2, Ngroups),
  delta1 = rep(0, Ngroups),
  delta2 = rep(0, Ngroups),
  kappa1 = rep(kappa_default, Ngroups),
  kappa2 = rep(kappa_default, Ngroups),
  phi = rep(0, Ngroups), 
  ta = rep(0, Ngroups),
  e1 = rep(.8, Ngroups),
  e2 = rep(.95, Ngroups),
  pdeath = default_pdeath,
  vrf = 1,
  vstop = rep(.8, Ngroups), # slow down when 80% vaccinated, may be modified in apap_*()
  constantrisk = 0
)

pars_fdf_fast <- list_modify(pars_fdf_slow,
                             y0 = y0_gen(13, 9, pre_immunity, infected0[["fast"]]),
                             q = rep(r0(2), Ngroups))
pars_fdf_linear <- list_modify(pars_fdf_slow,
                               y0 = y0_gen(13, 9, pre_immunity, infected0[["decreasing"]]),
                               q = rep(r0(.99), Ngroups))
pars_fdf_cr <- list_modify(pars_fdf_slow,
                           y0 = y0_gen(13, 9, pre_immunity, .1/30.5),
                           q = rep(0, Ngroups),
                           constantrisk = .01/30.5)
pars_fdf_end <- list_modify(pars_fdf_fast,
                            y0 = y0_gen(13, 9, rep(.5, Ngroups), infected0[["fast"]]))
# c(rep(12,1), 0) resets cumI (compartment 13) to 0:
set0 <- c(rep(1,4), 0, rep(1,7), 0)
pars_fdf_late <- list_modify(pars_fdf_fast, 
                             y0 = (sr(pars_fdf_fast, f = "2d_v2")["120", ,])*set0)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------





# ------------------------------------------------------------------------------
# 4. Two-vaccine model parameters
# ------------------------------------------------------------------------------
pars_le_slow <- list_modify(pars_fdf_slow,
                            e1 = 0.95, e2 = 0, 
                            ta1 = rep(0, Ngroups), 
                            ta2 = rep(0, Ngroups), 
                            tmore1 = rep(Inf, Ngroups),
                            tmore2 = rep(Inf, Ngroups),
                            ts1 = rep(Ndays, Ngroups))
pars_le_slow$ta  <- NULL
pars_le_fast <- list_modify(pars_le_slow,
                            y0 = y0_gen(13, 9, pre_immunity, infected0[["fast"]]),
                            q = rep(r0(2), Ngroups))
pars_le_cr <- list_modify(pars_le_slow,
                          y0 = y0_gen(13, 9, pre_immunity, .1/30.5),
                          q = rep(0, Ngroups),
                          constantrisk = .01/30.5)
pars_le_linear <- list_modify(pars_le_slow,
                              y0 = y0_gen(13, 9, pre_immunity, infected0[["decreasing"]]),
                              q = rep(r0(.99), Ngroups))

pars_le_late <- list_modify(pars_le_fast, 
                            y0 = (sr(pars_le_fast)["120", ,])*set0)




# ------------------------------------------------------------------------------
# 5. Selection of parameters
# ------------------------------------------------------------------------------


scenario_par_nms_2v <- c("pars_linear", "pars_le_slow", "pars_le_fast")#, "pars_le_late", "pars_le_cr")
scenario_nms_2v <- c("Slow-decrease epidemic", "Slow-growth epidemic", "Fast-growth epidemic")#, "Declining risk", "Constant risk")
scenario_list_2v <- lst(
  # "Constant risk of infection" = pars_le_cr,
  "Slow decrease (R = 0.99)" = pars_le_linear,
  "Slow growth (R = 1.1)" = pars_le_slow,
  "Fast growth (R = 2)" = pars_le_fast,
  # "Declining risk (R0 = 3, after peak)" = pars_le_late
) %>%
  setNames(scenario_nms_2v)


# Grab parameter sets (fast, slow, late, constant risk)
grab_2v_parms <- function(model){
  if(model == "pars_linear")  pars <- pars_le_linear
  if(model == "pars_le_cr")   pars <- pars_le_cr
  if(model == "pars_le_slow") pars <- pars_le_slow
  if(model == "pars_le_fast") pars <- pars_le_fast
  if(model == "pars_le_late") pars <- pars_le_late
  pars
}



# ------------------------------------------------------------------------------
# Unknown (seems important but don't know the use)
#-------------------------------------------------------------------------------

# Main FDF assumptions (Don't really know the use)
delay_default <- 28 - 10
delay_fdf <- 84 - 10
delay_hybrid_k <- sapply(1:8, function(k){
  c(rep(delay_fdf, k), rep(delay_default, 9-k))
})
colnames(delay_hybrid_k) <- 1:8
all_k <- TRUE

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



as.percent <- function(x, d=2, perc=FALSE){
  if (perc) paste0(format(round(100*x, d), nsmall = d),'%')
  else format(round(100*x, d), nsmall = 2)
}

