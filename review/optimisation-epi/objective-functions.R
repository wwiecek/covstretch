# ------------------------------------------------------------------------------
# Define objective functions for the dynamic or static optimization problem
# ------------------------------------------------------------------------------

# Dynamic objective function
# ------------------------------------------------------------------------------
# - scenario: model scenario such as "pars_fast", "pars_slow", etc.
# - length_campaign: length of the vaccination campaign in days, which regulates the supply constraint
# - fd: dose fraction
# - phi_x: dose response function
# - rm: whether the raw output from single_run will be returned
# - ret: whether the output compartment values are scaled by population
# - objective: which compartment value to be minimized, either "D" (death) or "cumI" (total infection)
# - homogen: whether the mixing (contacts across groups) is homogenous or heterogenous
# - rep: how many times vaccination recurs

model_fd_dynamic <- function(scenario, 
                             length_campaign, 
                             fd, 
                             phi_x = "covid_default",
                             pdeath = "ifr_hic",
                             rm = FALSE,
                             ret = 0,
                             objective = "D",
                             homogen = FALSE,
                             rep = c(0.8)) {
  if (phi_x == "covid_default") {phi_x <- function(x) -25.31701*x^1.037524 + 1.037524*25.31701*x
  } else if (phi_x == "flu_default") {phi_x <- function(x) 0}
  
  e1 <- phi_x(fd)
  
  if (pdeath == "ifr_hic") {pdeath <- ifr_hic} else {pdeath <- ifr_lic}
  
  pars <- apap_2v(grab_2v_parms(scenario), fractional_dose = fd, len = length_campaign)
  pars <- list_modify(pars, e1 = e1, pdeath = pdeath)
  if(homogen){
    pars$contacts <- t(replicate(Ngroups, pop))
    pars$q <- ev*pars$q
  }
  if (length(rep) == 1) {
    y <- sr(pars, "2v_v2")
  } else if (length(rep) > 1) {
    y <- multi_year_run(par, "2v_v2", rep)
  }
  if(rm) return(y)
  if(ret == 0)
    return(main_metrics(y, pop))
  if(ret == 1)
    y <- rescale_rcs(y, pop, TRUE)
    return(y[360 * length(rep),objective,1])
}
# ------------------------------------------------------------------------------


# Static objective function
# ------------------------------------------------------------------------------
# - scenario: model scenario such as "pars_le_fast", "pars_le_slow", etc.
# - fd: dose fraction
# - phi_x: dose response function
# - rm: whether the raw output from single_run will be returned
# - ret: whether the output compartment values are scaled by population
# - objective: which compartment value to be minimized, either "D" (death) or "cumI" (total infection)
# - homogen: whether the mixing (contacts across groups) is homogenous or heterogenous
# 
# Note that the supply constraint is not built into the objective function, 
# rather, it is built into the optimization problem as the optimization constraint. 
# Refer to line xxx of nlopt_general.R for more detail. 

model_fd_static <- function(scenario, 
                            fd, 
                            phi_x = "covid_default",
                            pdeath = "ifr_hic",
                            rm = FALSE,
                            ret = 0,
                            objective = "D",
                            homogen = FALSE,
                            full = 0,
                            rep) {

  if (phi_x == "covid_default") {phi_x <- function(x) -25.31701*x^1.037524 + 1.037524*25.31701*x
  } else if (phi_x == "flu_default") {phi_x <- function(x) 0}
  
  e_vector <- phi_x(fd)
  
  if (pdeath == "ifr_hic") {pdeath <- ifr_hic
  } else if (pdeath == "ifr_lic") {pdeath <- ifr_lic}
  
  if (full){
    e_vector <- fd*phi_x(1)
  }
  pars <- list_modify(
    grab_2v_parms(scenario),
    y0 = y0_gen(13, Ngroups, pre_immunity = pre_immunity + (1-pre_immunity)*e_vector))
  if(homogen){
    pars$contacts <- t(replicate(Ngroups, pop))
    pars$q <- ev*pars$q
  }
  pars <- list_modify(pars, pdeath = pdeath)
  # We do not update e1, because there is no vaccination past t=0 
  if (length(rep) == 1) {
    y <- sr(pars, "2v_v2")
  } else if (length(rep) > 1) {
    y <- multi_year_run(par, "2v_v2", rep)
  }
  if(ret == 0)
    return(main_metrics(y, pop)[1:2])
  if(ret == 1)
    y <- rescale_rcs(y, pop, TRUE)
    return(y[360 * length(rep),objective,1])
  
}
