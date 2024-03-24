library(abind)

sr <- function(pars, f = "2v_v2") { 
  gnames <- 1:pars$Ngroups
  times <- 1:(pars$Ndays)
  pars$Ndays <- NULL
  
  pars_to_expand <- c("e1", "e2")
  for(p in pars_to_expand){
    if(length(pars[[p]]) == 1)
      pars[[p]] <- rep(pars[[p]], pars$Ngroups)
  }
  
  
  if(pars$Ngroups == 9)
    gnames <- colnames(pbc_spread)
  
  # if(f == "2d") {
  #   mod <- odin_ode_2dose(user = pars)
  #   cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV", "cumI")
  # }
  # if(f == "2d_v2") {
  #   mod <- odin_ode_2dose_v2(user = pars)
  #   cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  # }
  # if(f == "2d_v3") {
  #   mod <- odin_ode_2dose_v3(user = pars)
  #   cnames <- c("S", "E0", "E1", "E2", "I0", "I1", "I2", 
  #               "R", "RV", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  # }
  # if(f == "2v") {
  #   mod <- odin_ode_2vaccines(user = pars)
  #   cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumI")
  # }
  if(f == "2v_v2") {
    mod <- odin_ode_2vaccines_v2$new(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  }
  y <- array(mod$run(times)[,-1], 
             dim = c(max(times), pars$Ngroups, pars$Nc), 
             dimnames = list(times, gnames, cnames)) %>%
    aperm(c(1, 3, 2))
}

# Chain a few models
multi_year_run <- function(par, f = "2v_v2", takeup = list(rep(0.8, 9))) {
  nperiods = length(takeup)
  y <- sr(par) #"joined" result (see below)
  cy <- y #current result
  if(nperiods == 1)
    return(y)
  
  for(i in 1:(nperiods-1)){
    # each new period just takes the last day of last period as starting condition
    last_state <- cy[par$Ndays,,]
    # reset vaccination indicators to zero ahead of the new season
    last_state["S",] <- last_state["S",] + last_state["N1",]
    last_state["R",] <- last_state["R",] + last_state["P1",]
    last_state["cumV1",] <-
      last_state["cumV2",] <-
      last_state["cumV",] <-
      last_state["P1",] <- 
      last_state["N1",] <- 
      0.0005 + 0*last_state["cumV1",]
    # simulate another period
    cpar <- list_modify(par, y0 = last_state, vstop = takeup[[i]])
    cy <- sr(cpar, f)
    # append new simulation at the end of previous one
    y <- abind(y, cy, along = 1)
  }
  dimnames(y)[[1]] <- as.character(1:dim(y)[1])
  y
}
