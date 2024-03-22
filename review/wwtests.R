# Demonstration of a model with loss of immunity
par <- apap_2v(grab_2v_parms("pars_le_fast"), len = 1000)
par <- list_modify(par, 
                   Ndays = 720,
                   kappa1 = rep(1/365, 9),
                   phi = rep(1/365, 9))
y <- sr(par)
y
plot_rcs(y, c("S", "I", "cumV", "R"))

# Chain a few models
multi_year_run <- function(par, nperiods = 2) {
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
    cpar <- list_modify(par, y0 = last_state, vstop = rep(0.5, 9))
    cy <- sr(cpar)
    # append new simulation at the end of previous one
    y <- abind(y, cy, along = 1)
  }
  dimnames(y)[[1]] <- as.character(1:dim(y)[1])
  y
}

library(abind)
par <- apap_2v(grab_2v_parms("pars_le_fast"), fractional_dose = rep(1, 9), len = 100)
par <- list_modify(par, 
                   phi = rep(1/365, 9),
                   kappa1 = rep(1/365, 9))
multi_year_run(par, nperiods = 5) %>% plot_rcs(c("I", "P1", "S", "cumV"))
