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
  
  if(f == "2d") {
    mod <- odin_ode_2dose(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV", "cumI")
  }
  if(f == "2d_v2") {
    mod <- odin_ode_2dose_v2(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  }
  if(f == "2d_v3") {
    mod <- odin_ode_2dose_v3(user = pars)
    cnames <- c("S", "E0", "E1", "E2", "I0", "I1", "I2", 
                "R", "RV", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  }
  if(f == "2v") {
    mod <- odin_ode_2vaccines(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumI")
  }
  if(f == "2v_v2") {
    mod <- odin_ode_2vaccines_v2(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  }
  
  y <- array(mod$run(times)[,-1], 
             dim = c(max(times), pars$Ngroups, pars$Nc), 
             dimnames = list(times, gnames, cnames)) %>%
    aperm(c(1, 3, 2))
}