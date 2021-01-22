

# Running the model function -----

# Single run = sr()
sr <- function(pars, f = "2v") { 
  gnames <- 1:pars$Ngroups
  if(pars$Ngroups == 9)
    gnames <- colnames(pbc_spread)
  
  if(f == "2d") {
    mod <- odin_ode_2dose(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV", "cumI")
  }
  if(f == "2v") {
    mod <- odin_ode_2vaccines(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumI")
  }
  times <- 1:(pars$Ndays)
  y <- array(mod$run(times)[,-1], 
             dim = c(pars$Ndays, pars$Ngroups, pars$Nc), 
             dimnames = list(times, gnames, cnames)) %>%
    aperm(c(1, 3, 2))
}


# Generate a starting population -----
# with varying levels of initial
# infection, pre-existing immunity, appropriate dimensions
y0_gen <- function(Nc, Ngroups, 
                   pre_immunity = rep(0, Ngroups), 
                   # v = rep(0, Ngroups), 
                   ii = 5e-03){
  y0_default <- matrix(0, Nc, Ngroups)
  y0_default[1,] <- 1-pre_immunity
  y0_default[2,] <- (y0_default[1,]*ii)/2
  y0_default[3,] <- (y0_default[1,]*ii)/2
  # subtract initial infections
  y0_default[1,] <- y0_default[1,] - y0_default[2,]
  y0_default[4,] <- pre_immunity
  y0_default
}


# Gathering outputs (burden of disease) -----

# Infection burden (only if no loss of immunity):
bi <- function(sr, pop, sub = 0) { 
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  y[dim(y)[1], "R", 1]
}
# Death burden:
bd <- function(sr, pop) {
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  y[dim(y)[1], "D", 1]
}
b_any <- function(sr, pop, comp) {
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  y[dim(y)[1], comp, 1]
}
