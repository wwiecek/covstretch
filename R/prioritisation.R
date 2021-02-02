# Function for deriving age prioritised vaccination parameters -----

# apap = adjust parameters for age prioritisation
# Derive values of delta and ta such that 
# len = len of vaccination campaign (in days), i.e. time to 100%
apap_2d <- function(pars, len, d2 = 18) {
  prop <- sum(pop[7:9])
  if(length(d2) == 1)
    d2 <- rep(d2, Ngroups)
  if(length(len) == 1)
    len <- rep(len, 9)
  list_modify(pars, 
              delta2 = 1/d2,
              ta = 10 + len*prop*c(rep(1,6),rep(0,3)),
              # delta1 = c(rep(1/len/(1-prop), 6), rep(1/len/prop, 3)))
              delta1 = 1/len/c(rep(1-prop, 6), rep(prop, 3)))
}


apap_2v <- function(pars, len, switch=Inf, delay = 10, n = pop) {
  prop <- sum(n[7:9])
  d1 <- c(rep(1/len/(1-prop), 6), rep(1/len/prop, 3))
  t1 <- delay + c(rep(len*prop, 6), 0, 0, 0)
  list_modify(pars, 
              ta1 = t1,
              ta2 = sapply(t1, function(x) max(x, switch+delay)),
              ts1 = rep(switch+delay, Ngroups),
              delta1 = d1,
              delta2 = d1)
}

# apap <- function(pars, len, d2 = 18) {
#   prop <- sum(pop[7:9])
#   list_modify(pars, 
#               delta2 = rep(1/d2, Ngroups),
#               ta = rep(10, 9),
#               delta1 = rep(1/len, Ngroups))
# }