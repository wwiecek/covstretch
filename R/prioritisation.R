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


apap_2v <- function(pars, len1, len2=Inf, n = pop) {
  prop <- sum(n[7:9])
  list_modify(pars, 
              ta1 = 10 + c(rep(len1*prop, 6), 0, 0, 0),
              delta1 = c(rep(1/len1/(1-prop), 6), rep(1/len1/prop, 3)))
}

# apap <- function(pars, len, d2 = 18) {
#   prop <- sum(pop[7:9])
#   list_modify(pars, 
#               delta2 = rep(1/d2, Ngroups),
#               ta = rep(10, 9),
#               delta1 = rep(1/len, Ngroups))
# }