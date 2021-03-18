# install.packages("nloptr")
library(nloptr)
library(tidyverse)
source("run-all.R")


# Objective Function
# eval_f <- function(x)
#   max((sr(apap_2v(pars_le_fast, 100, fractional_dose = c(1,1, x)), f = "2v_v2") %>% 
#          rescale_rcs(pop, TRUE))[360,outcome,1])

phi_x <- function(x) 
  3.49706*sqrt(x) - 1.74853*x  -0.798528



prop_adults <- sum(pop[3:9])/sum(pop) #for now I ignore children, so Q is scaled down to adult pop
outcome <- "D"
# q_seq <- seq(0.1, 1, 0.1)
q_seq <- c(0.5, 1) #just try two quantities at first

eval_f <- function(v_prop) {
  e_vector <- c(0, 0, sapply(v_prop, phi_x))
  # Pre-vaccinate individuals:
  pars <- list_modify(
    pars_le_fast,
    y0 = y0_gen(13, Ngroups, pre_immunity = pre_immunity + (1-pre_immunity)*e_vector*c(1,1,v_prop)))
  
  y <- rescale_rcs(sr(pars, "2v_v2"), pop, TRUE)
  y[360,outcome,1]
}

q_sol <- sapply(prop_adults*q_seq, function(Q) {

  # Equality constraints
  eval_g_eq <- function(x)
  {
    return ( c(sum(x*pop)-Q) )
  }
  
  # Lower and upper bounds
  lb <- rep(0,7)
  ub <- rep(1,7)
  
  # Initial values
  # x0 <- rep(1,7)
  x0 <- rep(0.2, 7)
  
  # Set optimization options.
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-4,
                "maxeval"= 1000,
                "local_opts" = list( "algorithm" = "NLOPT_LN_NELDERMEAD", 
                                     "xtol_rel" = 1.0e-4 ),
                "print_level" = 0 )
  
  # Solving
  res <- nloptr ( x0 = x0,
                  eval_f = eval_f,
                  lb = lb,
                  ub = ub,
                  eval_g_eq = eval_g_eq,
                  # eval_g_ineq = eval_g_eq,
                  opts = opts
  )
  c(res$objective, res$solution)
  # res$objective
})

colnames(q_sol) <- c(outcome, c(0,0,q_seq))

q_sol

# Constraints met?
apply(pop[3:9]*q_sol[-1,], 2, sum)