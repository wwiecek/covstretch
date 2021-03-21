# install.packages("nloptr")
library(nloptr)
library(tidyverse)
source("project-setup.R")


# Objective Function
# eval_f <- function(x)
#   max((sr(apap_2v(pars_le_fast, 100, fractional_dose = c(1,1, x)), f = "2v_v2") %>% 
#          rescale_rcs(pop, TRUE))[360,outcome,1])

phi_x <- function(x) 
  sapply(3.49706*sqrt(x) - 1.74853*x  -0.798528, function(x) max(c(0,x)))
# 0.95*x^0.25



prop_adults <- sum(pop[3:9])/sum(pop) #for now I ignore children, so Q is scaled down to adult pop
outcome <- "cumI"
# q_seq <- seq(0.05, 1, 0.1)
q_seq <- c(0.1, 0.25, 0.5, 0.75, 1) #just try two quantities at first


# Choose dimensionality of x here: -----
n_x <- 2
x_unroll <- function(x)
  c(0,0,x[1],x[1],x[1],x[1],x[2],x[2],x[2])


# Objective function (working with x of any length) -----
model_opt_static <- function(v_prop, homogen = FALSE) {
  e_vector <- sapply(x_unroll(v_prop), phi_x)
  pars <- list_modify(
    pars_le_fast,
    y0 = y0_gen(13, Ngroups, pre_immunity = pre_immunity + (1-pre_immunity)*e_vector))
  
  if(homogen){
    pars$contacts <- 1/Ngroups + 0*pars$contacts
    pars$q <- ev*pars$q
  }
  
  y <- rescale_rcs(sr(pars, "2v_v2"), pop, TRUE)
  y[360,outcome,1]
}


# Solve for various Q ------



opt_problem <- function(q_seq,h=FALSE){
  
  eval_f <- function(x) model_opt_static(x,h)
  
  q_sol <- sapply(prop_adults*q_seq, function(Q) {
    
    # Equality constraints
    eval_g_eq <- function(x) c(sum(x_unroll(x)*pop)-Q)
    
    # Lower and upper bounds
    lb <- rep(0,n_x)
    ub <- rep(1,n_x)
    
    # Initial values
    # x0 <- c(.215,.052)
    x0 <- rep(0, n_x)
    
    # Set optimization options.
    opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                  "xtol_rel"= 1.0e-4,
                  "maxeval"= 5e02,
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
  
  colnames(q_sol) <- c(q_seq)
  q_sol
}

q_seq <- rep(seq(0.1, 0.7, 0.1), each = 5)
h0 <- opt_problem(q_seq, 0)
h1 <- opt_problem(q_seq, 1)

save(h0, h1, file = "results/wip-nl-solutions.Rdata")

# as.data.frame(h0[-1,]) %>% 
#   mutate(age = 1:n_x) %>% 
#   gather(var, val, -age) %>% 
#   mutate(age = factor(age)) %>% 
#   ggplot(aes(color=age,y=val,x=var, group=age)) + 
#   geom_line()





