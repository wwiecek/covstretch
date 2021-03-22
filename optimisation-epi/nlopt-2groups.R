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
x_unroll <- function(x, sub=1) #set sub to 0 for static model where no doses go to children
  c(sub,sub,x[1],x[1],x[1],x[1],x[2],x[2],x[2])


# Objective function (working with x of any length) -----
model_nlopt_static <- function(v_prop, homogen = FALSE) {
  e_vector <- sapply(x_unroll(v_prop,0), phi_x)
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

model_nlopt_dynamic <- function(fd, d1, homogen = FALSE) {
  fd <- x_unroll(fd, 1)
  pars <- apap_2v(pars_le_fast, fractional_dose = fd, len = d1)
  pars$e1 <- sapply(fd, phi_x)
  if(homogen){
    pars$contacts <- 1/Ngroups + 0*pars$contacts
    pars$q <- ev*pars$q
  }
  
  y <- rescale_rcs(sr(pars, "2v_v2"), pop, TRUE)
  y[360,outcome,1]
}


# Solve for various Q's (static) or deltas (dynamic) ------



opt_problem <- function(q_seq,h=FALSE, static = TRUE){
  
  sol <- sapply(q_seq, function(Q) {
    
    # Optimisation problem
    if(static)
      eval_f <- function(x) model_nlopt_static(x,h)
    else
      eval_f <- function(x) model_nlopt_dynamic(x,d1 = Q, h)
    
    # Equality constraints
    if(static)
      eval_g_eq <- function(x) c(sum(x_unroll(x)*pop)-Q)
    else
      eval_g_eq <- NULL
    
    # Lower and upper bounds
    lb <- rep(0,n_x)
    ub <- rep(1,n_x)
    
    # Initial values
    # x0 <- c(.215,.052)
    if(static)
      x0 <- rep(0, n_x)
    else
      x0 <- rep(1, n_x)
    
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
                    opts = opts
    )
    c(res$objective, res$solution)
  })
  
  colnames(sol) <- c(q_seq)
  sol
}

nl_q_seq <- rep(seq(0.1, 0.7, 0.1), each = 5)
nl_d_seq <- rep(c(1000,400,200,100,50), each=5)
nlopt_s0 <- opt_problem(nl_q_seq, 0)
nlopt_s1 <- opt_problem(nl_q_seq, 1)
nlopt_d0 <- opt_problem(nl_d_seq, 0, static = F)
nlopt_d1 <- opt_problem(nl_d_seq, 1, static = F)

save(nl_q_seq, nl_d_seq,
     nlopt_s0, nlopt_s1, nlopt_d0, nlopt_d1, file = "results/wip-nl-solutions.Rdata")

# as.data.frame(h0[-1,]) %>% 
#   mutate(age = 1:n_x) %>% 
#   gather(var, val, -age) %>% 
#   mutate(age = factor(age)) %>% 
#   ggplot(aes(color=age,y=val,x=var, group=age)) + 
#   geom_line()
