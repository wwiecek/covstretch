# install.packages("nloptr")
library(nloptr)
library(tidyverse)
source("project-setup.R")

source("optimisation-epi/objective-functions.R")

prop_adults <- sum(pop[3:9])/sum(pop) #for now I ignore children, so Q is scaled down to adult pop
outcome_nlopt <- "D"
q_seq <- c(0.1, 0.25, 0.5, 0.75, 1) #just try two quantities at first


# Choose dimensionality of x here: -----
n_x <- 7
unroll_x <- function(x, sub=1) #set sub to 0 for static model where no doses go to children
  c(sub,sub,x[1],x[2],x[3],x[4],x[5],x[6],x[7])


# Solve for various Q's (static) or deltas (dynamic) ------

opt_problem <- function(q_seq,h=FALSE, static = TRUE){
  
  sol <- sapply(q_seq, function(Q) {
    # Optimisation problem: define objective function for a given Q, h:
    if(static)
      eval_f <- function(x) model_fd_static(unroll_x(x,sub=0), 
                                            homogen = h, 
                                            outcome = outcome_nlopt,
                                            ret = 1)
    else
      eval_f <- function(x) model_fd_dynamic(model = "pars_le_fast", 
                                             fd = unroll_x(x), 
                                             d1 = Q, 
                                             outcome = outcome_nlopt,
                                             ret = 1,
                                             homogen = h)
    
    # Equality constraints
    if(static)
      #GROUP NUMBER - Need to change this constraint according to the number of groups
      eval_g_ineq <- function(x) c(x[1]*pop[3]+x[2]*pop[4]+x[3]*pop[5]+x[4]*pop[6]+x[5]*pop[7]+x[6]*pop[8]+x[7]*pop[9]-Q,
                                   -x[1]*pop[3]-x[2]*pop[4]-x[3]*pop[5]-x[4]*pop[6]-x[5]*pop[7]-x[6]*pop[8]-x[7]*pop[9]+Q)
    else
      eval_g_ineq <- NULL
    
    # Lower and upper bounds
    ub <- rep(1,n_x)
    
    # Initial values
    if(static){
      lb <- rep(0,n_x)
      x0 <- rep(0.5, n_x)
    } else {
      lb <- rep(0.01,n_x)
      x0 <- rep(0.01, n_x)
    }
    # Set optimization options.
    opts <- list( "algorithm"= "NLOPT_LN_COBYLA",#"NLOPT_GN_ISRES"
                  "xtol_rel"= 1.0e-4,
                  "maxeval"= 100,
                  #"local_opts" = list( "algorithm" = "NLOPT_LN_NELDERMEAD",
                  #                     "xtol_rel" = 1.0e-4 ),
                  "print_level" = 0 )
    
    
    # Solving
    res <- nloptr ( x0 = x0,
                    eval_f = eval_f,
                    lb = lb,
                    ub = ub,
                    eval_g_ineq = eval_g_ineq,
                    opts = opts
    )
    c(res$objective, res$solution)
  })
  
  colnames(sol) <- c(q_seq)
  sol
}

nl_q_seq <- seq(0.1, 0.7, 0.1)#Variance analysis are unnecessary since the solution only changes with the initial values
nl_d_seq <- c(1000,400,200,100,50)
ptm <- proc.time()
nlopt_s0 <- opt_problem(nl_q_seq, 0)
proc.time() - ptm
nlopt_s1 <- opt_problem(nl_q_seq, 1)
ptm <- proc.time()
nlopt_d0 <- opt_problem(nl_d_seq, 0, static = F)
proc.time() - ptm
nlopt_d1 <- opt_problem(nl_d_seq, 1, static = F)

#save(nl_q_seq, nlopt_s0, file = "results/wip-nl-solutions5-s-h0-D.Rdata")
save(nl_q_seq, nl_d_seq,
     nlopt_s0, nlopt_s1, nlopt_d0, nlopt_d1, file = "results/wip-nl-solutions7-D.Rdata")
load(file = "results/wip-nl-solutions7-D.Rdata")
