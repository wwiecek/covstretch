# install.packages("nloptr")
library(nloptr)
library(tidyverse)
source("project-setup.R")

source("optimisation-epi/objective-functions.R")

prop_adults <- sum(pop[3:9])/sum(pop) #for now I ignore children, so Q is scaled down to adult pop
outcome_nlopt <- "D"
q_seq <- c(0.1, 0.25, 0.5, 0.75, 1) #just try two quantities at first


# Choose dimensionality of x here: -----
n_x <- 4
unroll_x <- function(x, sub=1) #set sub to 0 for static model where no doses go to children
  c(sub,sub,x[1],x[1],x[2],x[2],x[3],x[3],x[4])


# Solve for various Q's (static) or deltas (dynamic) ------

opt_problem <- function(q_seq,h=FALSE, static = TRUE){
  
  sol <- sapply(q_seq, function(Q) {
    
    # Optimisation problem: define objective function for a given Q, h:
    if(static)
      eval_f <- function(x) model_fd_static(unroll_x(x), 
                                            homogen = h, 
                                            outcome = outcome_nlopt)
    else
      eval_f <- function(x) model_fd_dynamic(model = "pars_le_fast", 
                                             fd = unroll_x(x), 
                                             d1 = Q, 
                                             outcome = outcome_nlopt,
                                             homogen = h)
    
    # Equality constraints
    if(static)
      eval_g_eq <- function(x) c(sum(unroll_x(x)*pop)-Q)
    else
      eval_g_eq <- NULL
    
    # Lower and upper bounds
    lb <- rep(0,n_x)
    ub <- rep(1,n_x)
    
    # Initial values
    if(static)
      x0 <- rep(0, n_x)
    else
      x0 <- rep(1, n_x)
    
    # Set optimization options.
    opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                  "xtol_rel"= 1.0e-4,
                  "maxeval"= 1e03,
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

# save(nl_q_seq, nl_d_seq,
     # nlopt_s0, nlopt_s1, nlopt_d0, nlopt_d1, file = "results/wip-nl-solutions4.Rdata")
