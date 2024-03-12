# ------------------------------------------------------------------------------
# Define the optimization problem of fractional dose
# ------------------------------------------------------------------------------

# Helper objects
prop_adults <- sum(pop[3:9])/sum(pop) #for now I ignore children, so Q is scaled down to adult pop
q_seq <- c(0.1, 0.25, 0.5, 0.75, 1) #just try two quantities at first


# Choose dimensionality of x here: -----
n_x <- 7
# these are our controls - how many doses go to each age bin
unroll_x <- function(x, sub=1) #set sub to 0 for static model where no doses go to children
  c(sub,sub,x[1],x[2],x[3],x[4],x[5],x[6],x[7])


# Function Arguments:
# 
# - q_seq: supply constraints (optimization is run for each constraint)
#          * (dynamic) length of vaccination campaign in days
#          * (static) total dose amount available (percentage of full dose for entire population)
# - initial_value: initial value of dose fraction, same across all age groups
# - objective: either "D" (death) or "cumI" (cumulative infection) 
# - dose_response: function phi(x) of efficacy response to dose fraction
# - supply_constraint:
# - homogen_mixing: TRUE or FALSE
# - static: TRUE or FALSE for static problem
# - pdeath: mortality risk vector, default to high-income country case
# - scenario: cases decreasing, slowly increasing, rapidly increasing, etc
# - recurring: TRUE or FALSE for recurring periods with re-vaccination
opt_general <- 
  function(q_seq,
           initial_value,
           objective = "D",
           dose_response = function(x) -25.31701*x^1.037524 + 1.037524*25.31701*x,
           homogen_mixing = F, 
           static = T,
           pdeath = ifr_hic,
           scenario = "pars_le_slow",
           recurring = T,
           iterations = 100) {
    default_pdeath <- pdeath
  
    # Run optimization for each supply constraint
    sol <- sapply(q_seq, function(Q) {
      # Select objective function
      if(static)
        eval_f <- function(x) model_fd_static(scenario = scenario,
                                              fd = unroll_x(x,sub=0), 
                                              phi_x = dose_response, 
                                              objective = objective,
                                              homogen = homogen_mixing,
                                              ret = 1)
      else
        eval_f <- function(x) model_fd_dynamic(scenario = scenario,
                                               fd = unroll_x(x), 
                                               length_campaign = Q,
                                               phi_x = dose_response,
                                               objective = objective,
                                               ret = 1,
                                               homogen = homogen_mixing)
      # Set equality constraints
      if(static)
        # In the static case, the total dose applied is equal to the total dose available
        # Need to change this constraint according to the number of groups
        eval_g_ineq <- function(x) c(x[1]*pop[3]+x[2]*pop[4]+x[3]*pop[5]+x[4]*pop[6]+x[5]*pop[7]+x[6]*pop[8]+x[7]*pop[9]-Q,
                                   -x[1]*pop[3]-x[2]*pop[4]-x[3]*pop[5]-x[4]*pop[6]-x[5]*pop[7]-x[6]*pop[8]-x[7]*pop[9]+Q)
      else
        eval_g_ineq <- NULL
    
    # Lower and upper bounds
    ub <- rep(1,n_x)
    
    # Initial values
    if(static){
      lb <- rep(0,n_x)
      x0 <- rep(initial_value, n_x)
    } else {
      lb <- rep(0.01,n_x)
      x0 <- rep(initial_value, n_x)
    }
    # Set optimization options.
    opts <- list( "algorithm"= "NLOPT_LN_COBYLA",#"NLOPT_GN_ISRES"
                  "xtol_rel"= 0,
                  "maxeval"= iterations,
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
