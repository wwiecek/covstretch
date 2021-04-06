# install.packages("nloptr")
library(nloptr)
library(tidyverse)
library(ggplot2)
source("project-setup.R")

source("optimisation-epi/objective-functions.R")

prop_adults <- sum(pop[3:9])/sum(pop) #for now I ignore children, so Q is scaled down to adult pop
q_seq <- c(0.1, 0.25, 0.5, 0.75, 1) #just try two quantities at first


# Choose dimensionality of x here: -----
n_x <- 7
unroll_x <- function(x, sub=1) #set sub to 0 for static model where no doses go to children
  c(sub,sub,x[1],x[2],x[3],x[4],x[5],x[6],x[7])

# Solve for various Q's (static) or deltas (dynamic) ------

opt_problem <- function(q_seq,h=FALSE, static = TRUE, x0=0){
  
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
      x0 <- rep(x0, n_x)
    } else {
      lb <- rep(0.01,n_x)
      x0 <- rep(x0, n_x)
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

for (initial_value in c(0.3,0.5,0.8,0.99,0.01)){
  for (outcome_nlopt in c("D","cumI")){
    nl_q_seq <- seq(0.1, 0.7, 0.1)#Variance analysis are unnecessary since the solution only changes with the initial values
    nl_d_seq <- c(1000,400,200,100,50)
    nlopt_s0 <- opt_problem(nl_q_seq, 0, x0=initial_value)
    nlopt_s1 <- opt_problem(nl_q_seq, 1, x0=initial_value)
    nlopt_d0 <- opt_problem(nl_d_seq, 0, static = F, x0=initial_value)
    nlopt_d1 <- opt_problem(nl_d_seq, 1, static = F, x0=initial_value)
    
    #save(nl_q_seq, nlopt_s0, file = "results/wip-nl-solutions5-s-h0-D.Rdata")
    save(nl_q_seq, nl_d_seq,
         nlopt_s0, nlopt_s1, nlopt_d0, nlopt_d1, file = paste0("results/nlopt/wip-nl-solutions7-",outcome_nlopt,"-",initial_value*100,".Rdata"))
  }
}

##Comparing results - auxiliary functions----
eval_constr <- function(i,x) x[i,1]*pop[3]+x[i,2]*pop[4]+x[i,3]*pop[5]+x[i,4]*pop[6]+x[i,5]*pop[7]+x[i,6]*pop[8]+x[i,7]*pop[9]-as.numeric(rownames(x)[i])

fmt_res <- function(res) {
  x0<-data.frame(t(res[2:8,]))
  x0['Q'] <- rownames(x0)
  x0['constr'] <- sapply(1:dim(x0)[1], eval_constr, x=x0)
  x0
}

fmt_df_nlopt <- function(x,dim_value,static=FALSE,dim='initial_value'){
  res0 <- x[1,]
  x0 <- data.frame(t(x[2:8,]))
  colnames(x0) <- c('age3','age4','age5','age6','age7','age8','age9')
  if (static){
    x0['Q'] <- rownames(x0)
    x0['constr'] <- sapply(1:dim(x0)[1], eval_constr, x=x0)
  } else {
    x0['d1'] <- rownames(x0)
  }
  x0['d1'] <- rownames(x0)
  x0['res'] <- res0
  x0[dim] <- dim_value
  x0['outcome'] <- outcome_nlopt
  x0
}

##Comparing results - initial values----
nlopt_s0_all <- data.frame()
nlopt_s1_all <- data.frame()
nlopt_d0_all <- data.frame()
nlopt_d1_all <- data.frame()
for (initial_value in c(0.3,0.5,0.8,0.99,0.01)){
  for (outcome_nlopt in c("D","cumI")){
    load(paste0("results/nlopt/wip-nl-solutions7-",outcome_nlopt,"-",initial_value*100,".Rdata"))
    nlopt_d0_all <- rbind(nlopt_d0_all,fmt_df_nlopt(nlopt_d0,initial_value))
    nlopt_d1_all <- rbind(nlopt_d1_all,fmt_df_nlopt(nlopt_d1,initial_value))
    nlopt_s0_all <- rbind(nlopt_s0_all,fmt_df_nlopt(nlopt_s0,initial_value,static=T))
    nlopt_s1_all <- rbind(nlopt_s1_all,fmt_df_nlopt(nlopt_s1,initial_value,static=T))
  }
}

#Dynamic
nlopt_d0_all.fig <- nlopt_d0_all %>%
  mutate(d1 = factor(d1, levels = c("50","100","200","400","1000"), labels = c("50","100","200","400","1000"))) %>%
  mutate(initial_value = factor(initial_value)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = d1, y = res, color = initial_value)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/d0_initial_values.pdf",nlopt_d0_all.fig)

nlopt_d1_all.fig <- nlopt_d1_all %>%
  mutate(d1 = factor(d1, levels = c("50","100","200","400","1000"), labels = c("50","100","200","400","1000"))) %>%
  mutate(initial_value = factor(initial_value)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = d1, y = res, color = initial_value)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/d1_initial_values.pdf",nlopt_d1_all.fig)

#Static - objective function
nlopt_s0_all.fig <- nlopt_s0_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(initial_value = factor(initial_value)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = res, color = initial_value)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s0_initial_values_obj.pdf",nlopt_s0_all.fig)

nlopt_s1_all.fig <- nlopt_s1_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(initial_value = factor(initial_value)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = res, color = initial_value)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s1_initial_values_obj.pdf",nlopt_s1_all.fig)

#Static - constraint
nlopt_s0_all.fig <- nlopt_s0_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(initial_value = factor(initial_value)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = constr, color = initial_value)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s0_initial_values_constr.pdf",nlopt_s0_all.fig)

nlopt_s1_all.fig <- nlopt_s1_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(initial_value = factor(initial_value)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = constr, color = initial_value)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s1_initial_values_constr.pdf",nlopt_s1_all.fig)

##Comparing results - stopping criteria----
nlopt_s0_all <- data.frame()
nlopt_s1_all <- data.frame()
nlopt_d0_all <- data.frame()
nlopt_d1_all <- data.frame()
for (tol in c('maxeval1000-tol-6','maxeval500-tol-5','maxeval100-tol-4')){
  for (outcome_nlopt in c("D","I")){
    load(paste0("results/nlopt/wip-nl-solutions7-",outcome_nlopt,"-",tol,".Rdata"))
    nlopt_d0_all <- rbind(nlopt_d0_all,fmt_df_nlopt(nlopt_d0,tol,dim='tol'))
    nlopt_d1_all <- rbind(nlopt_d1_all,fmt_df_nlopt(nlopt_d1,tol,dim='tol'))
    nlopt_s0_all <- rbind(nlopt_s0_all,fmt_df_nlopt(nlopt_s0,tol,static=T,dim='tol'))
    nlopt_s1_all <- rbind(nlopt_s1_all,fmt_df_nlopt(nlopt_s1,tol,static=T,dim='tol'))
  }
}

#Dynamic
nlopt_d0_all.fig <- nlopt_d0_all %>%
  mutate(d1 = factor(d1, levels = c("50","100","200","400","1000"), labels = c("50","100","200","400","1000"))) %>%
  mutate(initial_value = factor(tol)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = d1, y = res, color = tol)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/d0_convergence.pdf",nlopt_d0_all.fig)

nlopt_d1_all.fig <- nlopt_d1_all %>%
  mutate(d1 = factor(d1, levels = c("50","100","200","400","1000"), labels = c("50","100","200","400","1000"))) %>%
  mutate(initial_value = factor(tol)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = d1, y = res, color = tol)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/d1_convergence.pdf",nlopt_d1_all.fig)

#Static - objective function
nlopt_s0_all.fig <- nlopt_s0_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(tol = factor(tol)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = res, color = tol)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s0_convergence_obj.pdf",nlopt_s0_all.fig)

nlopt_s1_all.fig <- nlopt_s1_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(tol = factor(tol)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = res, color = tol)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s1_convergence_obj.pdf",nlopt_s1_all.fig)

#Static - constraint
nlopt_s0_all.fig <- nlopt_s0_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(tol = factor(tol)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = constr, color = tol)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s0_convergence_constr.pdf",nlopt_s0_all.fig)

nlopt_s1_all.fig <- nlopt_s1_all %>%
  mutate(Q = factor(Q)) %>%
  mutate(tol = factor(tol)) %>%
  mutate(outcome = factor(outcome)) %>%
  ggplot(aes(x = Q, y = constr, color = tol)) + 
  geom_point(size = 2.5) + facet_grid(outcome ~ .,scales="free")
ggsave("figures/nlopt/s1_convergence_constr.pdf",nlopt_s1_all.fig)
