# install.packages("nloptr")
library(nloptr)
library(tidyverse)
library(xtable)
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

opt_problem <- function(q_seq,h=FALSE, static = TRUE, model = "pars_le_fast"){
  
  sol <- sapply(q_seq, function(Q) {
    # Optimisation problem: define objective function for a given Q, h:
    if(static)
      eval_f <- function(x) model_fd_static(unroll_x(x,sub=0), 
                                            homogen = h, 
                                            outcome = outcome_nlopt,
                                            ret = 1)
    else
      eval_f <- function(x) model_fd_dynamic(model = model, 
                                             fd = unroll_x(x), 
                                             d1 = Q, 
                                             outcome = outcome_nlopt,
                                             ret = 1,
                                             homogen = h)
    
    # Equality constraints
    if(static){
      # if (h==0)
      #GROUP NUMBER - Need to change this constraint according to the number of groups and if vaccinating young population
      eval_g_ineq <- function(x) c(0*pop[1]+0*pop[2]+x[1]*pop[3]+x[2]*pop[4]+x[3]*pop[5]+x[4]*pop[6]+x[5]*pop[7]+x[6]*pop[8]+x[7]*pop[9]-Q,
                                  -0*pop[1]-0*pop[2]-x[1]*pop[3]-x[2]*pop[4]-x[3]*pop[5]-x[4]*pop[6]-x[5]*pop[7]-x[6]*pop[8]-x[7]*pop[9]+Q)
      # else
      #   eval_g_ineq <- function(x) c(0*pop[1]+0*pop[2]+x[1]*pop[3]+x[2]*pop[4]+x[3]*pop[5]+x[4]*pop[6]+x[5]*pop[7]+x[6]*pop[8]+x[7]*pop[9]-Q,
      #                                -0*pop[1]-0*pop[2]-x[1]*pop[3]-x[2]*pop[4]-x[3]*pop[5]-x[4]*pop[6]-x[5]*pop[7]-x[6]*pop[8]-x[7]*pop[9]+Q,
      #                                x[1:6]-x[2:7]+1e-4)
    } else {
      eval_g_ineq <- NULL
    }
    
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

nl_q_seq <- seq(0.1, 1, 0.1)#Variance analysis are unnecessary since the solution only changes with the initial values
# nl_q_seq <- seq(0.1, 0.7, 0.1)
nl_d_seq <- c(1000,400,200,100,50)
#Q is adjusted in the next lines to match the percentage of adults (comparison with non-epi results), to change this remove prop_adults
nlopt_s0 <- opt_problem(nl_q_seq*prop_adults, 0)
nlopt_s1 <- opt_problem(nl_q_seq*prop_adults, 1)
nlopt_d0 <- opt_problem(nl_d_seq, 0, static = F)
nlopt_d0_slow <- opt_problem(nl_d_seq, 0, static = F, model = "pars_le_slow")
nlopt_d0_decline <- opt_problem(nl_d_seq, 0, static = F, model = "pars_le_linear")
nlopt_d1 <- opt_problem(nl_d_seq, 1, static = F)

#save(nl_q_seq, nlopt_s0, file = "results/wip-nl-solutions5-s-h0-D.Rdata")
# save(nl_d_seq, nlopt_d0, nlopt_d1, file = "results/nlopt/wip-nl-solutions7-d-D-hic-newconstr.Rdata")
save(nl_q_seq, nl_d_seq,
    nlopt_s0, nlopt_s1, nlopt_d0, nlopt_d0_slow, nlopt_d0_decline,
    nlopt_d1, file = "results/nlopt/wip-nl-solutions7-D-lic-edit_contacts.Rdata")


#Generating tables for report
#LIC
prop_adults_lic <- sum(lic_pop[3:9])/sum(lic_pop)
load("results/nlopt/wip-nl-solutions7-D-lic-edit_contacts.Rdata")
table_s0_lic<-data.frame(nlopt_s0[2:8,as.character(rev(nl_q_seq*prop_adults_lic))])
colnames(table_s0_lic) <- as.character(rev(nl_q_seq))
table_s0_lic["Age Group"] = c("20-30","30-40","40-50","50-60","60-70","70-80","80+")
table_s0_lic["Population Share"] = lic_pop[3:9]/sum(lic_pop[3:9])
table_s0_lic["Infection Fatality Rate"] = ifr_lic[3:9]
table_s0_lic <- table_s0_lic[,c("Age Group","Population Share","Infection Fatality Rate",as.character(rev(nl_q_seq)))]
print(xtable(table_s0_lic, digits=2), include.rownames=FALSE)

table_s1_lic<-data.frame(nlopt_s1[2:8,as.character(rev(nl_q_seq*prop_adults_lic))])
colnames(table_s1_lic) <- as.character(rev(nl_q_seq))
table_s1_lic["Age Group"] = c("20-30","30-40","40-50","50-60","60-70","70-80","80+")
table_s1_lic["Population Share"] = lic_pop[3:9]/sum(lic_pop[3:9])
table_s1_lic["Infection Fatality Rate"] = ifr_lic[3:9]
table_s1_lic <- table_s1_lic[,c("Age Group","Population Share","Infection Fatality Rate",as.character(rev(nl_q_seq)))]
print(xtable(table_s1_lic, digits=2), include.rownames=FALSE)

#Comparison with non-epi results----
# table_s0_lic <- rbind(table_s0_lic,c('res_nlopt',1,1,round(rev(as.numeric(nlopt_s0[1,])),8)))
# test_res <- function(x) model_fd_static(unroll_x(x,sub=0),
#                                       homogen = 1,
#                                       outcome = 'D',
#                                       ret = 1)
# 
# lic_nonepi <- read.csv('lic_non_epi.csv')
# res_nonepi <- c()
# for (i in 1:dim(lic_nonepi)[2]){
#   res_nonepi <- c(res_nonepi,test_res(lic_nonepi[,i]))
# }
# table_s0_lic <- rbind(table_s0_lic,c('res_nonepi',1,1,round(as.numeric(res_nonepi),8)))
##--

#HIC
prop_adults_hic <- sum(hic_pop[3:9])/sum(hic_pop)
load("results/nlopt/wip-nl-solutions7-D-hic-edit_contacts.Rdata")
table_s0_hic<-data.frame(nlopt_s0[2:8,as.character(rev(nl_q_seq*prop_adults_hic))])
colnames(table_s0_hic) <- as.character(rev(nl_q_seq))
table_s0_hic["Age Group"] = c("20-30","30-40","40-50","50-60","60-70","70-80","80+")
table_s0_hic["Population Share"] = hic_pop[3:9]/sum(hic_pop[3:9])
table_s0_hic["Infection Fatality Rate"] = ifr_hic[3:9]
table_s0_hic <- table_s0_hic[,c("Age Group","Population Share","Infection Fatality Rate",as.character(rev(nl_q_seq)))]
print(xtable(table_s0_hic, digits=2), include.rownames=FALSE)

table_s1_hic<-data.frame(nlopt_s1[2:8,as.character(rev(nl_q_seq*prop_adults_hic))])
colnames(table_s1_hic) <- as.character(rev(nl_q_seq))
table_s1_hic["Age Group"] = c("20-30","30-40","40-50","50-60","60-70","70-80","80+")
table_s1_hic["Population Share"] = hic_pop[3:9]/sum(hic_pop[3:9])
table_s1_hic["Infection Fatality Rate"] = ifr_hic[3:9]
table_s1_hic <- table_s1_hic[,c("Age Group","Population Share","Infection Fatality Rate",as.character(rev(nl_q_seq)))]
print(xtable(table_s1_hic, digits=2), include.rownames=FALSE)

#Comparison with non-epi results----
# table_s0_hic <- rbind(table_s0_hic,c('res_nlopt',1,1,round(rev(as.numeric(nlopt_s0[1,])),8)))
# 
# hic_nonepi <- read.csv('hic_non_epi.csv')
# res_nonepi <- c()
# for (i in 1:dim(hic_nonepi)[2]){
#   res_nonepi <- c(res_nonepi,test_res(hic_nonepi[,i]))
# }
# table_s0_hic <- rbind(table_s0_hic,c('res_nonepi',1,1,round(as.numeric(res_nonepi),8)))
##--