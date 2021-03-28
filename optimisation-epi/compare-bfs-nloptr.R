# Compare the nlopt with brute force approach - 2 groups -----
load("results/explore-2groups.Rdata")
load("results/wip-nl-solutions.Rdata")

res_bfs <- rbind(
  res_static %>% rename(d1 = Q) %>% mutate(model_type = "static") ,
  res_dynamic %>%
    # mutate(optimal = paste(format(age3, nsmall=2), format(age9, nsmall=2))) %>%
    ungroup() %>%
    mutate(model_type = "dynamic") %>%
    select(-model, -value)
) %>% mutate(opt = "brute force")

res_nlopt <- rbind(
  data.frame(d1 = nl_q_seq, age3 = nlopt_s0[2,], age9 = nlopt_s0[3,], homogeneous = 0, model_type = "static",  variable = "i"),
  data.frame(d1 = nl_q_seq, age3 = nlopt_s1[2,], age9 = nlopt_s1[3,], homogeneous = 1, model_type = "static",  variable = "i"),
  data.frame(d1 = nl_d_seq, age3 = nlopt_d0[2,], age9 = nlopt_d0[3,], homogeneous = 0, model_type = "dynamic", variable = "i"),
  data.frame(d1 = nl_d_seq, age3 = nlopt_d1[2,], age9 = nlopt_d1[3,], homogeneous = 1, model_type = "dynamic", variable = "i")
) %>% mutate(opt = "nloptr")

rbind(
  res_bfs  %>% filter(variable == "i"),
  res_nlopt 
) %>%
  filter(model_type == "dynamic") %>%
  # filter(d1 %in% seq(0.2, 0.7, 0.1)) %>%
  filter(d1 %in% c(1000,400,200,100,50)) %>%
  mutate(d1 = factor(d1)) %>%
  mutate(mixing = factor(homogeneous, levels = c(0,1), labels = c("Heterogen.", "Homogen."))) %>%
  ggplot(aes(x = age3, y = age9, color = opt, pch = factor(d1))) + 
  geom_point(size = 2.5) + facet_grid(mixing ~ .)

# Comparison of nlopt and bfs - 7 groups ----
load("results/opt-bfs-result7.Rdata")
#Change the file below, the 'variable' filter, and the outcome in model_fd_static for deaths/infections
load("results/wip-nl-solutions7-I.Rdata")

#Static case

unroll_x <- function(x,sub=1) c(sub,sub,x[1],x[2],x[3],x[4],x[5],x[6],x[7])

eval_g_ineq <- function(x) c(x[1]*pop[3]+x[2]*pop[4]+x[3]*pop[5]+x[4]*pop[6]+x[5]*pop[7]+x[6]*pop[8]+x[7]*pop[9]-Q,
                             -x[1]*pop[3]-x[2]*pop[4]-x[3]*pop[5]-x[4]*pop[6]-x[5]*pop[7]-x[6]*pop[8]-x[7]*pop[9]+Q)

df_comp <- data.frame()
for (Q in unique(nl_q_seq)){
  res_bfs_0 <- as.numeric(res_static[(res_static$q==round(Q,2))&(res_static$homogeneous==0)&(res_static$variable=='i')
                                   ,c('age3','age4','age5','age6','age7','age8','age9')])
  res_bfs_1 <- as.numeric(res_static[(res_static$q==round(Q,2))&(res_static$homogeneous==1)&(res_static$variable=='i')
                                   ,c('age3','age4','age5','age6','age7','age8','age9')])
  eval_bfs_0 <- model_fd_static(unroll_x(res_bfs_0,sub=0), 
                                homogen = 0, 
                                outcome = "cumI",
                                ret = 1)
  eval_bfs_1 <- model_fd_static(unroll_x(res_bfs_1,sub=0), 
                                homogen = 1, 
                                outcome = "cumI",
                                ret = 1)
  eval_ineq_bfs_0 <- eval_g_ineq(res_bfs_0)[1]
  eval_ineq_bfs_1 <- eval_g_ineq(res_bfs_1)[1]
  
  eval_nlopt_0 <- nlopt_s0[1,as.character(Q)]
  res_nlopt_0 <- nlopt_s0[2:8,as.character(Q)]
  eval_ineq_nlopt_0 <- eval_g_ineq(res_nlopt_0)[1]
  eval_nlopt_1 <- nlopt_s1[1,as.character(Q)]
  res_nlopt_1 <- nlopt_s1[2:8,as.character(Q)]
  eval_ineq_nlopt_1 <- eval_g_ineq(res_nlopt_1)[1]
  
  df_comp <- rbind(df_comp, data.frame("eval_bfs_0"=eval_bfs_0,"eval_nlopt_0"=eval_nlopt_0,"constr_bfs_0"=eval_ineq_bfs_0,
                                       "constr_nlopt_0"=eval_ineq_nlopt_0,"eval_bfs_1"=eval_bfs_1,"eval_nlopt_1"=eval_nlopt_1,
                                       "constr_bfs_1"=eval_ineq_bfs_1,"constr_nlopt_1"=eval_ineq_nlopt_1,"q"=Q))
}
#Comparing objective function and checking constraints for nlopt - h=0
print (dim(df_comp %>% filter(eval_nlopt_0<=eval_bfs_0))[1]/dim(df_comp)[1])
print (dim(df_comp %>% filter(eval_ineq_nlopt_0<=0.0001))[1]/dim(df_comp)[1])
#Comparing objective function and checking constraints for nlopt - h=1
print (dim(df_comp %>% filter(eval_nlopt_1<=eval_bfs_1))[1]/dim(df_comp)[1])
print (dim(df_comp %>% filter(eval_ineq_nlopt_1<=0.0001))[1]/dim(df_comp)[1])

#Dynamic case

df_comp <- data.frame()
for (D in unique(nl_d_seq)){
  res_bfs_0 <- as.numeric(res_dynamic[(res_dynamic$d1==D)&(res_dynamic$homogeneous==0)&(res_dynamic$variable=='i')
                                    ,c('age3','age4','age5','age6','age7','age8','age9')])
  eval_bfs_0 <- as.numeric(res_dynamic[(res_dynamic$d1==D)&(res_dynamic$homogeneous==0)&(res_dynamic$variable=='i')
                                     ,'value'])
  res_bfs_1 <- as.numeric(res_dynamic[(res_dynamic$d1==D)&(res_dynamic$homogeneous==1)&(res_dynamic$variable=='i')
                                    ,c('age3','age4','age5','age6','age7','age8','age9')])
  eval_bfs_1 <- as.numeric(res_dynamic[(res_dynamic$d1==D)&(res_dynamic$homogeneous==1)&(res_dynamic$variable=='i')
                                     ,'value'])
  
  eval_nlopt_0 <- nlopt_d0[1,as.character(D)]
  res_nlopt_0 <- nlopt_d0[2:8,as.character(D)]
  eval_nlopt_1 <- nlopt_d1[1,as.character(D)]
  res_nlopt_1 <- nlopt_d1[2:8,as.character(D)]
  
  df_comp <- rbind(df_comp, data.frame("eval_bfs_0"=eval_bfs_0,"eval_nlopt_0"=eval_nlopt_0,
                                       "eval_bfs_1"=eval_bfs_1,"eval_nlopt_1"=eval_nlopt_1,"d1"=D))
}
#Comparing objective function - h=0
print (dim(df_comp %>% filter(eval_nlopt_0<=eval_bfs_0))[1]/dim(df_comp)[1])
#Comparing objective function - h=0
print (dim(df_comp %>% filter(eval_nlopt_1<=eval_bfs_1))[1]/dim(df_comp)[1])

