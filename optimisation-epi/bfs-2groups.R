source("project-setup.R")

w <- expand.grid(k1 = seq(0.1,1,0.05), k2 = seq(0.1,1,0.05))
brute_force_inc_vectors <- t(apply(w, 1, function(x) c(1,1,x[1],x[1],x[1],x[1],x[2],x[2],x[2])))
colnames(brute_force_inc_vectors) <- paste0("age", 1:9)


phi_x <- function(x) 
  sapply(3.49706*sqrt(x) - 1.74853*x  -0.798528, function(y) max(y,0))
# 3.49706*sqrt(x) - 1.74853*x  -0.798528

model_fd_dynamic <- function(model, d1, fd, default_e1 = 0.95, 
                             rm = FALSE,
                             homogen = FALSE) {
  
  e1 <- phi_x(fd)
  pars <- apap_2v(grab_2v_parms(model), fractional_dose = fd, len = d1)
  pars <- list_modify(pars, e1 = e1)
  if(homogen){
    pars$contacts <- 1/Ngroups + 0*pars$contacts
    pars$q <- ev*pars$q
  }
  y <- sr(pars, "2v_v2")
  if(rm) return(y)
  main_metrics(y, pop)
}


# Simulations for the dynamic problem -----
df_fd_dynamic <- expand.grid(model = scenario_par_nms_2v[3],
                             homogeneous = c(0,1),
                             # d1 = c(1000,400,200,100,50),
                             d1 = c(1000,750,500,400,300,200,133,100,50,25),
                             fd_v = 1:nrow(brute_force_inc_vectors)) %>%
  mutate(data = pmap(list(model, d1, fd_v,homogeneous), 
                     function(x,y,z,h) data.frame(value = model_fd_dynamic(x,y,brute_force_inc_vectors[z,], homogen=h), 
                                                  var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 



res_dynamic <- cbind(df_fd_dynamic, brute_force_inc_vectors[df_fd_dynamic$fd_v,]) %>%
  # filter(!homogeneous) %>%
  select(age3, age9, model, d1, d, i,homogeneous) %>%
  gather(variable, value, -age3,-age9,-model,-d1,-homogeneous) %>%
  group_by(model, d1,variable,homogeneous) %>%
  slice_min(value) %>%
  arrange(homogeneous)
# summarise(optimal_v = paste(age3[which.min(value)], age9[which.min(value)]))
res_dynamic %>%
  filter(d1 %in% c(1000,400,200,100,50))

cbind(df_fd_dynamic, brute_force_inc_vectors[df_fd_dynamic$fd_v,]) %>%
  filter(d1 %in% c(1000,400,200,100,50)) %>%
  filter(homogeneous == 0) %>%
  select(age3, age9, model, d1, d, i) %>%
  gather(variable, value, -age3,-age9,-model,-d1) %>%
  group_by(model, d1,variable) %>%
  mutate(value = value/max(value)) %>%
  mutate(value = value - .1*(value == min(value))) %>%
  mutate(variable = factor(variable, levels = c("i", "d", "harm"),
                           labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(d1 = paste0(100/d1, "% / day")) %>%
  ggplot(aes(x = age3, y = age9, fill = value)) + geom_tile() + facet_wrap(variable~d1, ncol = 5) +
  scale_fill_viridis_c(name = "Objective function: relative benefit (normalised by max harm in each panel)") +
  xlab("Dosing in 20-60 year olds (1 is full dose)") +
  ylab("Dosing in 60+ year olds (1 is full dose)") +
  theme(legend.position = "top") +
  labs(subtitle="Shaded tile is the optimal solution (but there is no discontinuity, it's just an extra visual flair)")





# Simulations for the static problem -----

prop_young <- sum(pop[1:2])
prop_old <- sum(pop[7:9])
prop_work <- 1-prop_old-prop_young
v <- c()
for(Q in seq(0.1, 1, .05)) 
  for(w1 in seq(0, min(prop_old, Q), .005))
    v <- rbind(v, c(Q, (Q-w1)/prop_work, w1/prop_old))
v <- v[v[,3] <= 1,]
# v <- v[v[,3] >= .1,]
v <- v[v[,2] <= 1,]
# v <- v[v[,2] >= .1,]
colnames(v) <- c("Q", "age3", "age9")
plot(v[,2] ~ v[,3])

model_fractional_static <- function(v_prop, rm = FALSE, homogen = FALSE) {
  e_vector <- phi_x(v_prop)
  pars <- list_modify(
    pars_le_fast,
    y0 = y0_gen(13, Ngroups, pre_immunity = pre_immunity + (1-pre_immunity)*e_vector))
  
  if(homogen){
    pars$contacts <- 1/Ngroups + 0*pars$contacts
    pars$q <- ev*pars$q
  }
  y <- sr(pars, "2v_v2")
  main_metrics(y, pop)[1:2]
}

unroll_x <- function(x) c(0,0,rep(x[2],4),rep(x[3],3))
df_fd_static_h0 <- cbind(v, t(apply(v, 1, function(x) model_fractional_static(unroll_x(x), homogen = FALSE))))
df_fd_static_h1 <- cbind(v, t(apply(v, 1, function(x) model_fractional_static(unroll_x(x), homogen = TRUE))))
df_fd_static <- rbind(as.data.frame(df_fd_static_h0) %>% mutate(homogeneous = 0),
                      as.data.frame(df_fd_static_h1) %>% mutate(homogeneous = 1))

res_static <- rbind(
  df_fd_static %>% group_by(Q, homogeneous) %>% slice_min(i) %>% mutate(variable = "i") %>% select(-i, -d),
  df_fd_static %>% group_by(Q, homogeneous) %>% slice_min(d) %>% mutate(variable = "d") %>% select(-i, -d)
)




rbind(
  res_static %>% rename(d1 = Q) %>% mutate(model_type = "static") %>% mutate(d1 = paste("Q =", d1)),
  res_dynamic %>%
    # mutate(optimal = paste(format(age3, nsmall=2), format(age9, nsmall=2))) %>%
    ungroup() %>%
    mutate(model_type = "dynamic") %>%
    # filter(d1 %in% c(100,200,400)) %>%
    mutate(d1 = paste("d1 =", d1)) %>%
    select(-model, -value)
) %>%
  mutate(variable = factor(variable, levels = c("i", "d", "harm"),
                           labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(mixing = factor(homogeneous, levels = c(0,1), labels = c("Heterogen.", "Homogen."))) %>%
  ggplot(aes(x = age3, y = age9, pch = mixing, color = model_type,
             group = interaction(mixing, model_type))) + 
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_wrap(~variable) +
  xlab("Dosing in 20-60 year olds (1 is full dose)") +
  ylab("Dosing in 60+ year olds (1 is full dose)")


save(res_static, res_dynamic, df_fd_dynamic, file = "results/explore-2groups.Rdata")
