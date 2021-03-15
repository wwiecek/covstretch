
model_fractional_d <- function(model, d1, alpha, fd, K, default_e1 = 0.95, rm = FALSE) {
  
  # Prepare parameters, taking fractional dosing into account:
  dose_response <- default_e1*fd^alpha
  N <- length(pop)
  
  if(K <= N){
    fd <- c(rep(fd, K), rep(1, N-K))
    e1 <- c(rep(dose_response, K), rep(default_e1, N-K))
  } else {
    fd <- rep(1, N)
    e1 <- rep(default_e1, N)
  }
  pars <- apap_2v(grab_2v_parms(model), fractional_dose = fd, len = d1)
  # Done.
  
  y <- sr(list_modify(pars, e1 = e1), "2v_v2")
  if(rm) return(y)
  main_metrics(y, pop)
}

df_fractional <- expand.grid(model = scenario_par_nms_2v,
                             d1 = c(d1_general, Inf),
                             alpha = c(0, 1, 0.25),
                             fractional_dose = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.75),
                             K = 3:10) %>%
  mutate(data = pmap(list(model, d1, alpha, fractional_dose, K), 
                     function(x,y,z,a,b) data.frame(value = model_fractional_d(x,y,z,a,b), 
                                                    var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 


df_fractional %>% select(model, d1, alpha, fractional_dose, K, d) %>%
  filter(d1 == 100) %>%
  mutate(K = factor(K)) %>%
  # spread(d1, i) %>% 
  ggplot(aes(x = fractional_dose, y = d, group = K, color = K)) + 
  geom_line() + facet_wrap(model ~ alpha, scales = "free")

# Find the optimal fractional dosing policy for each speed, scenario, dose response shape
df_fractional %>% select(model, d1, alpha, fractional_dose, K, d) %>%
  rename(i = d) %>%
  group_by(model, d1, alpha) %>%
  summarise(v = paste(K[which.min(i)], fractional_dose[which.min(i)])) %>%
  # summarise(v = min(i)/max(i)) %>%
  spread(d1, v)
# filter(alpha == 1)

filter(d1 == 100) %>%
  mutate(K = factor(K)) %>%
  # spread(d1, i) %>% 
  ggplot(aes(x = fractional_dose, y = i, group = K, color = K)) + geom_line() + facet_wrap(model ~ alpha, scales = "free")




# Brute force search for optimal solution -----

default_e1 <- 0.95
x <- seq(0.1, 1, length = 100)
y <- default_e1*x^0.25
plot(y~x, type = "l", ylim = c(0.5,1))
y <- default_e1*x^0.5
lines(y~x, type = "l")
y <- default_e1*x^0.75
lines(y~x, type = "l")


# m <- seq(0.1, 1.0, 0.1)
# w <- c()
# for(k1 in m) {
#   for(k2 in seq(k1, 1, 0.1)) {
#     for(k3 in seq(k2, 1, 0.1)) {
#       for(k4 in seq(k3, 1, 0.1)) {
#         for(k5 in seq(k4, 1, 0.1)) {
#           for(k6 in seq(k5, 1, 0.1)) {
#             w <- rbind(w, c(k1,k2,k3,k4,k5,k6))
#           }
#         }
#       }
#     }
#   }
# }
# brute_force_inc_vectors <- cbind(0, 0, w, 1)
# rm(w); rm(m)

step <- 0.1
m <- seq(0.1, 1.0, step)
w <- c()
for(k1 in m) {
  for(k2 in seq(k1, 1, step)) {
    for(k3 in seq(k2, 1, step)) {
      for(k4 in seq(k3, 1, step))
        w <- rbind(w, c(k1,k1,k2,k2,k3,k3,k4))
    }
  }
}
brute_force_inc_vectors <- cbind(1, 1, w) #by convention for under 20s f.d. is 1, not 0
                                          #otherwise division by 0
rm(w); rm(m); rm(step)
brute_force_inc_vectors <- brute_force_inc_vectors[brute_force_inc_vectors[,9] >= 0.5, ]

model_fractional_d <- function(model, d1, alpha, fd, K, default_e1 = 0.95, rm = FALSE) {
  
  # Prepare parameters, taking fractional dosing into account:
  e1 <- default_e1*fd^alpha
  pars <- apap_2v(grab_2v_parms(model), fractional_dose = fd, len = d1)
  # Done.
  y <- sr(list_modify(pars, e1 = e1), "2v_v2")
  if(rm) return(y)
  main_metrics(y, pop)
}

df_fractional <- expand.grid(model = scenario_par_nms_2v,
                             d1 = c(d1_general, Inf),
                             alpha = c(0.25, 0.5, 1.0),
                             fd_v = 1:nrow(brute_force_inc_vectors)) %>%
  mutate(data = pmap(list(model, d1, alpha, fd_v), 
                     function(x,y,z,a) data.frame(value = model_fractional_d(x,y,z,brute_force_inc_vectors[a,],b), 
                                                    var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 
save(df_fractional, file = "fractional_wip.Rdata")

bfv <- apply(format(brute_force_inc_vectors[,c(3,5,7,9)], nsmall = 1), 1, paste, collapse = ", ")

df_fractional %>% select(model, d1, alpha, fd_v, d) %>%
  rename(i = d) %>%
  group_by(model, d1, alpha) %>%
  summarise(bfv[fd_v[which.min(i)]], min(i))


df_fractional %>% select(model, d1, alpha, fd_v, d, i) %>%
  mutate(d1 = 1/d1) %>% 
  filter(d1 %in% c(0.001, 0.005, 0.02)) %>%
  gather(variable, value, -model, -d1, -fd_v, -alpha) %>%
  mutate(variable = factor(variable, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  group_by(model, d1, variable, alpha) %>%
  summarise(v = bfv[fd_v[which.min(value)]]) %>%
  spread(d1, v) %>%
  arrange(variable, model, alpha)
  
  
  
  mutate(delta1 = 1/d1) %>% 
  select(delta1, e, model, i,d) %>%
  gather(var, value, -delta1, -e, -model) %>%
  group_by(model,var) %>%
  mutate(ref = value[e == .95 & delta1 == default_delta_value]) %>%
  filter(e %in% seq(.5, .9, .1)) %>%
  ungroup() %>%
  mutate(r = (value/ref)) %>% 
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Less effective better by 5% or more",
                                    
                                    "Comparable (+-5%)", 
                                    "95% effective better by at least 5%"
                         ))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(speedup = factor(round(delta1/default_delta_value, 1))) %>%
  # mutate(speedup = factor(as.percent(delta1, 2))) %>%
  mutate(e = factor(e)) %>%
  filter(var != "Economic harm") %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(var~model) + ylab("e2 (efficacy for the less effective vaccine)") + 
  xlab("delta2/delta1 (speed-up factor vs 0.25% base case)") + 
  # scale_x_continuous(breaks = 1/d1_general,
  # labels = as.percent(1/d1_general))
  geom_text(aes(label = value), color = "white", size = 2.5)
  
