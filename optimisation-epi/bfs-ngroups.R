source("project-setup.R")

# Set up N groups and search grid -----
n_x <- 3
search_grid_dynamic <- seq(0.1, 1, 0.1)
unroll_x <- function(x,sub=1) c(sub,sub,x[1],x[1],x[2],x[2],x[3],x[3],x[3])
d1_opt <- c(1000,400,200,100,50) #dynamic model speeds
# d1_opt <- c(1000,750,500,400,300,200,133,100,50,25) #dynamic model speeds

# Generate all vectors to search over -----
k <- lapply(as.list(1:n_x), function(i) search_grid_dynamic)
names(k) <- paste0("x", 1:n_x)
w <- do.call(expand.grid, k)
brute_force_inc_vectors <- t(apply(w, 1, unroll_x)) %>% unique()
colnames(brute_force_inc_vectors) <- paste0("age", 1:9)


# Define the dose-response function -----
phi_x <- function(x) 
  # 3.49706*sqrt(x) - 1.74853*x  -0.798528
  # sapply(3.49706*sqrt(x) - 1.74853*x  -0.798528, function(y) max(y,0))
  -25.31701*x^1.037524 + 1.037524*25.31701*x


model_fd_dynamic("pars_le_fast", 100, c(0.8, 0.6, 0.8), h = 1)

# Define objective functions (static and dynamic cases) ------
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
df_fd_dynamic <- expand.grid(model = "pars_le_fast",
                             homogeneous = c(0,1),
                             # d1 = c(1000,400,200,100,50),
                             d1 = d1_opt,
                             fd_v = 1:nrow(brute_force_inc_vectors)) %>%
  mutate(data = pmap(list(model, d1, fd_v,homogeneous), 
                     function(x,y,z,h) data.frame(value = model_fd_dynamic(x,y,brute_force_inc_vectors[z,], homogen=h), 
                                                  var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 



res_dynamic <- cbind(df_fd_dynamic, brute_force_inc_vectors[df_fd_dynamic$fd_v,]) %>%
  select(age3, age5, age9, model, d1, d, i,homogeneous) %>%
  gather(variable, value, -age3,-age5,-age9,-model,-d1,-homogeneous) %>%
  group_by(model, d1,variable,homogeneous) %>%
  slice_min(value) %>%
  arrange(homogeneous) %>%
  ungroup()

res_dynamic %>%
  filter(d1 %in% c(1000,400,200,100,50)) %>%
  arrange(model, homogeneous, variable, d1) %>%
  mutate(opt = paste(age3, age5, age9))

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
search_grid_static <- seq(0.1, 1, 0.025)

k_s <- lapply(as.list(1:n_x), function(i) search_grid_static)
names(k_s) <- paste0("x", 1:n_x)
w_s <- do.call(expand.grid, k_s)
brute_force_inc_vectors_s <- t(apply(w_s, 1, unroll_x, sub = 0)) %>% unique()
colnames(brute_force_inc_vectors_s) <- paste0("age", 1:9)
qm <- t(apply(brute_force_inc_vectors_s, 1, function(x) {
  # Find the closest Q on the list:
  q <- sum(x*pop)
  qs <- floor(10*q)/10
  round(x*qs/q, 5)
})) %>% unique()
q <- round(apply(qm, 1, function(x) sum(x*pop)), 3)

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

df_fd_static_h0 <- cbind(q, qm, t(apply(qm, 1, function(x) model_fractional_static(unroll_x(x), homogen = FALSE))))
df_fd_static_h1 <- cbind(q, qm, t(apply(qm, 1, function(x) model_fractional_static(unroll_x(x), homogen = TRUE))))
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


save(res_static, res_dynamic, df_fd_dynamic, file = paste0("results/opt-bfs", n_x, ".Rdata"))
save.image(file = paste0("results/opt-bfs", n_x, "-all.Rdata"))
