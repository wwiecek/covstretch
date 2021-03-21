# Compare the results with brute force approach -----
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

