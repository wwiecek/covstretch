# Compare the results with brute force approach -----
load("results/explore-2groups.Rdata")
load("results/wip-nl-solutions.Rdata")

q_seq <- rep(seq(0.1, 0.7, 0.1), each = 5)

res_static_nl <- rbind(
  data.frame(Q = q_seq, age3 = h0[2,], age9 = h0[3,], homogeneous = 0, variable = "i"),
  data.frame(Q = q_seq, age3 = h1[2,], age9 = h1[3,], homogeneous = 1, variable = "i")
)

rbind(
  res_static %>% mutate(opt = "brute force") %>% filter(variable == "i"),
  res_static_nl %>% mutate(opt = "nloptr")
) %>%
  filter(Q %in% q_seq) %>%
  ggplot(aes(x = age3, y = age9, color = opt, pch = factor(Q))) + 
  geom_point(size = 2.5) + facet_grid(homogeneous ~ .)

