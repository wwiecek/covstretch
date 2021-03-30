library(nloptr)
library(tidyverse)
source("project-setup.R")
load("results/wip-nl-solutions7-D.Rdata")



rbind(
  as.data.frame(nlopt_d0)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "dynamic",  mixing = "heterogeneous"),
  as.data.frame(nlopt_d1)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "dynamic",  mixing = "homogeneous"),
  as.data.frame(nlopt_s0)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "static", mixing = "heterogeneous"),
  as.data.frame(nlopt_s1)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "static", mixing = "homogeneous")
) %>%
  mutate(agegr = factor(age, levels = paste0(1:9), labels = colnames(pbc_spread))) %>%
  mutate(mixing = factor(mixing)) %>%
  ggplot(aes(x = agegr, y= value, group = interaction(mixing, s), lty = mixing, color = mixing)) + 
  facet_grid(model ~ .) +
  geom_line(size = 1.5)
