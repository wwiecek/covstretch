library(nloptr)
library(tidyverse)
source("project-setup.R")
#The file below is used only for the dynamic case
load("results/nlopt/wip-nl-solutions7-D-hic.Rdata")
#Infection - load("results/nlopt/wip-nl-solutions7-I-maxeval100-tol-4.Rdata")
#The file below is used for the static case to be comparable with the non-epi table
load("results/nlopt/wip-nl-solutions7-s-D-hic.Rdata")

colnames(nlopt_s0) <- nl_q_seq
colnames(nlopt_s1) <- nl_q_seq

#Importing results from non-epi model

nonepi_s0 <- read.csv("model/hic_non_epi.csv")
colnames(nonepi_s0) <- rev(nl_q_seq)

# Single plot
rbind(
  as.data.frame(nlopt_d0)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "epi-dynamic",  mixing = "heterogeneous"),
  as.data.frame(nlopt_d1)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "epi-dynamic",  mixing = "homogeneous"),
  as.data.frame(nlopt_s0)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "epi-static", mixing = "heterogeneous"),
  as.data.frame(nlopt_s1)[-1,] %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "epi-static", mixing = "homogeneous"),
  as.data.frame(nonepi_s0) %>% mutate(age = 3:9) %>% 
    gather(s, value, -age) %>% 
    mutate(model = "non-epi-static", mixing = "NA")
) %>%
  filter(s %in% c("1000","400","100","0.3","0.6","0.9")) %>% 
  mutate(agegr = factor(age, levels = paste0(1:9), labels = colnames(pbc_spread))) %>%
  mutate(mixing = factor(mixing)) %>%
  ggplot(aes(x = agegr, y= value, group = interaction(mixing, s, model), lty = mixing, color = s)) + 
  facet_grid(model ~ .) +
  geom_line(size = 1.5)


#Multiple plots
plot1 <- 
  rbind(
  as.data.frame(nlopt_d0)[-1,] %>% mutate(age = 3:9) %>% 
  gather(s, value, -age) %>% 
  mutate(mixing = "heterogeneous"),
  as.data.frame(nlopt_d1)[-1,] %>% mutate(age = 3:9) %>% 
  gather(s, value, -age) %>% 
  mutate(mixing = "homogeneous")
  ) %>% 
  filter(s %in% c("1000","400","100")) %>% 
  mutate(delta1 = as.percent(round(1/as.numeric(s),4))) %>%
  mutate(agegr = factor(age, levels = paste0(1:9), labels = colnames(pbc_spread))) %>%
  mutate(mixing = factor(mixing)) %>%
  ggplot(aes(x = agegr, y= value, group = interaction(mixing, delta1), lty = mixing, color = delta1)) + 
  geom_line(size = 1.5) +
  xlab("Age group") + ylab("Dose fraction") + ggtitle("Dynamic Epidemiological - Homogeneous/Heterogeneous") +
  guides(color=guide_legend(title="vaccination speed (% per day)"))

plot2 <- 
  rbind(
    as.data.frame(nlopt_s0)[-1,] %>% mutate(age = 3:9) %>% 
      gather(Q, value, -age) %>% 
      mutate(model = "epi-static", mixing = "heterogeneous"),
    as.data.frame(nlopt_s1)[-1,] %>% mutate(age = 3:9) %>% 
      gather(Q, value, -age) %>% 
      mutate(model = "epi-static", mixing = "homogeneous")
  ) %>% 
  filter(Q %in% c("0.3","0.6","0.9")) %>% 
  mutate(agegr = factor(age, levels = paste0(1:9), labels = colnames(pbc_spread))) %>%
  mutate(mixing = factor(mixing)) %>%
  ggplot(aes(x = agegr, y= value, group = interaction(mixing, Q), lty = mixing, color = Q)) + 
  geom_line(size = 1.5) +
  xlab("Age group") + ylab("Dose fraction") + ggtitle("Static Epidemiological - Homogeneous/Heterogeneous") +
  guides(color=guide_legend(title="supply constraint (% above 20)"))

plot3 <- 
  rbind(
    as.data.frame(nlopt_s0)[-1,] %>% mutate(age = 3:9) %>% 
      gather(Q, value, -age) %>% 
      mutate(model = "epidemiological (homogeneous)"),
    as.data.frame(nonepi_s0) %>% mutate(age = 3:9) %>% 
      gather(Q, value, -age) %>% 
      mutate(model = "theoretical")
  ) %>% 
  filter(Q %in% c("0.3","0.6","0.9")) %>% 
  mutate(agegr = factor(age, levels = paste0(1:9), labels = colnames(pbc_spread))) %>%
  mutate(model = factor(model)) %>%
  ggplot(aes(x = agegr, y= value, group = interaction(model, Q), lty = model, color = Q)) + 
  geom_line(size = 1.5) +
  xlab("Age group") + ylab("Dose fraction") + ggtitle("Static Models - Epidemiological/Theoretical") +
  guides(color=guide_legend(title="supply constraint (% above 20)"))

fig_dose_sharing <-
  ggarrange(plotlist=list(plot1,plot2,plot3), common.legend = FALSE, ncol = 1, legend = "right")
ggsave("figures/dose_sharing_models.pdf", fig_dose_sharing+theme(text = element_text(size=3)), width = 6.5, height=7.4/5.55*6.5)
