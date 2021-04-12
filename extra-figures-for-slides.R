source("project-setup.R")
source("cases/prep-results.R")
source("cases/general-example.R")

gg1.df %>%
  # filter(var %in% c("No vaccination", "Vaccinate 1% per day")) %>%
  filter(var != "Vaccinate 0.5% per day") %>%
  ggplot(aes(x = time, y = value, color = var)) + geom_line() + 
  facet_wrap(scenario ~ ., scales = "free", ncol = 1) +
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + 
  ylab("fraction infected, I(t)") +
  scale_color_discrete(name = "") +
  theme(legend.position = "top", legend.direction = "vertical")

ggsave(width = 2.5, height = 6, file = "figures/presentation-scenarios.pdf")
