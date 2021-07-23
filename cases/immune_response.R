library(tidyverse)
library(ggplot2)
library(shades)
library(stats)
theme_set(theme_minimal(base_size = 10))

curve <- read.csv('data/curve.csv',header = FALSE)
colnames(curve) <- c("x","y")

points <- read.csv('data/doses.csv')

points['type'] <- 'Alternative doses'
points[points$dose==1,'type'] <- 'Status quo dose'

points <- points %>% mutate(dose=factor(dose))

#Points for variant factor
model.eff <- approxfun(curve$x,curve$y)
points['delta_efficacy'] <- model.eff(points$neutralization/5.8)
points[points$vaccine=='Convalescent','delta_efficacy'] <- NA

model.neutr <- approxfun(curve$y,curve$x)
model.neutr(92)/model.neutr(79)
model.neutr(93)/model.neutr(88)

frac_nab.fig <- ggplot() +
  geom_line(data = curve, aes(x=x,y=y), size=0.8) +
  geom_point(data=points,# %>% filter(age %in% c("18-55","12-65","18-59","All")),
             aes(x=neutralization,y=efficacy,color=vaccine,shape=type,size=dose)) +
  lightness(scale_color_manual(values = c("purple","deepskyblue3","brown","blue","deeppink","green4","orange","grey50","red"),name="Manufacturer"),scalefac(0.99))+
  scale_shape_manual(name="",values=c(16, 10)) +
  scale_size_discrete(range=c(3.5,6)) +
  ylab("Protective efficacy (%)") + xlab("Mean neutralization level (fold of convalescent)") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(), legend.text = element_text(size=10)) +
  guides(color = FALSE, size = FALSE, shape=guide_legend(override.aes = list(size = 3))) +
  scale_y_continuous(breaks = seq(35, 100, by = 5)) +
  scale_x_log10()
frac_nab.fig
ggsave(paste0("figures", "/fraction_imm_response_stddose.png"), frac_nab.fig, width = 6.5, height=6.5)

##Variant factor----
frac_nab.fig <- ggplot() +
  geom_line(data = curve, aes(x=x,y=y), size=0.8) +
  geom_point(data=points %>% filter(dose==1),
             aes(x=neutralization,y=efficacy,color=vaccine,size=dose),shape=10) +
  geom_point(data=points %>% filter(dose==1),
             aes(x=neutralization/5.8,y=delta_efficacy,color=vaccine,size=dose),shape=16) +
  lightness(scale_color_manual(values = c("purple","deepskyblue3","brown","blue","deeppink","green4","orange","grey50","red"),name="Manufacturer"),scalefac(0.99))+
  scale_size_discrete(range=c(3.5,6)) +
  ylab("Protective efficacy (%)") + xlab("Mean neutralization level (fold of convalescent)") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(), legend.text = element_text(size=10)) +
  guides(size = FALSE, shape=guide_legend(override.aes = list(size = 3))) +
  scale_y_continuous(breaks = seq(35, 100, by = 5)) +
  scale_x_log10()
frac_nab.fig

##With confidence intervals----
lb <- read.csv('data/lower_bound.csv',header = FALSE)
colnames(lb) <- c("x","y")
ub <- read.csv('data/upper_bound.csv',header = FALSE)
colnames(ub) <- c("x","y")

frac_nab.fig <- ggplot() +
  geom_line(data = curve, aes(x=x,y=y), size=0.8) +
  geom_line(data = lb, aes(x=x, y=y), fill = "grey70") +
  geom_line(data = ub, aes(x=x, y=y), fill = "grey70") +
  geom_point(data=points %>% filter(age %in% c("18-55","12-65","18-59","All")),
             aes(x=neutralization,y=efficacy,color=vaccine,shape=type,size=dose)) +
  # geom_point(data=points %>% filter(age %in% c("18-55","12-65","18-59","All")) %>% filter(dose!=1),
  #            aes(x=neutralization,y=efficacy,color=vaccine,size=dose),shape=16) +
  # lightness(scale_color_brewer(palette = "Paired",name="Manufacturer"),scalefac(0.99))+
  lightness(scale_color_manual(values = c("purple","deepskyblue3","brown","blue","deeppink","green4","orange","grey50","red"),name="Manufacturer"),scalefac(0.99))+
  scale_shape_manual(name="",values=c(16, 10)) +
  scale_size_discrete(range=c(3.5,6)) +
  ylab("Protective efficacy (%)") + xlab("Mean neutralization level (fold of convalescent)") +
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(), legend.text = element_text(size=10)) +
  guides(color = FALSE, size = FALSE, shape=guide_legend(override.aes = list(size = 3))) +
  scale_y_continuous(breaks = seq(35, 100, by = 5)) +
  scale_x_log10()
frac_nab.fig
