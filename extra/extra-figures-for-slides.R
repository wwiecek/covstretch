source("project-setup.R")
source("cases/prep-results.R")
source("cases/general-example.R")

##Slide Epidemic and Vaccination Scenarios----

gg1.presentation <- gg1.df %>%
  # filter(var %in% c("No vaccination", "Vaccinate 1% per day")) %>%
  filter(var != "Vaccinate 0.5% per day") %>%
  ggplot(aes(x = time, y = value, color = var)) + geom_line() + 
  facet_wrap(scenario ~ ., scales = "free", ncol = 3) +
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + 
  ylab("fraction infected, I(t)") +
  scale_color_discrete(name = "") +
  theme(legend.position = "top")

ggsave(file = "figures/presentation-scenarios.png",gg1.presentation+theme(text = element_text(size=14)),width = 7.5, height = 3.5)

##Slide Impact of Speed on Health Benefits From Vaccination----

benefits_curve_scenario <- function(scenario,e=0.95){
  params <- scenario_list_2v[[scenario]]
  ii <- 1e-04
  sapply(
    seq(0,1,length=21),
    function(p) {
      v <- vac_top_p(p, pop)
      y0_v <- y0_gen(13, Ngroups, pre_immunity = pre_immunity + (1-pre_immunity)*e*v)
      sr <- sr(list_modify(params, y0 = y0_v))
      c(p = p,
        bd = bd(sr, pop),
        bi = b_any(sr, pop, "cumI"))
    }) %>% t() %>% as.data.frame()
}

bf_scenarios <- rbind(
  benefits_curve_scenario("Slow decrease") %>% mutate(scenario = "Slow decrease"),
  benefits_curve_scenario("Slow growth") %>% mutate(scenario = "Slow growth"),
  benefits_curve_scenario("Fast growth") %>% mutate(scenario = "Fast growth")) 

bf_scenarios <-
  bf_scenarios %>%
  group_by(scenario) %>%
  mutate(bi = 1-(bi/max(bi))) %>%
  mutate(bd = 1-(bd/max(bd))) %>%
  ungroup()

benefits_gg_scenarios <- bf_scenarios %>% 
  select(p,bi,bd,scenario) %>% 
  setNames(c("p", "Infections", "Deaths", "scenario")) %>%
  gather(var, value, -p, -scenario) %>%
  arrange(var) %>% 
  mutate(scenario=factor(scenario,levels=c("Slow decrease", "Slow growth", "Fast growth"))) %>% 
  mutate(var=factor(var,levels=c("Infections", "Deaths"))) %>% 
  ggplot(aes(x = p, y = value, group = scenario, color = scenario, lty = scenario)) + 
  geom_line() + 
  # scale_color_viridis_d() +
  facet_wrap(~var, scales = "free") + 
  xlab("Fraction vaccinated before epidemic") + 
  ylab("Fraction of harm averted") +
  ylim(0, 1) +
  theme(legend.position = "right", legend.title = element_blank(),text=element_text(size=9))

g2b <- df_efficacy_delta  %>%
  filter(e == .95) %>%
  filter(round(d1,4) %in% round(d1_general,4)) %>%
  select(ri, rd) %>%
  gather(key, value, -d1, -model, -e) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(delta1 = 1/d1) %>%
  ggplot(aes(x = delta1, y = value, group = model, color = model)) + 
  geom_line(size=1.1) +
  geom_point(pch = 21, size = 2, fill = "white") +
  facet_wrap(~key, ncol = 3) +
  scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
  scale_color_discrete(name = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, size = 11), legend.position = "top",legend.direction = "vertical") +
  xlab(def_labels$speed) + 
  ylab("Fraction of harm averted")

fig_harm_averted<-ggarrange(benefits_gg_scenarios+ggtitle("Static (vaccination pre-epidemic)")+theme(text=element_text(size=12),
                                                                                                     plot.title = element_text(size=10)), 
                  g2b+ggtitle("Dynamic")+theme(text=element_text(size=12),
                                               plot.title = element_text(size=10)), 
                  common.legend = TRUE, ncol = 1, 
                  heights = c(2.6,3), legend = "bottom")

fig_harm_averted

ggsave(fig_harm_averted, width = 6.5, height = 5, file = "figures/presentation-harm-averted-speed.png")

ggsave(file = "figures/presentation-g2_reductions.png", g2b+theme(text=element_text(size=14)), width = 6.5, height=4.5)

##Slide Fractional Dosing - Optimal Allocation----
load("results/nlopt/wip-nl-solutions7-D-hic-edit_contacts.Rdata")

colnames(nlopt_s0) <- nl_q_seq
colnames(nlopt_s1) <- nl_q_seq

#Importing results from non-epi model

nonepi_s1 <- read.csv("model/hic_non_epi.csv")
colnames(nonepi_s1) <- rev(nl_q_seq)
nonepi_s1[nonepi_s1<1e-10] <- 0

nonepi_s1_lic <- read.csv("model/lic_non_epi.csv")
colnames(nonepi_s1_lic) <- rev(nl_q_seq)

fig_static <-
  rbind(
    as.data.frame(nlopt_s0)[-1,] %>% mutate(age = 3:9) %>% 
      gather(Q, value, -age) %>% 
      mutate(model = "Heterogeneous mixing"),
    as.data.frame(nlopt_s1)[-1,] %>% mutate(age = 3:9) %>% 
      gather(Q, value, -age) %>% 
      mutate(model = "Homogeneous mixing"),
    as.data.frame(nonepi_s1) %>% mutate(age = 3:9) %>% 
      gather(Q, value, -age) %>% 
      mutate(model = "Exogenous disease risk")
  ) %>% 
  filter(Q %in% c("0.3","0.8")) %>% 
  mutate(agegr = factor(age, levels = paste0(1:9), labels = colnames(pbc_spread))) %>%
  mutate(model = factor(model)) %>%
  ggplot(aes(x = agegr, y= value, group = interaction(model, Q), lty = model, color = Q)) + 
  geom_line(size = 1) +
  xlab("Age group") + ylab("Dose fraction") + ylim(0,1) +
  guides(color=guide_legend(title="Q (supply constraint)"))+
  theme(legend.position = "top",legend.direction = "vertical")
fig_static

ggsave("figures/presentation-dose_sharing_static_top.png", fig_static+theme(text = element_text(size=14)), width = 5.5, height=5)

##Different allocation strategies - benefits----
vac_top_p_adult <- function(p, pop) {
  pop <- pop/sum(pop)
  w <- rev(pop)
  ptemp <- p
  vrev <- vector(length = 7)
  for(j in 1:7){
    if(ptemp > 0)
      vrev[j] <- min(ptemp, w[j])
    ptemp <- ptemp - w[j]
  }
  rev(vrev)/pop
}

df_allocation_obj <- rbind(
  as.data.frame(nlopt_s0)[-1,] %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Optimal allocation",mixing="Heterogeneous"),
  as.data.frame(nlopt_s1)[-1,] %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Optimal allocation",mixing="Homogeneous"),
  as.data.frame(nonepi_s1) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Exogenous disease risk",mixing="Homogeneous"),
  as.data.frame(nonepi_s1) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Exogenous disease risk",mixing="Heterogeneous"),
  data_frame(
    "1"= rep(1,7),
    "0.9"= rep(0.9,7),
    "0.8"= rep(0.8,7),
    "0.7"= rep(0.7,7),
    "0.6"= rep(0.6,7),
    "0.5"= rep(0.5,7),
    "0.4"= rep(0.4,7),
    "0.3"= rep(0.3,7),
    "0.2"= rep(0.2,7),
    "0.1"= rep(0.1,7)) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Equal distribution",mixing="Homogeneous"),
  data_frame(
    "1"= rep(1,7),
    "0.9"= rep(0.9,7),
    "0.8"= rep(0.8,7),
    "0.7"= rep(0.7,7),
    "0.6"= rep(0.6,7),
    "0.5"= rep(0.5,7),
    "0.4"= rep(0.4,7),
    "0.3"= rep(0.3,7),
    "0.2"= rep(0.2,7),
    "0.1"= rep(0.1,7)) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Equal distribution",mixing="Heterogeneous"),
  data_frame(
    "1"= vac_top_p_adult(1,hic_pop[3:9]),
    "0.9"= vac_top_p_adult(0.9,hic_pop[3:9]),
    "0.8"= vac_top_p_adult(0.8,hic_pop[3:9]),
    "0.7"= vac_top_p_adult(0.7,hic_pop[3:9]),
    "0.6"= vac_top_p_adult(0.6,hic_pop[3:9]),
    "0.5"= vac_top_p_adult(0.5,hic_pop[3:9]),
    "0.4"= vac_top_p_adult(0.4,hic_pop[3:9]),
    "0.3"= vac_top_p_adult(0.3,hic_pop[3:9]),
    "0.2"= vac_top_p_adult(0.2,hic_pop[3:9]),
    "0.1"= vac_top_p_adult(0.1,hic_pop[3:9])) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Full dose (age prioritization)",mixing="Homogeneous"),
  data_frame(
    "1"= vac_top_p_adult(1,hic_pop[3:9]),
    "0.9"= vac_top_p_adult(0.9,hic_pop[3:9]),
    "0.8"= vac_top_p_adult(0.8,hic_pop[3:9]),
    "0.7"= vac_top_p_adult(0.7,hic_pop[3:9]),
    "0.6"= vac_top_p_adult(0.6,hic_pop[3:9]),
    "0.5"= vac_top_p_adult(0.5,hic_pop[3:9]),
    "0.4"= vac_top_p_adult(0.4,hic_pop[3:9]),
    "0.3"= vac_top_p_adult(0.3,hic_pop[3:9]),
    "0.2"= vac_top_p_adult(0.2,hic_pop[3:9]),
    "0.1"= vac_top_p_adult(0.1,hic_pop[3:9])) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "Full dose (age prioritization)",mixing="Heterogeneous"),
  data_frame(
    "1"= rep(0,7),
    "0.9"= rep(0,7),
    "0.8"= rep(0,7),
    "0.7"= rep(0,7),
    "0.6"= rep(0,7),
    "0.5"= rep(0,7),
    "0.4"= rep(0,7),
    "0.3"= rep(0,7),
    "0.2"= rep(0,7),
    "0.1"= rep(0,7)) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "No vaccine",mixing="Homogeneous"),
  data_frame(
    "1"= rep(0,7),
    "0.9"= rep(0,7),
    "0.8"= rep(0,7),
    "0.7"= rep(0,7),
    "0.6"= rep(0,7),
    "0.5"= rep(0,7),
    "0.4"= rep(0,7),
    "0.3"= rep(0,7),
    "0.2"= rep(0,7),
    "0.1"= rep(0,7)) %>% mutate(age = 3:9) %>% 
    gather(Q, value, -age) %>%
    mutate(model = "No vaccine",mixing="Heterogeneous")
)
df_allocation_obj$Q <- as.numeric(df_allocation_obj$Q)
df_allocation_obj <- df_allocation_obj %>% arrange(model,Q,age)

source("optimisation-epi/objective-functions.R")

unroll_x <- function(x, sub=0)
  c(sub,sub,x[1],x[2],x[3],x[4],x[5],x[6],x[7])

df_allocation_obj["obj"] <- 0
for (q in unique(df_allocation_obj$Q)){
  df_allocation_obj[(df_allocation_obj$model=="Optimal allocation")&
                      (df_allocation_obj$mixing=="Heterogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Optimal allocation")&
                                            (df_allocation_obj$mixing=="Heterogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 0,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="No vaccine")&
                      (df_allocation_obj$mixing=="Heterogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="No vaccine")&
                                            (df_allocation_obj$mixing=="Heterogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 0,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="Optimal allocation")&
                      (df_allocation_obj$mixing=="Homogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Optimal allocation")&
                                            (df_allocation_obj$mixing=="Homogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 1,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="No vaccine")&
                      (df_allocation_obj$mixing=="Homogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="No vaccine")&
                                            (df_allocation_obj$mixing=="Homogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 1,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="Equal distribution")&
                      (df_allocation_obj$mixing=="Homogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Equal distribution")&
                                            (df_allocation_obj$mixing=="Homogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 1,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="Equal distribution")&
                      (df_allocation_obj$mixing=="Heterogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Equal distribution")&
                                            (df_allocation_obj$mixing=="Heterogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 0,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="Exogenous disease risk")&
                      (df_allocation_obj$mixing=="Homogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Exogenous disease risk")&
                                            (df_allocation_obj$mixing=="Homogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 1,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="Exogenous disease risk")&
                      (df_allocation_obj$mixing=="Heterogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Exogenous disease risk")&
                                            (df_allocation_obj$mixing=="Heterogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 1,
                        outcome = 'D',
                        ret = 1)
  df_allocation_obj[(df_allocation_obj$model=="Full dose (age prioritization)")&
                      (df_allocation_obj$mixing=="Homogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Full dose (age prioritization)")&
                                            (df_allocation_obj$mixing=="Homogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 1,
                        outcome = 'D',
                        ret = 1,
                        full=1)
  df_allocation_obj[(df_allocation_obj$model=="Full dose (age prioritization)")&
                      (df_allocation_obj$mixing=="Heterogeneous")&
                      (df_allocation_obj$Q==q),"obj"] <- model_fd_static(unroll_x(
                        df_allocation_obj[(df_allocation_obj$model=="Full dose (age prioritization)")&
                                            (df_allocation_obj$mixing=="Heterogeneous")&
                                            (df_allocation_obj$Q==q),"value"],sub=0),
                        homogen = 0,
                        outcome = 'D',
                        ret = 1,
                        full=1)
}

#Homogeneous mixing
df_allocation_obj %>% 
  select(-age,-value) %>% 
  filter(mixing=="Homogeneous") %>% 
  group_by(Q) %>%
  mutate(obj=(1-obj/max(obj))) %>%
  ungroup() %>%
  unique() %>% 
  ggplot(aes(x = Q, y = obj, color = model)) +
  geom_line(size=1) +
  guides(color=guide_legend(title="Allocation Strategy",ncol=1))+
  xlab("Q (supply constraint)") + ylab("Fraction of Deaths Averted")+ylim(0,1)+
  # scale_color_viridis_d(begin = 0.1,end = 0.8)+
  theme(legend.position = "bottom")#,legend.direction = "vertical"

df_hom <- df_allocation_obj %>% 
  select(-age,-value) %>% 
  filter(mixing=="Homogeneous") %>% 
  group_by(Q) %>%
  mutate(obj=(1-obj/max(obj))) %>%
  ungroup() %>% 
  filter(Q %in% c(0.3,0.8)) %>% 
  unique() %>% 
  arrange(Q)
write.csv(df_hom,'benefits_allocation_homogeneous.csv',row.names = F)
#Heterogeneous mixing
df_allocation_obj %>% 
  select(-age,-value) %>% 
  filter(mixing=="Heterogeneous") %>% 
  group_by(Q) %>%
  mutate(obj=(1-obj/max(obj))) %>%
  ungroup() %>%
  unique() %>% 
  ggplot(aes(x = Q, y = obj, color = model)) +
  geom_line(size=1) +
  guides(color=guide_legend(title="Allocation Strategy",ncol=1))+
  xlab("Q (supply constraint)") + ylab("Fraction of Deaths Averted")+ylim(0,1)+
  # scale_color_viridis_d(begin = 0.1,end = 0.8)+
  theme(legend.position = "bottom")#,legend.direction = "vertical"

df_het <- df_allocation_obj %>%
  select(-age,-value) %>% 
  filter(mixing=="Heterogeneous") %>% 
  group_by(Q) %>%
  mutate(obj=(1-obj/max(obj))) %>%
  ungroup() %>% 
  filter(Q %in% c(0.3,0.8)) %>% 
  unique() %>% 
  arrange(Q)
write.csv(df_het,'benefits_allocation_heterogeneous.csv',row.names = F)
  
##Slide Fractional dosing with uniform dose----
le2.presentation <- df_efficacy_delta_raw %>%
  filter(d1 %in% le_speeds) %>%
  mutate(delta1 = 1/d1) %>% 
  select(delta1, e, model, i,d,harm) %>%
  gather(var, value, -delta1, -e, -model) %>%
  group_by(model,var) %>%
  mutate(ref = value[e == .95 & delta1 == default_delta_value]) %>%
  filter(e %in% c(seq(.5, .9, .1),0.95)) %>%
  ungroup() %>%
  mutate(r = (value/ref)) %>% 
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Fractional dosing better by >=5%",
                                    
                                    "Comparable (+-5%)", 
                                    "Base case better by >=5%"
                         ))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(speedup = factor(round(delta1/default_delta_value, 1))) %>%
  # mutate(speedup = factor(as.percent(delta1, 2))) %>%
  mutate(e = factor(e)) %>%
  filter(var == "Deaths" & model == "Fast growth") %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Efficacy") + 
  xlab("Fold increase in capacity (vs 0.25% base case)") + 
  # scale_x_continuous(breaks = 1/d1_general,
  # labels = as.percent(1/d1_general))
  geom_text(aes(label = value), color = "white", size = 5)
le2.presentation

ggsave("figures/presentation-le_optimal_deaths_fast.png", le2.presentation+theme(text = element_text(size=15),
                                                                                 legend.text = element_text(size=15),
                                                                                 legend.direction = "vertical"), 
       width = 7, height=5.5)

# Alternative figure with reciprocals of x on the x axis (MK sugestions 25 Apr 2021):
le2.presentation <- df_efficacy_delta_raw %>%
  filter(d1 %in% le_speeds) %>%
  mutate(delta1 = 1/d1) %>% 
  select(delta1, e, model, i,d,harm) %>%
  gather(var, value, -delta1, -e, -model) %>%
  group_by(model,var) %>%
  mutate(ref = value[e == .95 & delta1 == default_delta_value]) %>%
  filter(e %in% c(seq(.5, .9, .1),0.95)) %>%
  ungroup() %>%
  mutate(r = (value/ref)) %>% 
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Fractional dosing better by >=5%",
                                    "Comparable (+-5%)", 
                                    "Base case better by >=5%"
                         ))) %>%
  mutate(var = factor(var, levels = c("i", "d", "harm"),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(speedup = factor(round(delta1/default_delta_value, 1),
                          levels = c(1, 1.2, 1.4, 1.6, 2, 4, 8),
                          labels = c("1", "5/6", "5/7", "5/8", "1/2", "1/4", "1/8"))) %>%
  # mutate(speedup = factor(as.percent(delta1, 2))) %>%
  mutate(e = factor(e)) %>%
  filter(var == "Deaths" & model == "Fast growth") %>%
  ggplot(aes(x = speedup, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Efficacy") + 
  xlab("Dose (1 is full dose)") + 
  # scale_x_continuous(breaks = 1/d1_general,
  # labels = as.percent(1/d1_general))
  geom_text(aes(label = value), color = "white", size = 5)
le2.presentation
ggsave("figures/presentation-le_optimal_deaths_fast.png", le2.presentation+
         theme(text = element_text(size=15),
         legend.text = element_text(size=15),
         legend.direction = "vertical"), 
       width = 7, height=5.5)


##Slide Importance of starting early----
death_rates <- data.frame(delay = (1:nrow(w))/30, imax = w[,"D",]/max(w[,"D",]), imin = 0)

fig_delay_presentation <- select(df_delay_vac, d1, delay, d) %>%
  group_by(delay) %>%
  mutate(d = 1 - d/max(d)) %>%
  ungroup() %>%
  filter(d1 %in% c(d1_general)) %>%
  mutate(d1 = factor(d1, levels = c(d1_general), labels = c(as.percent(1/d1_general)))) %>%
  ggplot() + 
  geom_ribbon(data = inf_rates, aes(x=delay, ymin=imin, ymax=imax), alpha = .10) +
  geom_line(aes(x=delay, y=d, group=d1, color=d1)) +
  theme(legend.position = "top") +
  scale_color_discrete(name = "Vaccinated per day") +
  xlab("Vaccination start date (months until infection peak)") + 
  ylab("Fraction of deaths averted (vs no vaccination)") +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6, 7, 8), labels = -(0:8) + 4, limits = c(0, 9))

fig_delay_presentation

ggsave("figures/presentation-delay_impact_deaths_fast.png", fig_delay_presentation+theme(text=element_text(size=13)), width = 6.5, height=4)

##Slide Using less effective vaccine earlier----
delay_optimal_presentation <- 
  gg_delay %>% 
  filter(delay %in% (30.5*(0:6))) %>%
  spread(scenario, value) %>%
  mutate(r = `No switching (V2 only)`/`No early vaccine (V1 only)`) %>%
  select(delay, e, model, var, r) %>%
  mutate(le_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                         labels = c("Immediately available vaccine better by >=5%",
                                    "Comparable (+-5%)", 
                                    "Delayed vaccine better by >=5%"
                         ))) %>%
  mutate(value = round(r, 2)) %>%
  mutate(delay = factor(delay/30.5)) %>%
  mutate(e = factor(e)) %>%
  filter(var == "Deaths" & model == "Fast growth") %>%
  ggplot(aes(x = delay, y = e, fill = le_better)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "top") +
  ylab("Efficacy") + 
  xlab("Months until delayed vaccine available") +
  geom_text(aes(label = value), color = "white", size = 5)  

delay_optimal_presentation

ggsave("figures/presentation-delay_optimal_deaths_fast.png", delay_optimal_presentation+theme(text = element_text(size=15),
                                                                                              legend.text = element_text(size=12),
                                                                                              legend.direction = "vertical"), 
       width = 7, height=5.5)

##Timing between doses for 2 dose vaccines----
fdf_presentation <- df_fdf %>% 
  filter(d1 > 40, d1 <= 1000) %>%
  filter(e >= .5) %>%
  select(d1, model, e, policy, d, harm, i) %>%
  gather(var, value, -d1, -model, -e, -policy) %>%
  group_by(d1, model, e, var) %>% 
  summarise(value_m = min(value[policy != "default"])/value[policy == "default"], 
            policy = policy[which.min(value)]
  ) %>%
  mutate(better = ifelse(value_m<=1&value_m>=0.95,1,NA)) %>% 
  # mutate(better = factor(better, levels=c(1), labels = c("Less than 5% better"))) %>% 
  mutate(policy = factor(policy, levels = c("fdf", "hybrid_8","hybrid_7",
                                            "hybrid_6","hybrid_5","hybrid_4","hybrid_3", "default"),
                         labels = c("FDF for all", "FDF under 80","FDF under 70",
                                    "FDF under 60","FDF under 50","FDF under 40",
                                    "FDF under 30", "SDF"))) %>%
  mutate(var = factor(var, levels = c("i", "d", ""),
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(delta1 = factor(1/d1,
                         levels = rev(1/fdf_speeds),
                         labels = as.percent(rev(1/fdf_speeds)))) %>%
  mutate(e = factor(e)) %>%
  mutate(value_m = round(value_m, 2)) %>%
  filter(var == "Deaths" & model == "Fast growth") %>%
  ggplot(aes(x = delta1, y = e, fill = policy)) + #, alpha=better)) + 
  geom_tile() +
  # scale_alpha_manual(values = c(0.7,1), name = "",labels=NULL,guide = 'none') +
  ylab("e1 (efficacy following 1st dose)") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),text=element_text(size=8)) +
  # scale_fill_viridis_d(option = "magma") +
  scale_fill_viridis_d(end = 0.7) +
  xlab(paste0(def_labels, " (1st dose, default policy)")) + 
  geom_text(aes(label = value_m), color = "white", alpha=1, size = 5)

fdf_presentation

ggsave("figures/presentation-fdf_deaths_fast.png", fdf_presentation+theme(text = element_text(size=15)), width = 7, height=5.5)

