library(kableExtra)
library(tidyverse)
source("project-setup.R")

fig_folder <- "figures"

country_case <- list(list('hic',1),list('lic',1))#,list('hic',0.25),list('lic',0.25)

gg1.out.ceiling <- data.frame()
gg2.out.ceiling <- data.frame()
df_efficacy_delta_raw.out.ceiling <- data.frame()
df_fdf.out.ceiling <- data.frame()

for (d in country_case){
  country <- as.character(d[1])
  default_supply_ceiling <- as.numeric(d[2])
  if (country=='hic'){
    pop <- hic_pop/sum(hic_pop)
    default_pdeath <- ifr_hic
  } 
  if (country=='lic'){
    pop <- lic_pop/sum(lic_pop)
    default_pdeath <- ifr_lic
  }
  
  source("R/setup.R")
  
  source("cases/fdf-prep-delta.R")
  source("cases/fdf-results.R")
  
  source("cases/prep-results.R")
  source("cases/general-example.R")
  #source("cases/lower-efficacy.R")
  #source("cases/lower-efficacy-delay.R")
  
  #source("cases/kappa-impact.R")
  
  df_efficacy_delta_raw.ceiling <- df_efficacy_delta_raw
  
  df_efficacy_delta_raw.ceiling$ceiling <- default_supply_ceiling
  df_efficacy_delta_raw.ceiling$country <- country
  
  df_efficacy_delta_raw.out.ceiling <- rbind(df_efficacy_delta_raw.out.ceiling,df_efficacy_delta_raw.ceiling)
  
  ##Dataframes for figure g1 (evolution of infection/death/vaccination over time)
  gg1.ceiling <- gg1.df
  
  gg1.ceiling$ceiling <- default_supply_ceiling
  gg1.ceiling$country <- country
  
  gg1.out.ceiling <- rbind(gg1.out.ceiling,gg1.ceiling)
  
  gg2.ceiling <- gg2.df
  
  gg2.ceiling$ceiling <- default_supply_ceiling
  gg2.ceiling$country <- country
  
  gg2.out.ceiling <- rbind(gg2.out.ceiling,gg2.ceiling)
  
  ##Dataframes for the fdf analysis
  df_fdf.ceiling <- df_fdf
  
  df_fdf.ceiling$ceiling <- default_supply_ceiling
  df_fdf.ceiling$country <- country
  df_fdf.out.ceiling <- rbind(df_fdf.out.ceiling,df_fdf.ceiling)
}

save(df_efficacy_delta_raw.out.ceiling,df_fdf.out.ceiling,gg1.out.ceiling,gg2.out.ceiling,
     file = "results/results_lic_hic_analysis.Rdata")
# load("results/results_ceiling_lic_analysis.Rdata")

df_efficacy_delta.out.ceiling <- 
  df_efficacy_delta_raw.out.ceiling %>%
  ungroup() %>%
  # filter(d1 %in% c(90, 180, 360, 730, 1460, Inf)) %>%
  filter(round(d1,4) %in% c(round(d1_general,4), Inf)) %>%
  group_by(model, e,ceiling,country) %>% 
  # mutate(benefit_vr = 1 - harm_vr) %>%
  mutate(ref_e = harm[d1 > 1460]) %>%
  mutate(ref_i = i[d1 > 1460]) %>%
  mutate(ref_d = d[d1 > 1460]) %>%
  ungroup() %>% 
  group_by(d1, model, e,ceiling,country) %>% 
  mutate(re = 1-(harm/ref_e)) %>%
  mutate(ri = 1-(i/ref_i)) %>%
  mutate(diffi = (i-ref_i)) %>%
  mutate(rd = 1-(d/ref_d)) %>%
  mutate(diffd = (d-ref_d))

####Preparing merged dataframes####

g2b.out.ceiling <- df_efficacy_delta.out.ceiling  %>%
  filter(e == .95) %>%
  filter(round(d1,4) %in% round(d1_general,4)) %>%
  select(ri, rd, ceiling, country) %>%
  gather(key, value, -d1, -model, -e, -ceiling, -country) %>%
  mutate(key = factor(key, levels = c("ri", "rd", "re"), 
                      labels = c("Infections", "Deaths", "Economic harm"))) %>%
  mutate(delta1 = 1/d1)

g3.out.ceiling <- df_efficacy_delta.out.ceiling %>% 
  filter(e == .95) %>%
  filter(round(d1,4) %in% c(round(d1_general,4), Inf)) %>% 
  mutate(d1 = as.percent(1/d1)) %>%
  select(d1, model, i, ri, d, rd, diffi, diffd, ceiling, country) %>%
  mutate(d = round(d*1e05,2)) %>%
  mutate(diffd = round(diffd*1e05,2)) %>%
  # mutate(i = i*1e05) %>%
  gather(variable, value, -d1, -model, -e, -ceiling, -country) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "harm", "ri", "rd", "re", "diffd", "diffi"),
                           labels = c("Infections", "Deaths per 100,000", "Economic harm",
                                      "Fraction of infections averted", "Fraction of deaths averted", "RE", "Difference in deaths", "Difference in infections"))) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)

####Fraction of Harm averted (plots)####
g2b.out.ceiling <- g2b.out.ceiling %>% 
  mutate(ceiling = factor(ceiling,
                          levels = c(0.25, 1),
                          labels = c("25%", "100%"))) %>% 
  mutate(country = factor(country,
                          levels = c("hic", "lic"),
                          labels = c("HIC", "LIC")))

g2.ceiling.fig <- ggarrange(g2b.out.ceiling %>% filter(delta1!=0&country=="HIC") %>%
                              ggplot(aes(x = delta1, y = value, color = model, linetype = ceiling)) + 
                              geom_line(size=0.8) +
                              geom_point(pch = 21, size = 1.5, fill = "white") +
                              facet_wrap(~key, scales = "free", ncol = 3) +
                              scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
                              scale_color_discrete(name = "Scenario") +
                              scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
                              theme(axis.text.x = element_text(angle = 45), legend.position = "top",text = element_text(size=8.3)) +
                              ylim(0, 1)+
                              xlab(def_labels$speed) + 
                              ylab("Fraction of harm averted") + ggtitle("HIC"),
                            g2b.out.ceiling %>% filter(delta1!=0&country=="LIC") %>% 
                              ggplot(aes(x = delta1, y = value, color = model, linetype = ceiling)) + 
                              geom_line(size=0.8) +
                              geom_point(pch = 21, size = 1.5, fill = "white") +
                              facet_wrap(~key, scales = "free", ncol = 3) +
                              scale_x_continuous(breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
                              scale_color_discrete(name = "Scenario") +
                              scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
                              theme(axis.text.x = element_text(angle = 45), legend.position = "top",text = element_text(size=8.3)) +
                              xlab(def_labels$speed) + 
                              ylim(0, 1)+
                              ylab("Fraction of harm averted") + ggtitle("LIC"),
                            common.legend = TRUE, ncol = 1)

head(g2b.out.ceiling)

g2.lic_hic.fig <- rbind(g2b.out.ceiling %>% filter(delta1!=0&ceiling=="100%"),
                         data.frame(d1=rep(Inf,12),model=as.factor(rep(c("Slow decrease","Slow growth","Fast growth"),4)),
                                    e=rep(0.95,12),ceiling=rep("100%",12),country=as.factor(rep(c("HIC","LIC"),6)),key=c(rep("Infections",6),rep("Deaths",6)),
                                    value=rep(0,12),delta1=rep(0,12))) %>%
  ggplot(aes(x = delta1, y = value*100, color = model, linetype = country)) + 
  geom_line(size=0.5) +
  # geom_point(size = 1) +
  facet_wrap(~key, scales = "free", ncol = 3) +
  scale_x_continuous(breaks = c(0.00,1/d1_general), labels = c("",as.percent(1/d1_general))) +
  lightness(scale_color_brewer(name = "Epidemic scenario",palette = "YlOrRd",direction = 1,labels = c("Slow-decrease", "Slow-growth", "Fast-growth")),scalefac(0.95))+
  # scale_color_discrete(name = "Epidemic scenario",labels = c("Slow-decrease", "Slow-growth", "Fast-growth")) +
  scale_fill_discrete(name = "Parameters") +
  scale_linetype_manual(values=c("solid", "dashed"),name = "Supply constraint") +
  theme(axis.text.x = element_text(angle = 45), legend.position = "top", legend.title = element_text(size=8),legend.text = element_text(size=7),
        text = element_text(size=9)) +
  ylim(0, 100)+
  xlab("Percentage of population vaccinated daily") + 
  ylab("% burden averted")
g2.lic_hic.fig

ggsave(paste0(fig_folder, "/harm_lic_hic.pdf"), g2.lic_hic.fig, width=width, height=0.6*width)

####Final number of infections and death cases (table - G3 equivalent)####
g3.out.ceiling <- g3.out.ceiling %>% 
  mutate(country = factor(country,
                          levels = c("hic", "lic"),
                          labels = c("HIC", "LIC"))) %>% 
  rename(region=country)

table.ceiling <- g3.out.ceiling %>% filter(variable %in% c("Fraction of infections averted","Fraction of deaths averted")) %>% 
  arrange(region,desc(ceiling),variable,model) %>% 
  select(variable,model,`0%`,`0.1%`,`0.25%`,`0.5%`,`0.75%`,`1%`,`2%`)
# write.csv(table.ceiling,'results/ceiling.csv')
# print (xtable(table.ceiling))

fmt.table <- kbl(table.ceiling, "latex", align = "r") %>%
  kable_styling(full_width = F, font_size = 14) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1, valign = "top")%>%
  pack_rows("High Income Countries Demographics",1,12,label_row_css="text-align: center; font-size: medium")%>%#,latex_align="c"
  pack_rows("Low Income Countries Demographics",13,24,label_row_css="text-align: center; font-size: medium")%>%
  pack_rows("100% Supply",1,6,label_row_css="text-align: center; font-size: small")%>%
  pack_rows("25% Supply",7,12,label_row_css="text-align: center; font-size: small")%>%
  pack_rows("100% Supply",13,18,label_row_css="text-align: center; font-size: small")%>%
  pack_rows("25% Supply",19,24,label_row_css="text-align: center; font-size: small")

print(fmt.table)

fmt.table %>% save_kable('results/ceiling.png')

# kbl(table.ceiling, align = "c") %>%
#   kable_paper(full_width = F) %>%
#   column_spec(1, bold = T) %>%
#   collapse_rows(columns = 1:2, valign = "top") %>% 
#   save_kable('results/ceiling.png')

####Number of vaccinated, infections and death cases in time (plot)####
gg1.out.ceiling <- gg1.out.ceiling %>% 
  mutate(ceiling = factor(ceiling))
gg1.plot.ceiling <- gg1.out.ceiling %>% filter(country=="hic") %>% 
  ggplot(aes(x = time, y = value, color = var, linetype = ceiling)) + geom_line(size=0.5) + facet_wrap(.~scenario, scales = "free") +
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + ylab("fraction infected") +
  scale_color_discrete(name = "") +
  scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
  theme(legend.position = "top", legend.box="vertical",legend.margin=margin())


gg2.out.ceiling <- gg2.out.ceiling %>% 
  mutate(ceiling = factor(ceiling,
                          levels = c(0.25, 1),
                          labels = c("25%", "100%")))
gg2.plot.ceiling <- gg2.out.ceiling %>% filter(country=="hic") %>% 
  ggplot(aes(x = time, y = value, color = var, linetype = ceiling)) + geom_line(size=0.5) + 
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + 
  ylab("fraction currently vaccinated (I)") +
  scale_color_discrete(name = "") +
  scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
  theme(legend.position = "top", legend.box="vertical",legend.margin=margin())

g1_joint.ceiling <- ggarrange(gg2.plot.ceiling + ggtitle("Vaccinations")+ theme(legend.spacing.x = unit(0.15, 'in'),legend.box="vertical",legend.margin=margin(), text = element_text(size=7)), common.legend = TRUE, 
                              gg1.plot.ceiling + ggtitle("Infections")+ theme(text = element_text(size=7)), 
                              widths = c(1,2.5))

ggsave(paste0(fig_folder, "/g1_joint_ceiling.pdf"),g1_joint.ceiling, width = 5.55, height=2)

##Analysis lower efficacy x speed----
#LIC
le2_lic <- df_efficacy_delta_raw.out.ceiling %>%
  filter(ceiling==0.25) %>% 
  filter(round(d1,4) %in% round(le_speeds,4)) %>%
  mutate(delta1 = 1/d1) %>% 
  select(delta1, e, model, i,d,harm,country) %>%
  gather(var, value, -delta1, -e, -model,-country) %>%
  group_by(model,var,country) %>%
  mutate(ref = value[e == .95 & delta1 == default_delta_value]) %>%
  filter(e %in% seq(.5, .9, .1)) %>%
  ungroup() %>%
  mutate(ref=as.numeric(ref),value=as.numeric(value)) %>% 
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
  facet_grid(country+var~model) + ylab("e2 (efficacy for the less effective vaccine)") + 
  xlab("delta2/delta1 (speed-up factor vs 0.25% base case)") + 
  # scale_x_continuous(breaks = 1/d1_general,
  # labels = as.percent(1/d1_general))
  geom_text(aes(label = value), color = "white", size = 2.5)

ggsave(paste0(fig_folder, "/le_optimal_lic_hic_supply25.pdf"), le2_lic, width = 6.5, height = 7.7)

##FDF analysis----
fig2.all_k.lic <- df_fdf.out.ceiling %>% 
  filter(ceiling==0.25) %>% 
  filter(d1 > 40, d1 <= 1000) %>%
  filter(e >= .5) %>%
  select(d1, model, e, policy, d, harm, i, country) %>%
  gather(var, value, -d1, -model, -e, -policy, -country) %>%
  group_by(d1, model, e, var, country) %>% 
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
  filter(var != "Economic harm") %>%
  ggplot(aes(x = delta1, y = e, fill = policy)) + #, alpha=better)) + 
  geom_tile() +
  # scale_alpha_manual(values = c(0.7,1), name = "",labels=NULL,guide = 'none') +
  facet_grid(country+var~model) + 
  ylab("e1 (efficacy following 1st dose)") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),text=element_text(size=8)) +
  # scale_fill_viridis_d(option = "magma") +
  scale_fill_viridis_d(end = 0.7) +
  xlab(paste0(def_labels, " (1st dose, default policy)")) + 
  geom_text(aes(label = value_m), color = "white", alpha=1, size = 2)

ggsave(paste0(fig_folder, "/fdf_allk_lic_hic_supply25.pdf"), fig2.all_k.lic, width = 6.5, height = 7.7)

# save(df_efficacy_delta_raw.out.ceiling,df_fdf.out.ceiling,gg1.out.ceiling,gg2.out.ceiling,
     # file = "results/results_ceiling_analysis_pt1of4.Rdata")
