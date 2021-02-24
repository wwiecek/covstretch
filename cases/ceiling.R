library(kableExtra)

country_case <- list(list('hic',1),list('hic',0.25),list('lic',1),list('lic',0.25))

gg1.out.ceiling <- data.frame()
gg2.out.ceiling <- data.frame()
df_efficacy_delta.out.ceiling <- data.frame()
df_fdf.out.ceiling <- data.frame()

####Running all cases####
for (d in country_case){
  country <- as.character(d[1])
  default_supply_ceiling <- as.numeric(d[2])
  if (country=='hic'){
    pop <- hic_pop/sum(hic_pop)
  } 
  if (country=='lic'){
    pop <- lic_pop/sum(lic_pop)
  }
  #I thought redefining defaul_supply_ceiling would be enough but I get different results when calling prioritisation at every iteration
  source("R/prioritisation.R")
  
  ##Dataframe for table g3 and the plot for harm reduction (can also be used for the lower efficacy analysis)
  df_efficacy_delta_raw.ceiling <- expand_grid(d1 = default_speeds,
                                               e = seq(.5, .95, .05),
                                               model = c("pars_le_cr", "pars_le_slow", "pars_le_fast")) %>%
    mutate(data = pmap(list(model, d1, e), function(x,y, z) data.frame(value = model_i(x,y, z), 
                                                                       var = metric_nms))) %>%
    unnest(data) %>%
    spread(var, value) %>%
    group_by(model, e)  %>%
    mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                          labels = c("Constant risk", "Slow growth", "Fast growth"))) 
  
  df_efficacy_delta.ceiling <- df_efficacy_delta_raw.ceiling %>%
    filter(round(d1,4) %in% c(round(d1_general,4), Inf)) %>%
    mutate(ref_e = harm[d1 > 1460]) %>%
    mutate(ref_i = i[d1 > 1460]) %>%
    mutate(ref_d = d[d1 > 1460]) %>%
    ungroup() %>%
    mutate(re = 1-(harm/ref_e)) %>%
    mutate(ri = 1-(i/ref_i)) %>%
    mutate(diffi = (i-ref_i)) %>%
    mutate(rd = 1-(d/ref_d)) %>%
    mutate(diffd = (d-ref_d)) %>%
    group_by(d1, model, e)
  
  df_efficacy_delta.ceiling$ceiling <- default_supply_ceiling
  df_efficacy_delta.ceiling$country <- country
  
  df_efficacy_delta.out.ceiling <- rbind(df_efficacy_delta.out.ceiling,df_efficacy_delta.ceiling)
  
  ##Dataframes for figure g1 (evolution of infection/death/vaccination over time)
  gg1.ceiling <- lapply(mlist, function(pars) {
    ll <- list(
      "Vaccinate 0.3% per day" = apap_2v(pars, 100/.3),
      "Vaccinate 0.5% per day" = apap_2v(pars, 100/.5),
      "Vaccinate 1% per day"  = apap_2v(pars, 100/1),
      "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
      lapply(sr, f = "2v_v2") %>%
      lapply(rescale_rcs, pop, merge=T) %>% 
      abind::abind()
    as.data.frame(ll[,"I",]) %>%
      mutate(time = 1:nrow(.))
  }) %>%
    bind_rows(.id = "scenario") %>%
    gather(var, value, -time, -scenario)
  
  gg1.ceiling$ceiling <- default_supply_ceiling
  gg1.ceiling$country <- country
  
  gg1.out.ceiling <- rbind(gg1.out.ceiling,gg1.ceiling)
  
  pars <- mlist[[1]]
  ll <- list(
    "Vaccinate 0.3% per day" = apap_2v(pars, 100/.3),
    "Vaccinate 0.5% per day" = apap_2v(pars, 100/.5),
    "Vaccinate 1% per day"  = apap_2v(pars, 100/1),
    "No vaccination" = list_modify(pars, delta1 = rep(0, Ngroups))) %>%
    lapply(sr, f = "2v_v2") %>%
    lapply(rescale_rcs, pop, merge=T) %>% 
    abind::abind()
  gg2.ceiling <- as.data.frame(ll[,"cumV",]) %>%
    mutate(time = 1:nrow(.)) %>%
    gather(var, value, -time)
  
  gg2.ceiling$ceiling <- default_supply_ceiling
  gg2.ceiling$country <- country
  
  gg2.out.ceiling <- rbind(gg2.out.ceiling,gg2.ceiling)
  
  ##Dataframes for the fdf analysis
  df_fdf.ceiling <- expand.grid(d1 = c(fdf_speeds, Inf),
                                model = c("pars_le_cr", "pars_le_slow", "pars_le_fast"), 
                                e = c(.4, .5, .6, .7, .8, .9, .95), 
                                policy = c("default", "fdf", "hybrid")) %>%
    mutate(data = pmap(list(model, d1, e, policy), 
                       function(x,y,z,a) data.frame(value = model_fdf(x,y,z,a), 
                                                    var = metric_nms))) %>%
    unnest(data) %>%
    spread(var, value) %>%
    mutate(model = factor(model, levels = c("pars_le_cr", "pars_le_slow", "pars_le_fast"),
                          labels = c("Constant risk", "Slow growth", "Fast growth"))) %>%
    group_by(d1, model, e, policy)
  
  df_fdf.ceiling$ceiling <- default_supply_ceiling
  df_fdf.ceiling$country <- country
  df_fdf.out.ceiling <- rbind(df_fdf.out.ceiling,df_fdf.ceiling)
}

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
                                      "RI", "RD", "RE", "Difference in deaths", "Difference in infections"))) %>%
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
                              scale_x_continuous(trans = 'log10', breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
                              scale_color_discrete(name = "Scenario") +
                              scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
                              theme(axis.text.x = element_text(angle = 45, size = 9), legend.position = "top") +
                              xlab(def_labels$speed) + 
                              ylab("Fraction of harm averted") + ggtitle("HIC"),
                            g2b.out.ceiling %>% filter(delta1!=0&country=="LIC") %>% 
                              ggplot(aes(x = delta1, y = value, color = model, linetype = ceiling)) + 
                              geom_line(size=0.8) +
                              geom_point(pch = 21, size = 1.5, fill = "white") +
                              facet_wrap(~key, scales = "free", ncol = 3) +
                              scale_x_continuous(trans = 'log10', breaks = 1/d1_general, labels = as.percent(1/d1_general)) +
                              scale_color_discrete(name = "Scenario") +
                              scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
                              theme(axis.text.x = element_text(angle = 45, size = 9), legend.position = "top") +
                              xlab(def_labels$speed) + 
                              ylab("Fraction of harm averted") + ggtitle("LIC"),
                            common.legend = TRUE, ncol = 1)

####Final number of infections and death cases (table - G3 equivalent)####
g3.out.ceiling <- g3.out.ceiling %>% 
  mutate(country = factor(country,
                          levels = c("hic", "lic"),
                          labels = c("HIC", "LIC"))) %>% 
  rename(region=country)

table.ceiling <- g3.out.ceiling %>% filter(variable %in% c("Infections","Deaths per 100")) %>% 
  arrange(region,desc(ceiling),variable,model) %>% 
  select(variable,model,`0%`,`0.25%`,`0.5%`,`0.75%`,`1%`,`2%`)
# write.csv(table.ceiling,'results/ceiling.csv')
# print (xtable(table.ceiling))

fmt.table <- kbl(table.ceiling,"latex", align = "r") %>%
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
# fmt.table %>% save_kable('results/ceiling.png')

# kbl(table.ceiling, align = "c") %>%
#   kable_paper(full_width = F) %>%
#   column_spec(1, bold = T) %>%
#   collapse_rows(columns = 1:2, valign = "top") %>% 
#   save_kable('results/ceiling.png')

####Number of vaccinated, infections and death cases in time (plot)####
gg1.out.ceiling <- gg1.out.ceiling %>% 
  mutate(ceiling = factor(ceiling))
gg1.plot.ceiling <- gg1.out.ceiling %>% filter(country=="hic") %>% 
  ggplot(aes(x = time, y = value, color = var, linetype = ceiling)) + geom_line(size=0.8) + facet_wrap(.~scenario, scales = "free") +
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + ylab("fraction infected") +
  scale_color_discrete(name = "") +
  scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
  theme(legend.position = "top", legend.box="vertical",legend.margin=margin())


gg2.out.ceiling <- gg2.out.ceiling %>% 
  mutate(ceiling = factor(ceiling,
                          levels = c(0.25, 1),
                          labels = c("25%", "100%")))
gg2.plot.ceiling <- gg2.out.ceiling %>% filter(country=="hic") %>% 
  ggplot(aes(x = time, y = value, color = var, linetype = ceiling)) + geom_line(size=0.8) + 
  xlab("time [days]") + scale_x_continuous(breaks = seq(0, 360, 120)) + 
  ylab("fraction currently vaccinated (I)") +
  scale_color_discrete(name = "") +
  scale_linetype_manual(values=c("dotted", "solid"),name = "Supply constraint") +
  theme(legend.position = "top", legend.box="vertical",legend.margin=margin())

g1_joint.ceiling <- ggarrange(gg2.plot.ceiling + ggtitle("Vaccinations")+ theme(legend.spacing.x = unit(0.1, 'in'),legend.box="vertical",legend.margin=margin(), text = element_text(size=7)), common.legend = TRUE, 
                              gg1.plot.ceiling + ggtitle("Infections")+ theme(text = element_text(size=7)), 
                              widths = c(1,2.5))

ggsave(paste0(fig_folder, "/ceiling.pdf"), g2.ceiling.fig, width = 5.55, height=5.55)
ggsave(paste0(fig_folder, "/g1_joint_ceiling.pdf"),g1_joint.ceiling, width = 5.55, height=2)
