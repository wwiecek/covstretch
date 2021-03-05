#Plot - Increasing Relative Speed
df.vax_rate <- data.frame()
for (vax_rate in c(0.001,0.0025,0.005,0.0075,0.01,0.02)){
  df.vax_rate <- rbind(df.vax_rate, df_efficacy_delta_raw %>%
                         mutate(delta1 = round(1/d1,5)) %>% 
                         select(delta1, e, model, i,d,harm) %>%
                         gather(var, value, -delta1, -e, -model) %>%
                         mutate(ref_rate = vax_rate) %>%
                         group_by(model,var) %>%
                         mutate(ref = value[e == .95 & delta1 == vax_rate]) %>%
                         filter(delta1 %in% c(vax_rate*c(1.5,2,3,4)) & e == .95) %>%
                         ungroup() %>%
                         mutate(r = (value/ref)) %>% 
                         mutate(vax_rate_better = cut(r, c(-Inf, .95, 1.05, Inf), 
                                                      labels = c("Faster vaccination better by 5% or more",
                                                                 
                                                                 "Comparable (+-5%)", 
                                                                 "Slower vaccination better by at least 5%"
                                                      ))) %>%
                         mutate(var = factor(var, levels = c("i", "d", "harm"),
                                             labels = c("Infections", "Deaths", "Economic harm"))) %>%
                         mutate(value = round(r, 2)) %>%
                         mutate(speedup = factor(round(delta1/vax_rate, 1))) %>%
                         filter(var != "Economic harm"))
}

table.vax_rate <- kbl(df.vax_rate %>% select(model, speedup, var, ref_rate, r) %>% mutate(ref_rate=as.percent(ref_rate)) %>% 
                        pivot_wider(names_from = c("speedup"), values_from = r) %>% 
                        arrange(model, var, ref_rate) %>% select(var,ref_rate,`1.5`,`2`,`3`,`4`),
                      "latex", digits=3, align = "r", col.names = c('','Base rate','1.5','2','3','4')) %>%
  kable_styling(full_width = F, font_size = 14) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = T) %>%
  add_header_above(c(" " = 2,"Rate multiplier" = 4)) %>% 
  collapse_rows(columns = 1, valign = "top")%>%
  pack_rows("Constant risk",1,12,label_row_css="text-align: center; font-size: medium")%>%#,latex_align="c"
  pack_rows("Slow growth",13,24,label_row_css="text-align: center; font-size: medium")%>%
  pack_rows("Fast growth",25,36,label_row_css="text-align: center; font-size: medium")

print(table.vax_rate)

#table.vax_rate %>% save_kable("results/vax_rate_table.html")
table.vax_rate %>% save_kable("results/vax_rate_table.tex")


fig.vax_rate <- df.vax_rate %>% mutate(ref_rate = factor(ref_rate)) %>%
  ggplot(aes(x = speedup, y = ref_rate)) + geom_tile() +
  scale_fill_manual(values = c("grey60", "grey40", "grey20"), 
                    name = "") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(var~model) + ylab("baseline vaccination rate (delta1)") + 
  xlab("delta2/delta1 (speed-up factor vs base case)") + 
  scale_y_discrete(breaks = c(0.001,0.0025,0.005,0.0075,0.01,0.02),labels = as.percent(c(0.001,0.0025,0.005,0.0075,0.01,0.02))) +
  geom_text(aes(label = value), color = "white", size = 2)