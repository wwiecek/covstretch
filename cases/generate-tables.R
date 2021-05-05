library(xtable)
library(kableExtra)

table_g3.kappa <- df_kappa %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) %>%
  filter(d1 %in% c(d1_general, Inf)) %>%
  mutate(ref_i = i[d1 > 1460]) %>%
  mutate(ref_d = d[d1 > 1460]) %>%
  ungroup() %>%
  mutate(ri = 1-(i/ref_i)) %>%
  mutate(rd = 1-(d/ref_d)) %>%
  group_by(d1, model, e, kappa) %>%
  filter(e == .95) %>%
  mutate(d1 = as.percent(1/d1)) %>%
  select(d1, model, i, ri, d, rd) %>%
  mutate(d = round(d*1e05)) %>%
  gather(variable, value, -d1, -model, -e, -kappa) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "ri", "rd"),
                           labels = c("Infections", "Deaths per 100,000", "Fraction of infections averted", "Fraction of deaths averted"))) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)


table_g3 <- df_efficacy_delta %>% 
  filter(e == .95) %>%
  filter(round(d1,4) %in% c(round(d1_general,4), Inf)) %>% 
  mutate(d1 = as.percent(1/d1)) %>%
  select(d1, model, i, ri, d, rd) %>%
  mutate(d = round(d*1e05)) %>%
  # mutate(i = i*1e05) %>%
  gather(variable, value, -d1, -model, -e) %>%
  mutate(variable = factor(variable, 
                           levels = c("i", "d", "ri", "rd"),
                           labels = c("Fraction infected", "Deaths per 100,000", "Fraction of infections averted", "Fraction of deaths averted"))) %>%
  mutate(value = round(value, 2)) %>%
  spread(d1, value) %>% 
  arrange(variable, model) %>% ungroup() %>%
  select(-e)

# table_g3[table_g3$variable=="d",] <- mapply(  # Iterate over each column
#   function(df) {
#     formatC(df, format="f", digits=0)
#   },
#   df=table_g3[table_g3$variable=="d",]
# )
  
#formatC(table_g3[table_g3$variable=="d",], format="f", digits=0)


print ("Generating tables G3 and the analagous for kappa")
print (xtable(table_g3))

kbl(table_g3,"latex",align = "r", vline = "", row.names = FALSE) %>%
  kable_styling(full_width = F) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = F) %>%
  collapse_rows(columns = 2, valign = "top")

print (xtable(table_g3.kappa,digits=c(0,5,0,0,2,2,2,2,2,2)))

kbl(table_g3.kappa,"latex",align = "r", vline = "", row.names = FALSE) %>%
  kable_styling(full_width = F) %>%
  kable_paper(full_width = F) %>%
  column_spec(1, bold = F) %>%
  collapse_rows(columns = 3, valign = "top")

# This has to be ran after vrf.R
print (xtable(tab_screening))
