visualize_general <- 
  function(result, 
           ) {
  as.data.frame(result)[-1, , drop = F] %>% 
    mutate(age = 3:9) %>% 
    gather(q, value, -age) %>% 
    mutate(mixing = "heterogeneous") %>% 
    mutate(delta1 = as.percent(round(1/as.numeric(s),4))) %>%
    mutate(agegr = factor(age, levels = paste0(1:9), labels = colnames(pbc_spread))) %>%
    ggplot(aes(x = agegr, y= value, group = delta1, color = delta1)) + 
    geom_line(size = 1) +
    xlab("Age group") + ylab("Dose fraction") + ylim(0,1) +
    guides(color=guide_legend(title="vaccination speed with\nfull dose (% per day)"))
  
  
  
}






fig_dynamic <-


ggsave("figures/dose_sharing_dynamic.pdf", fig_dynamic+theme(text = element_text(size=10)), width = 6.5, height=5)



result %>% as.data.frame() %>% 
  








