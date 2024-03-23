# ------------------------------------------------------------------------------
# Visualize optimization results
# 
# # As a reminder, the options are:
# - objective = c("D", "cumI"), 
# - dose_response = "covid_default",
# - homogen_mixing = c(T, F),
# - static = c(T, F),
# - pdeath = c("ifr_hic", "ifr_lic"), 
# - scenario = c("pars_le_slow", "pars_le_fast"),
# - recurring = T,
# - iterations = 100)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. General plotting function for age-group-specific dose fractions. 
#    Specify four options to be fixed and two options to be varied. 
# ------------------------------------------------------------------------------

plot_frac_age <- 
  function(results_all,
           fix_1_name, fix_1_value, 
           fix_2_name, fix_2_value,
           fix_3_name, fix_3_value,
           fix_4_name, fix_4_value,
           vary_1_name, vary_2_name,
           show_frac_uni = F) {
    
    results_all_filter <- results_all %>% 
      filter(!!sym(fix_1_name) == fix_1_value, 
             !!sym(fix_2_name) == fix_2_value, 
             !!sym(fix_3_name) == fix_3_value, 
             !!sym(fix_4_name) == fix_4_value) %>% 
      mutate(q = as.factor(q))
    
    plot <- results_all_filter %>% 
      unnest(frac_age_sol) %>% 
      ggplot(aes(x = age_group, y = solution)) + 
      geom_point(aes(color = !!sym(vary_1_name))) + 
      geom_line(aes(group = interaction(!!sym(vary_1_name), !!sym(vary_2_name)), 
                    color = !!sym(vary_1_name),
                    linetype = !!sym(vary_2_name)))
    if (show_frac_uni) {
      plot <- plot + 
        geom_hline(data = results_all_filter,
                   aes(yintercept = frac_uni_sol,
                       group = interaction(!!sym(vary_1_name), !!sym(vary_2_name)), 
                       color = !!sym(vary_1_name),
                       linetype = !!sym(vary_2_name)))
    }
    plot + theme_bw() + 
      labs(x = "Age Group", y = "Optimal Dose Fraction") +
      ylim(0, 1) + 
      scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73", "#882E72"))
  }


# ------------------------------------------------------------------------------
# 2. Radar plot for comparing objective values from the three general strategies
# ------------------------------------------------------------------------------

# Subset the all-result data as desired and adjust axis label accordingly. 
# 
# As an example:
#
results_all_filter <- results_all %>%
  filter(static == T, objective == "D", scenario == "pars_le_slow") %>%
  mutate(setting = paste("q = ", q,
                         ", dose response = ", dose_response,
                         ", homogen mixing = ", homogen_mixing,
                         ", mortality risk = ", pdeath,
                         ", recurring = ", recurring,
                         sep = "")) %>% 
  select(setting, full_dose, frac_uni, frac_age) %>% 
  pivot_longer(cols = c(full_dose, frac_uni, frac_age), 
               names_to = "group", 
               values_to = "value") %>% 
  pivot_wider(names_from = setting, values_from = value)
  
radar_plot <- function(results_all_filter) {
    
  results_all_filter %>% 
    select(setting, full_dose, frac_uni, frac_age) %>% 
    pivot_longer(cols = c(full_dose, frac_uni, frac_age), 
                 names_to = "group", 
                 values_to = "value") %>% 
    pivot_wider(names_from = setting, values_from = value) %>% 
    ggradar(background.circle.colour = "white",
            values.radar = c("", "", ""),
            axis.label.size = 1,
            grid.label.size = 5,
            group.point.size = 1,
            legend.text.size = 5,
            group.line.width = 0.5,
            grid.line.width = 0.1,
            grid.min = 0, 
            grid.max = 0.005 / 8,
            grid.mid = 0.0025 / 8,
            group.colours = c("#E69F00", "#56B4E9", "#009E73")) + 
    theme_bw(base_size = 5) +
    theme(legend.position = "bottom")
}
  

















