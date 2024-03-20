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

load("results/results_full_dose.Rdata")
load("results/results_frac_uni.Rdata")
load("results/results_frac_age.Rdata")

# General plotting function for age-group-specific dose fractions. 
# Specify four options to be fixed and two options to be varied. 
plot_frac_age <- 
  function(fix_1_name, fix_1_value, 
           fix_2_name, fix_2_value,
           fix_3_name, fix_3_value,
           fix_4_name, fix_4_value,
           vary_1_name, vary_2_name) {
    results_frac_age %>% 
      filter(!!sym(fix_1_name) == fix_1_value, 
             !!sym(fix_2_name) == fix_2_value, 
             !!sym(fix_3_name) == fix_3_value, 
             !!sym(fix_4_name) == fix_4_value) %>% 
      mutate(q = as.factor(q)) %>% 
      unnest(result_solution) %>% 
      ggplot(aes(x = age_group, y = solution)) + 
      geom_point(aes(color = !!sym(vary_1_name))) + 
      geom_line(aes(group = interaction(!!sym(vary_1_name), !!sym(vary_2_name)), 
                    color = !!sym(vary_1_name),
                    linetype = !!sym(vary_2_name))) +
      theme_bw() + 
      labs(x = "Age Group", y = "Optimal Dose Fraction")
  }

plot_frac_age("static", F, 
              "pdeath", "ifr_hic",
              "homogen_mixing", F, 
              "scenario", "pars_le_slow", 
              "q", "objective")

plot_frac_age("static", F, 
              "objective", "D",
              "pdeath", "ifr_hic", 
              "scenario", "pars_le_slow", 
              "q", "homogen_mixing")

plot_frac_age("static", F, 
              "objective", "D",
              "homogen_mixing", F, 
              "scenario", "pars_le_slow", 
              "q", "pdeath")

# ------------------------------------------------------------------------------

















