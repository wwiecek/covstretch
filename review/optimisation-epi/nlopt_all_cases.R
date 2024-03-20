# ------------------------------------------------------------------------------
# Generate the aggregate list of optimization results with different options. 
# 
# Three general categories are: 
# 1. Full dose
# 2. Optimal universal fractional dose
# 3. Optimal fractional dose for each age group
# 
# Output: three lists of results for the three categories. 
# The first 10 columns are different combinations of options, the later 4 are:
#   - result_solution: 
#       - 
#       - 
#       - 7 x 2 tibble of age group and optimal dose fraction
#   - result_objective: objective value evaluated at the solution
#   - result_n_iterations: the number of iterations in the solving process
#   - result_message: message from nloptr about the status of the optimization
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Full dose
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# 2. Optimal universal fractional dose
# ------------------------------------------------------------------------------





# ------------------------------------------------------------------------------
# 3. Optimal fractional dose for each age group
# ------------------------------------------------------------------------------

results_frac_age_dynamic <- 
  expand_grid(q = c(1000, 700, 300, 100),
              initial_value = 0.01, 
              objective = c("D", "cumI"), 
              dose_response = "covid_default",
              homogen_mixing = c(T, F),
              static = F,
              pdeath = c("ifr_hic", "ifr_lic"), 
              scenario = c("pars_le_slow", "pars_le_fast"),
              recurring = T,
              iterations = 100
              ) %>% 
  mutate(result = pmap(., opt_general)) %>% 
  unnest_wider(result, names_sep = "_")

# save(results_frac_age_dynamic, file = "results/results_frac_age_dynamic.Rdata")

results_frac_age_static <- 
  expand_grid(q = c(0.1, 0.3, 0.7, 1),
              initial_value = 0.01, 
              objective = c("D", "cumI"), 
              dose_response = "covid_default",
              homogen_mixing = c(T, F),
              static = T,
              pdeath = c("ifr_hic", "ifr_lic"), 
              scenario = c("pars_le_slow", "pars_le_fast"),
              recurring = T,
              iterations = 100
              ) %>% 
  mutate(result = pmap(., opt_general)) %>% 
  mutate(solution = pluck(result, "solution"), 
         objective_value = pluck(result, "objective"), 
         n_iterations = pluck(result, "n_iterations"), 
         message = pluck(result, "message")) %>% 
  unnest_wider(result, names_sep = "_")

# save(results_frac_age_static, file = "results/results_frac_age_static.Rdata")

results_frac_age <- results_frac_age_dynamic %>% 
  rbind(results_frac_age_static)

save(results_frac_age, file = "results/results_frac_age.Rdata")
