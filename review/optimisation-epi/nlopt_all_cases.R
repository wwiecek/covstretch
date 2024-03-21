# ------------------------------------------------------------------------------
# Generate the aggregate list of optimization results with different options. 
# 
# Three general categories are: 
# 1. Full dose (no optimization needed)
# 2. Optimal universal fractional dose (1-dimensional optimization)
# 3. Optimal fractional dose for each age group (7-dimensional optimization)
# 
# Output: three lists of results for the three categories. 
# The first 10 columns are different combinations of options, the later 4 are:
#   - result_solution: 
#       1. N.A.
#       2. numeric value of the optimal universal dose fraction 
#       3. 7 x 2 tibble of age groups and corresponding optimal dose fractions
#   - result_objective: 
#       1 & 2 & 3: objective value evaluated at the solution
#   - result_n_iterations:
#       1: N.A.
#       2 & 3: the number of iterations in the solving process
#   - result_message: 
#       1: N.A.
#       2 & 3: message from nloptr about the status of the optimization
# ------------------------------------------------------------------------------


# First, generate all combinations of options. 
# Note that q is defined differently for dynamic and static problems, so we 
# run them separately and rbind later. 

cases_dynamic <- expand_grid(q = c(1000, 700, 300, 100),
                             initial_value = 0.01, 
                             objective = c("D", "cumI"), 
                             dose_response = "covid_default",
                             homogen_mixing = c(T, F),
                             static = F,
                             pdeath = c("ifr_hic", "ifr_lic"), 
                             scenario = c("pars_le_slow", "pars_le_fast"),
                             recurring = T,
                             iterations = 100)

cases_static <- expand_grid(q = c(0.1, 0.3, 0.7, 1),
                            initial_value = 0.01, 
                            objective = c("D", "cumI"), 
                            dose_response = "covid_default",
                            homogen_mixing = c(T, F),
                            static = T,
                            pdeath = c("ifr_hic", "ifr_lic"), 
                            scenario = c("pars_le_slow", "pars_le_fast"),
                            recurring = T,
                            iterations = 100)

# ------------------------------------------------------------------------------
# 1. Full dose
# ------------------------------------------------------------------------------

# Define wrapper function for the objective function with constant full dose. 
# The purpose of the wrapper is to fit the data frame structure defined above. 
obj_full <- function(q, initial_value, objective, dose_response, homogen_mixing,
                     static, pdeath, scenario, recurring, iterations) {
  dose_response <- ifelse(dose_response == "covid_default", 
                          function(x) -25.31701*x^1.037524 + 1.037524*25.31701*x, 
                          function(x) 0)
  
  default_pdeath <- ifelse(pdeath == "ifr_hic", ifr_hic, ifr_lic)
  
  if (static){
    model_fd_static(scenario = scenario, 
                    fd = c(0, 0, rep(1, 7)), 
                    phi_x = dose_response, 
                    ret = 1, 
                    objective = objective, 
                    homogen = homogen_mixing)
  } else {
    model_fd_dynamic(scenario = scenario, 
                     length_campaign = q, 
                     fd = rep(1, 9),  # try (00111111)
                     phi_x = dose_response, 
                     ret = 1, 
                     objective = objective, 
                     homogen = homogen_mixing)
  }
}

# Apply all combinations of options to the wrapper function.
results_full_dose_dynamic <- cases_dynamic %>% 
  mutate(result_objective = pmap(., obj_full))

results_full_dose_static <- cases_static %>% 
  mutate(result_objective = pmap(., obj_full))

results_full_dose <- results_full_dose_dynamic %>% rbind(results_full_dose_static)

save(results_full_dose, file = "results/results_full_dose.Rdata")

# ------------------------------------------------------------------------------
# 2. Optimal universal fractional dose
# ------------------------------------------------------------------------------

# Since we the dose fraction is universal, we redefine the dimension of the 
# problem to 1 and adjust the unroll_x function accordingly. 
n_x <- 1
unroll_x <- function(x, sub = 1) {c(sub, sub, rep(x, 7))} 

results_frac_uni_dynamic <- cases_dynamic %>% 
  mutate(result = pmap(., opt_general)) %>% 
  unnest_wider(result, names_sep = "_") %>% 
  mutate(result_solution = pull(result_solution, solution)[1])

results_frac_uni_static <- cases_static %>% 
  mutate(result = pmap(., opt_general)) %>% 
  unnest_wider(result, names_sep = "_") %>% 
  mutate(result_solution = pull(result_solution, solution)[1])

results_frac_uni <- results_frac_uni_dynamic %>% rbind(results_frac_uni_static)

save(results_frac_uni, file = "results/results_frac_uni.Rdata")

# ------------------------------------------------------------------------------
# 3. Optimal fractional dose for each age group
# ------------------------------------------------------------------------------

# Restore the original definition of n_x and unroll_x in nlopt_general.R. 
n_x <- 7
unroll_x <- function(x, sub = 1) {
  c(sub, sub, x[1], x[2], x[3], x[4], x[5], x[6], x[7])
} 

results_frac_age_dynamic <- cases_dynamic %>% 
  mutate(result = pmap(., opt_general)) %>% 
  unnest_wider(result, names_sep = "_")

# save(results_frac_age_dynamic, file = "results/results_frac_age_dynamic.Rdata")

results_frac_age_static <- cases_static %>% 
  mutate(result = pmap(., opt_general)) %>% 
  unnest_wider(result, names_sep = "_")

# save(results_frac_age_static, file = "results/results_frac_age_static.Rdata")

results_frac_age <- results_frac_age_dynamic %>% 
  rbind(results_frac_age_static)

save(results_frac_age, file = "results/results_frac_age.Rdata")
