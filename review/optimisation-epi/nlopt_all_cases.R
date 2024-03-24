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
                             recurring = list(list(rep(0.8, 9)),
                                           list(rep(0.8, 9),
                                                c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5),
                                                c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5),
                                                c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5),
                                                c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5))),
                             iterations = 100)

cases_static <- expand_grid(q = c(0.1, 0.3, 0.7, 1),
                            initial_value = 0.01, 
                            objective = c("D", "cumI"), 
                            dose_response = "covid_default",
                            homogen_mixing = c(T, F),
                            static = T,
                            pdeath = c("ifr_hic", "ifr_lic"), 
                            scenario = c("pars_le_slow", "pars_le_fast"),
                            recurring = list(list(rep(0.8, 9)),
                                          list(rep(0.8, 9),
                                               c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5),
                                               c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5),
                                               c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5),
                                               c(0, 0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5))),
                            iterations = 100)

# ------------------------------------------------------------------------------
# 1. Full dose
# ------------------------------------------------------------------------------

# Define wrapper function for the objective function with constant full dose. 
# The purpose of the wrapper is to fit the data frame structure defined above. 
obj_full <- function(q, initial_value, objective, dose_response, homogen_mixing,
                     static, pdeath, scenario, recurring, iterations) {

  if (static){
      fd <- c(0, 0, rep(1, 7))
      phi_x <- dose_response
      vac_perc <- rep(0, 9)
      
      for (i in c(9:1)) {
        q <- q - pop[i]
        if (q > 0) {
          vac_perc[i] <- 1
        } else {
          vac_perc[i] <- q + pop[i]
          break
        }
      }
      
      if (phi_x == "covid_default") {phi_x <- function(x) -25.31701*x^1.037524 + 1.037524*25.31701*x
      } else if (phi_x == "flu_default") {phi_x <- function(x) 0}
      
      e_vector <- vac_perc * phi_x(fd)
      
      if (pdeath == "ifr_hic") {pdeath <- ifr_hic
      } else if (pdeath == "ifr_lic") {pdeath <- ifr_lic}
      
      pars <- list_modify(
        grab_2v_parms(scenario),
        y0 = y0_gen(13, Ngroups, pre_immunity = pre_immunity + (1-pre_immunity)*e_vector))
      if(homogen_mixing){
        pars$contacts <- t(replicate(Ngroups, pop))
        pars$q <- ev*pars$q
      }
      pars <- list_modify(pars, pdeath = pdeath)
      
      y <- sr(pars, "2v_v2")
      y <- rescale_rcs(y, pop, TRUE)
      return(y[360,objective,1])
    } else {
    model_fd_dynamic(scenario = scenario, 
                     length_campaign = q, 
                     fd = rep(1, 9),  # try (00111111)
                     phi_x = dose_response, 
                     pdeath = pdeath,
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
  mutate(result_solution = map(result_solution, pluck, 2)) %>% 
  mutate(result_solution = map(result_solution, pluck, 1))

results_frac_uni_static <- cases_static %>% 
  mutate(result = pmap(., opt_general)) %>% 
  unnest_wider(result, names_sep = "_") %>% 
  mutate(result_solution = map(result_solution, pluck, 2)) %>% 
  mutate(result_solution = map(result_solution, pluck, 1))

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

# ------------------------------------------------------------------------------
# 4. Merge 
# ------------------------------------------------------------------------------

load("results/results_frac_age.Rdata")
load("results/results_frac_uni.Rdata")

results_frac_uni <- results_frac_uni %>% 
  select(-result_message, -result_n_iterations)

results_frac_age <- results_frac_age %>% 
  select(-result_message, -result_n_iterations)

results_all <- 
  results_full_dose %>% 
  left_join(results_frac_uni, 
            by = c("q", "initial_value", "objective","dose_response", 
                   "homogen_mixing", "static", "pdeath", "scenario", 
                   "recurring", "iterations")) %>% 
  rename(full_dose = result_objective.x, 
         frac_uni = result_objective.y,
         frac_uni_sol = result_solution) %>% 
  left_join(results_frac_age,
            by = c("q", "initial_value", "objective","dose_response", 
                   "homogen_mixing", "static", "pdeath", "scenario", 
                   "recurring", "iterations")) %>% 
  rename(frac_age = result_objective,
         frac_age_sol = result_solution) %>% 
  mutate(across(c(full_dose, frac_uni, frac_age, frac_uni_sol), ~map_dbl(., 1)))

save(results_all, file = "results/results_all.Rdata")
