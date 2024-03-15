# ------------------------------------------------------------------------------
# Add all cases to the (unnested) data frame of final results
# Can be easily modified to add particular cases. 
# ------------------------------------------------------------------------------

result_all <- 
  data.frame(age_group = character(), 
             solution = numeric(),
             q = numeric(),
             initial_value = numeric(),
             objective = character(),
             dose_response = character(), 
             homogen_mixing = logical(), 
             static = logical(), 
             pdeath = list(), 
             scenario = character(), 
             recurring = logical(), 
             objective_value = numeric(), 
             status = numeric(), 
             n_iterations = numeric(), 
             message = character())

for (q in c(1000, 400, 200, 100)){
  for (initial_value in c(0.5, 0.3, 0.1, 0.01)){
    for (objective in c("D", "cumI")){
      for (dose_response in list(function(x) -25.31701*x^1.037524 + 1.037524*25.31701*x)){
        for (homogen_mixing in c(T, F)){
          for (pdeath in list(ifr_hic, ifr_lic)){
            for (scenario in c("pars_le_slow", "pars_le_fast")){
              for (recurring in c(T, F)){
                result_all <- 
                  result_all %>% 
                  add_opt(q = q, 
                          initial_value = initial_value, 
                          objective = objective, 
                          dose_response = dose_response, 
                          homogen_mixing = homogen_mixing, 
                          static = static, 
                          pdeath = pdeath, 
                          scenario = scenario, 
                          recurring = recurring, 
                          iterations = iterations)
                save(result_all, file = "results/result_all_dynamic.RData")
                }
              }
            }
          }
        }
      }
    }
  }


















