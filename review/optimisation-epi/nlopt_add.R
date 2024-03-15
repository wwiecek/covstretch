# ------------------------------------------------------------------------------
# Define the function to add optimization result to a general data frame
# ------------------------------------------------------------------------------

add_opt <- 
  function(result_all,
           q,
           initial_value,
           objective = "D",
           dose_response = function(x) -25.31701*x^1.037524 + 1.037524*25.31701*x,
           homogen_mixing = F, 
           static = T,
           pdeath = ifr_hic,
           scenario = "pars_le_slow",
           recurring = T,
           iterations = 100) {
    
    opt_result <- opt_general(q,
                              initial_value,
                              objective = objective,
                              dose_response = dose_response, 
                              homogen_mixing = homogen_mixing, 
                              static = static,
                              pdeath = pdeath,
                              scenario = scenario,
                              recurring = recurring,
                              iterations = iterations)
    
    result_all %>% rbind(opt_result)
  }
