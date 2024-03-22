source("optimisation-epi/nlopt_general.R")

time_elapsed_vec <- c()
results_list <- list()

Q <- 365

for (n_iter in c(1,10,100,1000,10000)) {
  a <- Sys.time()
  results <- opt_general(q = c(0.5), initial_value = 0.2, iterations = n_iter, homogen_mixing = FALSE, static = FALSE, objective = "cumI")[2:8,]  
  b <- Sys.time()
  c <- b - a
  print(c)
  time_elapsed_vec[length(time_elapsed_vec) + 1] <- difftime(time1 = b, time2 = a, units = "secs")
  results_list[[length(results_list)+1]] <- results
}
