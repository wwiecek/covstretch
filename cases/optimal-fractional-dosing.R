default_e1 <- 0.95
x <- seq(0.1, 1, length = 100)
y <- default_e1*x^0.25
plot(y~x, type = "l")
y <- default_e1*x^0.2
lines(y~x, type = "l")
y <- default_e1*x^0
lines(y~x, type = "l")



# K = youngest age group where fractional dosing applies
pars_for_fractional_d <- function(pars, d1, alpha, fd, K) {
  
}


model_fractional_d <- function(model, d1, alpha, fd, K, default_e1 = 0.95, rm = FALSE) {
  
  # Prepare parameters, taking fractional dosing into account:
  dose_response <- default_e1*fd^alpha
  N <- length(pop)
  
  if(K <= N){
    fd <- c(rep(fd, K), rep(1, N-K))
    e1 <- c(rep(dose_response, K), rep(default_e1, N-K))
  } else {
    fd <- rep(1, N)
    e1 <- rep(default_e1, N)
  }
  pars <- apap_2v(grab_2v_parms(model), fractional_dose = fd, len = d1)
  # Done.
  
  y <- sr(list_modify(pars, e1 = e1), "2v_v2")
  if(rm) return(y)
  main_metrics(y, pop)
}

df_fractional <- expand.grid(model = scenario_par_nms_2v,
                             d1 = c(d1_general, Inf),
                             alpha = c(0, 1, 0.25),
                             fractional_dose = c(0.1, 0.2, 0.5),
                             K = 3:9) %>%
  mutate(data = pmap(list(model, d1, alpha, fractional_dose, K), 
                     function(x,y,z,a,b) data.frame(value = model_fractional_d(x,y,z,a,b), 
                                                    var = metric_nms))) %>%
  unnest(data) %>%
  spread(var, value) %>%
  mutate(model = factor(model, levels = scenario_par_nms_2v,
                        labels = scenario_nms_2v)) 

