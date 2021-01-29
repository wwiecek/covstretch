# Construct risk functions for different age groups -----
epi <- sr(list_modify(pars_le_fast, delta1 = rep(0, Ngroups), 
                      y0 = y0_gen(13, Ngroups, c(0.5, 0.5, rep(0,7)), 1e-03))) 
risk_age <- epi[360, c("cumI", "D"),]
risk_age_i <- as.numeric(risk_age[1,])
risk_age_d <- as.numeric(risk_age[2,])
normalising_f_i <- sum(risk_age_i*pop)
normalising_f_d <- sum(risk_age_d*pop)

linear_harm <- function(atrisk) {
  frac <- 1 - (sum(pop*atrisk)/sum(pop))
  if(frac < .7)
    frac <- frac/.7
  else
    frac <- 1
  return(1 - frac)
}
harm <- function(y) {
  yast <- y[,"S",] + y[,"N1",] + y[,"N2",] + y[,"I",] + y[,"E",]
  sum(apply(yast, 1, linear_harm))/nrow(yast)
}

harm_at_t <- function(atrisk, lambda = 0.5) {
  lambda*sum(risk_age_i*atrisk*pop)/normalising_f_i +
    (1-lambda)*sum(risk_age_d*atrisk*pop)/normalising_f_d
}
harm <- function(y) {
  yast <- y[,"S",] + y[,"N1",] + y[,"N2",] + y[,"I",] + y[,"E",]
  sum(apply(yast, 1, harm_at_t, lambda = 1))/nrow(yast)
}

