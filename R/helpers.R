

# Running the model function -----

# Vaccinate top p of the population, starting from the oldest
# Returns a vector of what proportion of each age group was vaccinated
vac_top_p <- function(p, pop) {
  pop <- pop/sum(pop)
  w <- rev(pop)
  ptemp <- p
  vrev <- vector(length = 9)
  for(j in 1:9){
    if(ptemp > 0)
      vrev[j] <- min(ptemp, w[j])
    ptemp <- ptemp - w[j]
  }
  rev(vrev)/pop
}

# Grabbing parameter sets (fast, slow, late, constant risk) -----

grab_2v_parms <- function(model){
  # ugly! but avoids some Env problems
  if(model == "pars_linear")  pars <- pars_le_linear
  if(model == "pars_le_cr")   pars <- pars_le_cr
  if(model == "pars_le_slow") pars <- pars_le_slow
  if(model == "pars_le_fast") pars <- pars_le_fast
  if(model == "pars_le_late") pars <- pars_le_late
  pars
}
grab_2d_parms <- function(model){
  if(model == "pars_linear")    pars <- pars_fdf_linear
  if(model == "pars_le_cr")     pars <- pars_fdf_cr
  if(model == "pars_le_slow")   pars <- pars_fdf_slow
  if(model == "pars_le_fast")   pars <- pars_fdf_fast
  if(model == "pars_le_late")   pars <- pars_fdf_late
  pars
}


# Single run = sr() ----
sr <- function(pars, f = "2v_v2") { 
  gnames <- 1:pars$Ngroups
  times <- 1:(pars$Ndays)
  pars$Ndays <- NULL
  
  pars_to_expand <- c("e1", "e2")
  for(p in pars_to_expand){
    if(length(pars[[p]]) == 1)
      pars[[p]] <- rep(pars[[p]], pars$Ngroups)
  }
    
  
  if(pars$Ngroups == 9)
    gnames <- colnames(pbc_spread)
  
  if(f == "2d") {
    mod <- odin_ode_2dose(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV", "cumI")
  }
  if(f == "2d_v2") {
    mod <- odin_ode_2dose_v2(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  }
  if(f == "2d_v3") {
    mod <- odin_ode_2dose_v3(user = pars)
    cnames <- c("S", "E0", "E1", "E2", "I0", "I1", "I2", 
                "R", "RV", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  }
  if(f == "2v") {
    mod <- odin_ode_2vaccines(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumI")
  }
  if(f == "2v_v2") {
    mod <- odin_ode_2vaccines_v2(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumV", "cumI")
  }
  
  y <- array(mod$run(times)[,-1], 
             dim = c(max(times), pars$Ngroups, pars$Nc), 
             dimnames = list(times, gnames, cnames)) %>%
    aperm(c(1, 3, 2))
}


# Generate a starting population -----
# with varying levels of initial
# infection, pre-existing immunity, appropriate dimensions
y0_gen <- function(Nc, Ngroups, 
                   pre_immunity = rep(0, Ngroups),
                   # v = rep(0, Ngroups), 
                   ii = 5e-03, 
                   S = 1, E = 2, I = 3, R = 4){
  y0_default <- matrix(0, Nc, Ngroups)
  y0_default[S,] <- 1-pre_immunity
  y0_default[E,] <- (y0_default[1,]*ii)/2
  y0_default[I,] <- (y0_default[1,]*ii)/2
  # subtract initial infections
  y0_default[S,] <- y0_default[1,] - y0_default[2,]
  y0_default[R,] <- pre_immunity
  y0_default
}


# Gathering outputs (derived metrics of burden) -----

# Infection burden (only if no loss of immunity):
bi <- function(sr, pop, sub = 0) { 
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  y[dim(y)[1], "R", 1]
}
# Death burden:
bd <- function(sr, pop) {
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  y[dim(y)[1], "D", 1]
}
b_any <- function(sr, pop, comp, t = dim(sr)[1]) {
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  y[t, comp, 1]
}

# Time to herd immunity
tthi <- function(sr, pop, th=.5, r=FALSE) {
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  im <- y[, "P1", 1] + y[, "P2", 1]
  if(r)
    im <- im + y[, "R", 1]
  x <- im > th
  if(any(x))
    return(min(which(x)))
  else
    return(Inf)
}

# Benefit
benefit_p <- function(p, x1 = 0, x2=.7, y1=0) {
  if(p < x1)
    return(y1*(p/x1))
  if(p < x2)
    return((1-y1)*((p-x1)/(x2-x1)) + y1)
  return(1)
}
benefit <- function(sr, pop, r = FALSE, benefit_f = benefit_p) {
  y <- sr %>% rescale_rcs(merge = T, pop_sizes = pop/sum(pop))
  im <- y[, "P1", 1] + y[, "P2", 1]
  if(r)
    im <- im + y[, "R", 1]
  d <- dim(y)[1]
  sum(sapply(im, benefit_f))/d
}

main_metrics <- function(y, pop, vat = 31) {
  v1 <- b_any(y, pop, "cumV", vat)
  c("i" = b_any(y, pop, "cumI"), 
    "d" = b_any(y, pop, "D"), 
    "v1" = v1, 
    "tt50" = tthi(y, pop),
    # Benefit integral, over vaccinations only:
    # "harm_v"  = 1-benefit(y, pop),
    "harm_vr" = harm(y))
}
metric_nms <- c("i", "d", "v1", "tt50", "harm")


check0sums <- function(ode, maxC=14) {
  apply(x, c(1,2), \(x) sum(x)/9) %>% apply(1, \(x) sum(x[1:maxC]))
}
