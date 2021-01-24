

# Running the model function -----

# Single run = sr()
sr <- function(pars, f = "2v") { 
  gnames <- 1:pars$Ngroups
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
  if(f == "2v") {
    mod <- odin_ode_2vaccines(user = pars)
    cnames <- c("S", "E", "I", "R", "D", "P1", "N1", "P2", "N2", "cumV1", "cumV2", "cumI")
  }
  times <- 1:(pars$Ndays)
  y <- array(mod$run(times)[,-1], 
             dim = c(pars$Ndays, pars$Ngroups, pars$Nc), 
             dimnames = list(times, gnames, cnames)) %>%
    aperm(c(1, 3, 2))
}


# Generate a starting population -----
# with varying levels of initial
# infection, pre-existing immunity, appropriate dimensions
y0_gen <- function(Nc, Ngroups, 
                   pre_immunity = rep(0, Ngroups), 
                   # v = rep(0, Ngroups), 
                   ii = 5e-03){
  y0_default <- matrix(0, Nc, Ngroups)
  y0_default[1,] <- 1-pre_immunity
  y0_default[2,] <- (y0_default[1,]*ii)/2
  y0_default[3,] <- (y0_default[1,]*ii)/2
  # subtract initial infections
  y0_default[1,] <- y0_default[1,] - y0_default[2,]
  y0_default[4,] <- pre_immunity
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
  min(which(im > th))
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


main_metrics <- function(y, pop) {
  if("cumV1" %in% dimnames(y)[[2]])
    v1 <- b_any(y, pop, "cumV1", 31)
  else
    v1 <- b_any(y, pop, "cumV", 31)
  c("i" = 100*b_any(y, pop, "cumI"), 
    "d" = 100*b_any(y, pop, "D"), 
    "v1" = v1, 
    "tt50" = tthi(y, pop),
    # Benefit integral, over vaccinations only:
    "harm_v"  = 1-benefit(y, pop),
    "harm_vr" = 1-benefit(y, pop, r = TRUE))
}
metric_nms <- c("i", "d", "v1", "tt50", "harm_v", "harm_vr")
