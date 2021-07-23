# Define the dose-response function -----
phi_x <- function(x) 
  # 3.49706*sqrt(x) - 1.74853*x  -0.798528
  # sapply(3.49706*sqrt(x) - 1.74853*x  -0.798528, function(y) max(y,0))
  -25.31701*x^1.037524 + 1.037524*25.31701*x



# Define objective functions (static and dynamic cases) ------
model_fd_dynamic <- function(model, d1, fd, default_e1 = 0.95, 
                             rm = FALSE,
                             ret = 0,
                             outcome = "d",
                             homogen = FALSE) {
  e1 <- phi_x(fd)
  pars <- apap_2v(grab_2v_parms(model), fractional_dose = fd, len = d1)
  pars <- list_modify(pars, e1 = e1)
  if(homogen){
    # pars$contacts <- 1/Ngroups + 0*pars$contacts
    pars$contacts <- t(replicate(Ngroups, pop))
    pars$q <- ev*pars$q
  }
  y <- sr(pars, "2v_v2")
  if(rm) return(y)
  if(ret == 0)
    return(main_metrics(y, pop))
  if(ret == 1)
    y <- rescale_rcs(y, pop, TRUE)
    return(y[360,outcome,1])
}

model_fd_static <- function(v_prop, rm = FALSE, homogen = FALSE,
                                    ret = 0,
                                    outcome = "d",
                                    full=0) {
  e_vector <- phi_x(v_prop)
  if (full){
    e_vector <- v_prop*phi_x(1)
  }
  pars <- list_modify(
    pars_le_fast,
    y0 = y0_gen(13, Ngroups, pre_immunity = pre_immunity + (1-pre_immunity)*e_vector))
  if(homogen){
    # pars$contacts <- 1/Ngroups + 0*pars$contacts
    pars$contacts <- t(replicate(Ngroups, pop))
    pars$q <- ev*pars$q
  }
  # We do not update e1, because there is no vaccination past t=0 
  y <- sr(pars, "2v_v2")
  if(ret == 0)
    return(main_metrics(y, pop)[1:2])
  if(ret == 1)
    y <- rescale_rcs(y, pop, TRUE)
    return(y[360,outcome,1])
  
}