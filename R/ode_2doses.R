odin_ode_2dose <- odin::odin({
  
  initial(S[])      <- y0[1,i]
  initial(E[])      <- y0[2,i]
  initial(I[])      <- y0[3,i]
  initial(R[])      <- y0[4,i]
  initial(D[])      <- y0[5,i]
  initial(P1[])     <- y0[6,i]
  initial(N1[])     <- y0[7,i]
  initial(P2[])     <- y0[8,i]
  initial(N2[])     <- y0[9,i]
  initial(cumV[])   <- y0[10,i]
  initial(cumI[])   <- y0[11,i]
  
  beta_matrix[,] <- contacts[i,j]*q[j]*I[j]
  beta[] <- sum(beta_matrix[i,])
  
  doses_available <- interpolate(doses_x, doses_y, "linear")
  tot_treated_pop[] <- (cumV[i]*pop_size[i])
  trt_avl <- 1*(sum(tot_treated_pop) < doses_available)
  # trt_avl <- 1*(sum(tot_treated_pop) < .2*(t/30))
  t_group[] <- (t > ta[i])
    
  # ODE equations are here:
  deriv(S[])      <- -(beta[i] + constantrisk)*S[i] + phi[i]*R[i] - trt_avl*delta1[i]*S[i]*t_group[i]
  deriv(E[])      <-  (beta[i] + constantrisk)*(S[i]+N1[i]+N2[i]) - gamma1[i]*E[i]
  deriv(I[])      <-  gamma1[i]*E[i] - gamma2[i]*I[i]
  deriv(D[])      <-  pdeath[i]*gamma2[i]*I[i]
  deriv(R[])      <-  (1-pdeath[i])*gamma2[i]*I[i] - phi[i]*R[i] - trt_avl*delta1[i]*R[i]*t_group[i]*v_recovered_flag
  
  deriv(P1[])     <-  trt_avl*t_group[i]*delta1[i]*(e1*S[i]+R[i]*v_recovered_flag) - trt_avl*delta2[i]*P1[i] - kappa1[i]*P1[i]
  deriv(N1[])     <-  (1-e1)*trt_avl*t_group[i]*delta1[i]*(S[i]) - trt_avl*delta2[i]*N1[i] - (beta[i] + constantrisk)*N1[i] + kappa1[i]*N1[i]
  deriv(P2[])     <-  trt_avl*t_group[i]*e2*delta2[i]*(P1[i]+N1[i]) - kappa2[i]*P2[i]
  deriv(N2[])     <-  trt_avl*t_group[i]*(1-e2)*delta2[i]*(P1[i]+N1[i]) + kappa2[i]*N2[i] - (beta[i] + constantrisk)*N2[i] 
  
  deriv(cumV[])   <-  trt_avl*t_group[i]*(delta1[i]*(S[i] + R[i]*v_recovered_flag) + delta2[i]*(P1[i]+N1[i]))
  deriv(cumI[])   <-  gamma1[i]*E[i]
  
  # PARAMETERS: general
  Ndays            <- user()
  Ngroups          <- user()
  Nc               <- user()
  v_recovered_flag <- user()
  constantrisk     <- user()
  y0[,]            <- user()
  pop_size[]       <- user()
  doses_x[]        <- user()
  doses_y[]        <- user()
  contacts[,]      <- user()
  ta[]             <- user()
  q[]              <- user()
  
  # Parameters: ODE system
  gamma1[]       <- user()
  gamma2[]       <- user()
  delta1[]       <- user()
  delta2[]       <- user()
  kappa1[]       <- user()
  kappa2[]       <- user()
  e1             <- user()
  e2             <- user()
  phi[]          <- user()
  pdeath[]       <- user()

  # Parameters: ODE system
  dim(gamma1)       <- Ngroups
  dim(gamma2)       <- Ngroups
  dim(delta1)       <- Ngroups
  dim(delta2)       <- Ngroups
  dim(kappa1)       <- Ngroups
  dim(kappa2)       <- Ngroups
  dim(phi)          <- Ngroups
  dim(pdeath)       <- Ngroups
  
  dim(doses_x)     <- Ndays
  dim(doses_y)     <- Ndays
  dim(ta)          <- Ngroups
  dim(t_group)     <- Ngroups
  dim(pop_size)    <- Ngroups
  dim(q)           <- Ngroups
  dim(beta)        <- Ngroups
  dim(tot_treated_pop) <- Ngroups
  dim(beta_matrix) <- c(Ngroups, Ngroups)
  dim(y0)          <- c(Nc, Ngroups)
  dim(contacts)    <- c(Ngroups, Ngroups)
  
  dim(S)    <- Ngroups
  dim(E)    <- Ngroups
  dim(I)    <- Ngroups
  dim(R)    <- Ngroups
  dim(N1)   <- Ngroups
  dim(N2)   <- Ngroups
  dim(P1)   <- Ngroups
  dim(P2)   <- Ngroups
  dim(D)    <- Ngroups
  dim(cumV) <- Ngroups
  dim(cumI) <- Ngroups
})
