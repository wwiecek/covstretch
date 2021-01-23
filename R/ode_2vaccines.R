odin_ode_2vaccines <- odin::odin({
  
  initial(S[])      <- y0[1,i]
  initial(E[])      <- y0[2,i]
  initial(I[])      <- y0[3,i]
  initial(R[])      <- y0[4,i]
  initial(D[])      <- y0[5,i]
  initial(V1[])     <- y0[6,i]
  initial(N1[])     <- y0[7,i]
  initial(V2[])     <- y0[8,i]
  initial(N2[])     <- y0[9,i]
  initial(cumV1[])  <- y0[10,i]
  initial(cumV2[])  <- y0[11,i]
  initial(cumI[])   <- y0[12,i]
  
  # kappa1; kappa2; phi; e1; e2; delta1; delta2; gamma1; gamma2; pdeath; contacts; q;
  # doses_x; doses_y1; doses_y2
  
  beta_matrix[,] <- contacts[i,j]*q[j]*I[j]
  beta[] <- sum(beta_matrix[i,])
  
  doses_available1 <- interpolate(doses_x, doses_y1, "linear")
  doses_available2 <- interpolate(doses_x, doses_y2, "linear")
  tot_treated_pop1[] <- (cumV1[i]*pop_size[i])
  tot_treated_pop2[] <- (cumV2[i]*pop_size[i])
  trt_avl1 <- 1*(sum(tot_treated_pop1) < doses_available1)
  trt_avl2 <- 1*(sum(tot_treated_pop2) < doses_available2)
  t_group1[] <- (t > ta1[i])*trt_avl1*(t < ts1[i])
  t_group2[] <- (t > ta2[i])*trt_avl2
  
  # ODE equations are here:
  deriv(S[])      <- -(beta[i] + constantrisk)*S[i] + phi[i]*R[i] - 
    t_group1[i]*delta1[i]*S[i] -
    t_group2[i]*delta2[i]*S[i] 
  deriv(E[])      <-  (beta[i] + constantrisk)*(S[i]+N1[i]+N2[i]) - gamma1[i]*E[i]
  deriv(I[])      <-  gamma1[i]*E[i] - gamma2[i]*I[i]
  deriv(R[])      <-  (1-pdeath[i])*gamma2[i]*I[i] - phi[i]*R[i] - 
    R[i]*v_recovered_flag*(t_group1[i]*delta1[i] + t_group2[i]*delta2[i])
  deriv(D[])      <-  pdeath[i]*gamma2[i]*I[i]
  deriv(V1[])     <-  t_group1[i]*delta1[i]*(e1*S[i] + R[i]*v_recovered_flag) - kappa1[i]*V1[i]
  deriv(N1[])     <-  t_group1[i]*(1-e1)*delta1[i]*(S[i]) + kappa1[i]*V1[i]  - (beta[i] + constantrisk)*N1[i]
  deriv(V2[])     <-  t_group2[i]*delta2[i]*(e2*S[i]+R[i]*v_recovered_flag) - kappa2[i]*V2[i]
  deriv(N2[])     <-  t_group2[i]*(1-e2)*delta2[i]*(S[i]) + kappa2[i]*V2[i]  - (beta[i] + constantrisk)*N2[i]
  deriv(cumV1[])  <-  t_group1[i]*delta1[i]*(S[i] + R[i]*v_recovered_flag)
  deriv(cumV2[])  <-  t_group2[i]*delta2[i]*(S[i] + R[i]*v_recovered_flag)
  deriv(cumI[])   <-  (beta[i] + constantrisk)*(S[i]+N1[i]+N2[i])
  
  # PARAMETERS: general
  Ndays          <- user()
  Ngroups        <- user()
  Nc             <- user()
  constantrisk   <- user()
  v_recovered_flag <- user()
  y0[,]          <- user()
  pop_size[]     <- user()
  doses_x[]      <- user()
  doses_y1[]     <- user()
  doses_y2[]     <- user()
  contacts[,]    <- user()
  q[]            <- user()
  ta1[]          <- user()
  ta2[]          <- user()
  ts1[]          <- user()
  
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
  
  dim(doses_x)           <- Ndays
  dim(doses_y1)          <- Ndays
  dim(doses_y2)          <- Ndays
  dim(pop_size)          <- Ngroups
  dim(t_group1)          <- Ngroups
  dim(t_group2)          <- Ngroups
  dim(ta1)               <- Ngroups
  dim(ts1)               <- Ngroups
  dim(ta2)               <- Ngroups
  dim(q)                 <- Ngroups
  dim(beta)              <- Ngroups
  dim(tot_treated_pop1)  <- Ngroups
  dim(tot_treated_pop2)  <- Ngroups
  dim(beta_matrix)       <- c(Ngroups, Ngroups)
  dim(y0) <- c(Nc, Ngroups)
  dim(contacts) <- c(Ngroups, Ngroups)
  
  dim(S) <- Ngroups
  dim(E) <- Ngroups
  dim(I) <- Ngroups
  dim(R) <- Ngroups
  dim(V1) <- Ngroups
  dim(N1) <- Ngroups
  dim(V2) <- Ngroups
  dim(N2) <- Ngroups
  dim(D) <- Ngroups
  dim(cumV1) <- Ngroups
  dim(cumV2) <- Ngroups
  dim(cumI)  <- Ngroups
})
