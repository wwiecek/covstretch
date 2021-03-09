odin_ode_2vaccines_v2 <- odin::odin({
  
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
  initial(cumV[])   <- y0[12,i]
  initial(cumI[])   <- y0[13,i]
  
  beta_matrix[,] <- contacts[i,j]*q[j]*I[j]
  beta[] <- sum(beta_matrix[i,])
  
  va1[] <- 1/((1 + exp(ta1[i] - t))*(1 + exp(1e03*cumV[i] - vstop[i]*1e03))*(1 + exp(t - ts1[i])))
  va2[] <- 1/((1 + exp(ta2[i] - t))*(1 + exp(1e03*cumV[i] - vstop[i]*1e03)))
  # Expansion factor:
  va1[] <- (1 + 1/(1+exp(1e03*(tmore1[i] - t))))*va1[i]
  va2[] <- (1 + 1/(1+exp(1e03*(tmore2[i] - t))))*va2[i]
  
  # ODE equations are here:
  deriv(S[])      <- -(beta[i] + constantrisk)*S[i] + phi[i]*R[i] - 
    (va1[i]*delta1[i] + va2[i]*delta2[i])*S[i]/(S[i]+R[i]*vrf) 
  deriv(E[])      <-  (beta[i] + constantrisk)*(S[i]+N1[i]+N2[i]) - gamma1[i]*E[i]
  deriv(I[])      <-  gamma1[i]*E[i] - gamma2[i]*I[i]
  deriv(R[])      <-  (1-pdeath[i])*gamma2[i]*I[i] - phi[i]*R[i] - 
    R[i]*vrf*(va1[i]*delta1[i] + va2[i]*delta2[i])/(S[i]+R[i])
  deriv(D[])      <-  pdeath[i]*gamma2[i]*I[i]
  deriv(V1[])     <-  va1[i]*delta1[i]*(e1*S[i] + R[i]*vrf)/(S[i]+R[i]*vrf) - kappa1[i]*V1[i]
  deriv(N1[])     <-  va1[i]*delta1[i]*(1-e1)*S[i]/(S[i]+R[i]*vrf) + kappa1[i]*V1[i]  - (beta[i] + constantrisk)*N1[i]
  deriv(V2[])     <-  va2[i]*delta2[i]*(e2*S[i]+R[i]*vrf)/(S[i]+R[i]*vrf) - kappa2[i]*V2[i]
  deriv(N2[])     <-  va2[i]*delta2[i]*(1-e2)*S[i]/(S[i]+R[i]*vrf) + kappa2[i]*V2[i]  - (beta[i] + constantrisk)*N2[i]
  deriv(cumV1[])  <-  va1[i]*delta1[i]
  deriv(cumV2[])  <-  va2[i]*delta2[i]
  deriv(cumV[])   <-  va1[i]*delta1[i] + va2[i]*delta2[i]
  deriv(cumI[])   <-  gamma1[i]*E[i]
  
  # PARAMETERS: general
  Ngroups        <- user()
  Nc             <- user()
  constantrisk   <- user()
  vrf            <- user()
  y0[,]          <- user()
  contacts[,]    <- user()
  vstop[]        <- user()
  q[]            <- user()
  ta1[]          <- user()
  ta2[]          <- user()
  ts1[]          <- user()
  tmore1[]       <- user()
  tmore2[]       <- user()
  
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
  
  dim(vstop)        <- Ngroups
  dim(va1)          <- Ngroups
  dim(va2)          <- Ngroups
  dim(ta1)          <- Ngroups
  dim(ts1)          <- Ngroups
  dim(ta2)          <- Ngroups
  dim(tmore1)       <- Ngroups
  dim(tmore2)       <- Ngroups
  dim(q)            <- Ngroups
  dim(beta)         <- Ngroups
  dim(beta_matrix)  <- c(Ngroups, Ngroups)
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
  dim(cumV) <- Ngroups
  dim(cumI)  <- Ngroups
})
