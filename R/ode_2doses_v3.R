odin_ode_2dose_v3 <- odin::odin({
  
  initial(S[])      <- y0[1,i]
  initial(E0[])     <- y0[2,i]
  initial(E1[])     <- y0[3,i]
  initial(E2[])     <- y0[4,i]
  initial(I0[])     <- y0[5,i]
  initial(I1[])     <- y0[6,i]
  initial(I2[])     <- y0[7,i]
  initial(R[])      <- y0[8,i]
  initial(RV[])     <- y0[9,i]
  initial(D[])      <- y0[10,i]
  initial(P1[])     <- y0[11,i]
  initial(N1[])     <- y0[12,i]
  initial(P2[])     <- y0[13,i]
  initial(N2[])     <- y0[14,i]
  initial(cumV1[])  <- y0[15,i]
  initial(cumV2[])  <- y0[16,i]
  initial(cumV[])   <- y0[17,i]
  initial(cumI[])   <- y0[18,i]
  
  beta_matrix[,] <- contacts[i,j]*q[j]*(I1[j]+I2[j]+I0[j])
  beta[] <- sum(beta_matrix[i,])
  
  va1[] <- 1/((1 + exp(ta[i] - t))*(1 + exp(50*(cumV1[i] - vstop[i]))))
  
  # ODE equations are here:
  deriv(S[])       <- -(beta[i] + constantrisk)*S[i] + phi[i]*R[i] - va1[i]*delta1[i]*S[i]/(S[i]+R[i]*vrf)
  deriv(E0[])      <-  (beta[i] + constantrisk)*S[i]  - gamma1[i]*E0[i]
  deriv(E1[])      <-  (beta[i] + constantrisk)*N1[i] - gamma1[i]*E1[i]
  deriv(E2[])      <-  (beta[i] + constantrisk)*N2[i] - gamma1[i]*E2[i]
  deriv(I0[])      <-  gamma1[i]*E0[i] - gamma2[i]*I0[i]
  deriv(I1[])      <-  gamma1[i]*E1[i] - gamma2[i]*I1[i]
  deriv(I2[])      <-  gamma1[i]*E2[i] - gamma2[i]*I2[i]
  deriv(R[])       <-  gamma2[i]*((1-pd0[i])*I0[i]) - 
    phi[i]*R[i] - va1[i]*delta1[i]*R[i]*vrf/(S[i]+R[i])
  deriv(RV[])       <- gamma2[i]*((1-pd1[i])*I1[i] + (1-pd2[i])*I2[i]) - 
    phi[i]*RV[i]
  deriv(D[])       <-  gamma2[i]*(pd0[i]*I0[i] + pd1[i]*I1[i] + pd2[i]*I2[i])
  
  deriv(P1[])      <-  va1[i]*delta1[i]*(e1[i]*S[i]+R[i]*vrf)/(S[i]+R[i]*vrf) - delta2[i]*P1[i] - kappa1[i]*P1[i]
  deriv(N1[])      <-  va1[i]*delta1[i]*(1-e1[i])*S[i]/(S[i]+R[i]*vrf) - delta2[i]*N1[i] - (beta[i] + constantrisk)*N1[i] + kappa1[i]*P1[i]
  deriv(P2[])      <-  delta2[i]*(P1[i]+e2[i]*N1[i]) - kappa2[i]*P2[i]
  deriv(N2[])      <-  delta2[i]*((1-e2[i])*N1[i]) + kappa2[i]*P2[i] - (beta[i] + constantrisk)*N2[i] + phi[i]*RV[i]
  
  deriv(cumV1[])   <-  va1[i]*delta1[i]
  deriv(cumV2[])   <-  (delta2[i]*(P1[i]+N1[i]))
  deriv(cumV[])    <-  va1[i]*delta1[i] + (delta2[i]*(P1[i]+N1[i]))
  deriv(cumI[])    <-  gamma1[i]*(E0[i]+E1[i]+E2[i])
  
  # PARAMETERS: general
  Ngroups          <- user()
  Nc               <- user()
  constantrisk     <- user()
  vrf              <- user()
  y0[,]            <- user()
  contacts[,]      <- user()
  vstop[]          <- user()
  ta[]             <- user()
  q[]              <- user()
  
  # Parameters: ODE system
  gamma1[]       <- user()
  gamma2[]       <- user()
  delta1[]       <- user()
  delta2[]       <- user()
  kappa1[]       <- user()
  kappa2[]       <- user()
  e1[]           <- user()
  e2[]          <- user()
  phi[]          <- user()
  pd0[]          <- user()
  pd1[]          <- user()
  pd2[]          <- user()
  
  # Parameters: ODE system
  dim(e1)           <- Ngroups
  dim(e2)          <- Ngroups
  dim(gamma1)       <- Ngroups
  dim(gamma2)       <- Ngroups
  dim(delta1)       <- Ngroups
  dim(delta2)       <- Ngroups
  dim(kappa1)       <- Ngroups
  dim(kappa2)       <- Ngroups
  dim(phi)          <- Ngroups
  dim(pd0)          <- Ngroups
  dim(pd1)          <- Ngroups
  dim(pd2)          <- Ngroups
  
  dim(vstop)       <- Ngroups
  dim(ta)          <- Ngroups
  dim(va1)         <- Ngroups
  dim(q)           <- Ngroups
  dim(beta)        <- Ngroups
  dim(beta_matrix) <- c(Ngroups, Ngroups)
  dim(y0)          <- c(Nc, Ngroups)
  dim(contacts)    <- c(Ngroups, Ngroups)
  
  dim(S)     <- Ngroups
  dim(E0)    <- Ngroups
  dim(E1)    <- Ngroups
  dim(E2)    <- Ngroups
  dim(I0)    <- Ngroups
  dim(I1)    <- Ngroups
  dim(I2)    <- Ngroups
  dim(R)     <- Ngroups
  dim(RV)    <- Ngroups
  dim(N1)    <- Ngroups
  dim(N2)    <- Ngroups
  dim(P1)    <- Ngroups
  dim(P2)    <- Ngroups
  dim(D)     <- Ngroups
  dim(cumV)  <- Ngroups
  dim(cumV1) <- Ngroups
  dim(cumV2) <- Ngroups
  dim(cumI)  <- Ngroups
})
