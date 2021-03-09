# Generate all results

# A few inputs needed to generate other parameters -----
# We use pbc_spread and default_cm from the default data inputs file
library(tidyverse)
load("data/default_inputs.Rdata")


# Sensitivity analyses (global parameters to modify) ------

# Main FDF assumptions
delay_default <- 28 - 10
delay_fdf <- 84 - 10
delay_hybrid <- c(rep(delay_fdf, 6), rep(delay_default, 3))

# Demographics (for comparing HIC vs LIC)
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
lic_pop <- pbc_spread[countries["Low-income countries"],] %>% as.numeric()
ifr_hic <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100
ifr_lic <- ifr_hic*(3.2/2)^(5:(-3))
default_pdeath <- ifr_hic

#concave function
A=.95
beta=(log(0.8)- log(A))/log(.5)
Q=0.5
#harm<-ifr_hic
#harm<-ifr_lic
pop <- hic_pop/sum(hic_pop)
#pop <- lic_pop/sum(lic_pop)

Q_temp<-Q

for (i in 1:length(harm)){
  x<-( Q_temp*harm[1:(length(harm)-i +1) ]^(1/(1-beta)) / sum(pop[1:(length(harm)-i +1)]*harm[1:(length(harm)-i +1)]^(1/(1-beta)) ) )
  print(x)
  i_star<-i
  if( max( A*x^beta) <0.95 ) {
    
    break
  }
  Q_temp<-Q_temp-pop[(length(harm)-i +1)]
}

for (i in 1:i_star){
  x<-( Q_temp*harm[i:(length(harm)-i_star+ 1) ]^(1/(1-beta)) / sum(pop[i:(length(harm)-i_star+ 1)]*harm[i:(length(harm)-i_star+ 1)]^(1/(1-beta)) ) )
  print(x)
  i_min<-i
  if( min( A*x^beta) >0.5 ) {
    break
  }
}

x_star<- c( rep(0,i_min-1),x,rep(1,i_star-1))
if (i_star==9){x_star<-rep(1,i_min)}
x_star


