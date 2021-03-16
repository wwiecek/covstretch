library(nleqslv)

library(tidyverse)
load("data/default_inputs.Rdata")

# Demographics (for comparing HIC vs LIC)
hic_pop <- pbc_spread[countries["High-income countries"],] %>% as.numeric()
lic_pop <- pbc_spread[countries["Low-income countries"],] %>% as.numeric()
ifr_hic <- c(0.002, 0.006, 0.03, 0.08, 0.15, 0.60, 2.2, 5.1, 9.3)/100
ifr_lic <- ifr_hic*(3.2/2)^(5:(-3))


#concave function
A=3.49706
B=1.74853
C=-0.798528
lbound<-0.2427465959040336529990383979950291386649616526142562882369696689

harm<-ifr_hic
#harm<-ifr_lic
pop <- hic_pop/sum(hic_pop)
#pop <- lic_pop/sum(lic_pop)



model<-function(p,Q,pop,harm){
  
  x<-c(rep(1,length(pop)),0)
  system_eq <- function(x) {
    y <- numeric(length(pop)+1)
    y[1]<- sum(pop*x[1:(length(pop))]) -Q
    y[2:length(y)]<-((B/sqrt(x[1:(length(pop))])) -B)*harm - p-x[length(pop)+1]
    return(y)
  }
  system_eq(x)
  
  out<-nleqslv(x, system_eq, control=list(btol=.01))
  #system_eq(out$x)
  out$x[1:(length(pop))]
  
  
}

optimal_policy<-function(pop,harm,Q,p=0){
  out<-model(p=0,Q,pop,harm)
  i<-1
  while (min(out)<lbound){
    i<-i+1
    out<-model(p=0,Q,pop[i:length(pop)],harm[i:length(pop)])
  }
  c(rep(0,length(pop)-length(out)), out)
}

lic_results<-data_frame(
  "Age Group" = c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80+"),
  "Population Share" = lic_pop/sum(lic_pop),
  "Infection Fatality Rate " = ifr_lic,
  "1.0"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,1.0),
  "0.9"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.9),
  "0.8"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.8),
  "0.7"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.7),
  "0.6"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.6),
  "0.5"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.5),
  "0.4"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.4),
  "0.3"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.3),
  "0.2"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.2),
  "0.1"= optimal_policy(lic_pop/sum(lic_pop),ifr_lic,0.1))

hic_results<-data_frame(
  "Age Group" = c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80+"),
  "Population Share" = hic_pop/sum(hic_pop),
  "Infection Fatality Rate " = ifr_hic,
  "1.0"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,1.0),
  "0.9"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.9),
  "0.8"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.8),
  "0.7"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.7),
  "0.6"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.6),
  "0.5"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.5),
  "0.4"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.4),
  "0.3"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.3),
  "0.2"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.2),
  "0.1"= optimal_policy(hic_pop/sum(hic_pop),ifr_hic,0.1))
lic_results
hic_results
