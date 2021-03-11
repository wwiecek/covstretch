# Function for deriving age prioritised vaccination parameters -----

# apap = adjust parameters for age prioritisation
# Derive values of delta and ta such that 
# len = length of vaccination campaign (in days), i.e. time to 100% vaccinated
# vhes = when vaccine hesitancy kicks in (in each age group)
# vsupply = total supply (fraction of total population)
# group_seq = whether to prioritise all of the groups sequentially (80+, 70-80, ..., "new" way)
#             or just the top 3 age groups before the rest ("old" way)

delay_by_age <- c(1,1,1,1,1,1,0.0001,0.0001,0.0001) #I use 0.0001 instead of 0 to avoid 0*Inf
avail_by_age <- c(0,0,1,1,1,1,1,1,1)
prop_young <- sum(pop[1:2])
prop_old <- sum(pop[7:9])
prop_adults <- 1-prop_old-prop_young
prop_all <- c(rep(prop_young, 2), rep(1-prop_young-prop_old, 4), rep(prop_old, 3))

apap_2d <- function(pars, len, 
                    d2 = 18,
                    vhes = .8, vsupply = default_supply_ceiling,
                    group_seq=default_group_seq) {
  if(length(d2) == 1)
    d2 <- rep(d2, Ngroups)
  if(length(len) == 1)
    len <- rep(len, 9)
  
  allocation_vector <- vac_top_p(vsupply/vhes, pop)
  
  if (!group_seq){
    ta <- 10 + len*prop_old*vhes*delay_by_age
    delta1 <- avail_by_age/len/prop_all
  } else {
    ta <- 10 + c(rev(cumsum(rev(len*pop*vhes*(allocation_vector+0.00001))))[-1],0)
    delta1 <- avail_by_age/len/pop
  }
  
  list_modify(pars, 
              vstop = vhes*allocation_vector,
              delta2 = 1/d2,
              ta = ta,
              delta1 = delta1)
}




# Same as above but for the model with 2 vaccines, not 2 doses -----

apap_2v <- function(pars, len, 
                    expand_from = Inf,
                    expansion_factor = 2,
                    fractional_dose = rep(1, length(pop)),
                    switch=Inf, 
                    delay = 10, 
                    vhes = .8, vsupply = default_supply_ceiling,
                    group_seq=default_group_seq) {
  allocation_vector <- vac_top_p(vsupply/vhes, pop)
  d1 <- avail_by_age/len/prop_all
  
  if(length(delay)==1){
    ts1 <- rep(switch+delay, Ngroups)
  } else {
    ts1 <- switch+delay
  }
  
  if (!group_seq){
    d1 <- avail_by_age/len/prop_all/fractional_dose
    t1 <- delay + len*prop_old*vhes*delay_by_age*fractional_dose
    
    # Find when the vaccinations would be completed in priority group if there was 
    # an expansion along the way (this is to adjust t1 in non-priority groups)
    
    if(expand_from < t1[1]){
      v <- 1/360/prop_old
      # x = e + (s - (e-d)v/Lv)
      t1[1:6] <- expand_from + (vhes - (expand_from - delay)*v)/(expansion_factor*v)
    }
  } else {
    d1 <- avail_by_age/len/pop/fractional_dose
    time_to_vaccinate_k <- len*pop*vhes*(allocation_vector+0.00001)*fractional_dose
    t1 <- delay + c(rev(cumsum(rev(time_to_vaccinate_k)))[-1],0)
    
    if ((expand_from>delay[length(delay)])&(expand_from!=Inf)){
      t1.expand <- t1 - expand_from
      t1.expand[t1.expand<0] <- 0
      t1.expand[t1.expand>min(t1.expand[t1.expand>0])] <- 1
      switch_index <- which(!t1.expand %in% c(0,1))
      t1.expand[switch_index] <- t1.expand[switch_index]/(time_to_vaccinate_k)[switch_index+1]
      pop.expand <- t1.expand*c((pop*vhes*allocation_vector)[-1],0.0001)
      t1.delta <- pop.expand*len*(1-1/expansion_factor)
      t1 <- t1 - rev(cumsum(rev(t1.delta)))
    }  
    if ((expand_from<=delay[length(delay)])){
      t1 <- delay + c(rev(cumsum(rev((len/expansion_factor)*pop*vhes*(allocation_vector+0.00001))))[2:length(pop)],0)
    }
  }
  
  list_modify(pars, 
              vstop = vhes*allocation_vector,
              ta1 = t1,
              tmore1 = rep(expand_from, Ngroups),
              ta2 = sapply(t1, function(x) max(x, switch+delay)),
              ts1 = ts1,
              delta1 = d1,
              delta2 = d1)
}


# sr(apap_2d(pars_fdf_slow, 360, vsupply = .2), f = "2d_v2") %>% plot_rcs(c("cumV", "cumV1", "cumV2", "I"))
# sr(apap_2d(pars_fdf_slow, 360, vsupply = .15), f = "2d_v2") %>% rescale_rcs(merge=TRUE, pop) %>% plot_rcs(c("cumV", "cumV1", "cumV2", "I"))