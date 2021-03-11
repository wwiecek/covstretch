# Function for deriving age prioritised vaccination parameters -----

# apap = adjust parameters for age prioritisation
# Derive values of delta and ta such that 
# len = len of vaccination campaign (in days), i.e. time to 100%
# vhes = when vaccine hesitancy kicks in (in each age group)
# vsupply = total supply (fraction of total population)

delay_by_age <- c(1,1,1,1,1,1,0.0001,0.0001,0.0001) #I use .0001 instead of 0 to avoid 0*Inf
avail_by_age <- c(0,0,1,1,1,1,1,1,1)
prop_young <- sum(pop[1:2])
prop_old <- sum(pop[7:9])
prop_adults <- 1-prop_old-prop_young
prop_all <- c(rep(prop_young, 2), rep(1-prop_young-prop_old, 4), rep(prop_old, 3))

apap_2d <- function(pars, len, 
                    d2 = 18,
                    vhes = .8, vsupply = default_supply_ceiling,
                    group_seq=FALSE) {
  if(length(d2) == 1)
    d2 <- rep(d2, Ngroups)
  if(length(len) == 1)
    len <- rep(len, 9)
  
  allocation_vector <- vac_top_p(vsupply/vhes, pop)
  
  if (group_seq==FALSE){
    ta <- 10 + len*prop_old*vhes*delay_by_age
  } else {
    ta <- 10 + c(rev(cumsum(rev(len*pop*vhes*(allocation_vector+0.00001))))[2:length(pop)],0)
  }
  
  list_modify(pars, 
              vstop = vhes*allocation_vector,
              delta2 = 1/d2,
              ta = ta,
              delta1 = avail_by_age/len/prop_all)
}

apap_2v <- function(pars, len, 
                    expand_from = Inf,
                    expansion_factor = 2,
                    switch=Inf, 
                    delay = 10, 
                    vhes = .8, vsupply = default_supply_ceiling,
                    group_seq=FALSE) {
  allocation_vector <- vac_top_p(vsupply/vhes, pop)
  d1 <- avail_by_age/len/prop_all
  
  if(length(delay)==1){
    ts1 <- rep(switch+delay, Ngroups)
  } else {
    ts1 <- switch+delay
  }
  
  if (group_seq==FALSE){
    
    t1 <- delay + len*prop_old*vhes*delay_by_age
    
    # Find when the vaccinations would be completed in priority group if there was 
    # an expansion along the way (this is to adjust t1 in non-priority groups)
    if(expand_from < t1[1]){
      v <- 1/360/prop_old
      # x = e + (s - (e-d)v/Lv)
      t1[1:6] <- expand_from + (vhes - (expand_from - delay)*v)/(expansion_factor*v)
    }
  } else {
    
    t1 <- delay + c(rev(cumsum(rev(len*pop*vhes*(allocation_vector+0.00001))))[2:length(pop)],0)
    
    if (expand_from>delay[length(delay)]){
      t1.expand <- t1 - expand_from
      t1.expand[t1.expand<0] <- 0
      t1.expand[t1.expand>min(t1.expand[t1.expand>0])] <- 1
      switch_index <- which(!t1.expand %in% c(0,1))
      t1.expand[switch_index] <- t1.expand[switch_index]/(len*pop*vhes*(allocation_vector+0.00001))[switch_index+1]
      pop.expand <- t1.expand*c((pop*vhes*allocation_vector)[2:length(pop)],0.0001)
      
      t1.delta <- pop.expand*len*(1-1/expansion_factor)
      t1 <- t1 - rev(cumsum(rev(t1.delta)))
    } else {
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
