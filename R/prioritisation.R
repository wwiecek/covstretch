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
                    vhes = .8, vsupply = default_supply_ceiling) {
  if(length(d2) == 1)
    d2 <- rep(d2, Ngroups)
  if(length(len) == 1)
    len <- rep(len, 9)
  
  allocation_vector <- vac_top_p(vsupply/vhes, pop)
    
  list_modify(pars, 
              vstop = vhes*allocation_vector,
              delta2 = 1/d2,
              ta = 10 + len*prop_old*vhes*delay_by_age,
              delta1 = avail_by_age/len/prop_all)
}

apap_2v <- function(pars, len, 
                    expand_from = Inf,
                    expansion_factor = 2,
                    switch=Inf, 
                    delay = 10, 
                    vhes = .8, vsupply = default_supply_ceiling) {
  d1 <- avail_by_age/len/prop_all
  t1 <- delay + len*prop_old*vhes*delay_by_age
  allocation_vector <- vac_top_p(vsupply/vhes, pop)
  
  # Find when the vaccinations would be completed in priority group if there was 
  # an expansion along the way (this is to adjust t1 in non-priority groups)
  if(expand_from < t1[1]){
    v <- 1/360/prop_old
    # x = e + (s - (e-d)v/Lv)
    t1[1:6] <- expand_from + (vhes - (expand_from - delay)*v)/(expansion_factor*v)
  }
  
  list_modify(pars, 
              vstop = vhes*allocation_vector,
              ta1 = t1,
              tmore1 = rep(expand_from, Ngroups),
              ta2 = sapply(t1, function(x) max(x, switch+delay)),
              ts1 = rep(switch+delay, Ngroups),
              delta1 = d1,
              delta2 = d1)
}


# sr(apap_2d(pars_fdf_slow, 360, vsupply = .2), f = "2d_v2") %>% plot_rcs(c("cumV", "cumV1", "cumV2", "I"))
# sr(apap_2d(pars_fdf_slow, 360, vsupply = .15), f = "2d_v2") %>% rescale_rcs(merge=TRUE, pop) %>% plot_rcs(c("cumV", "cumV1", "cumV2", "I"))
