# FINDING VALUES FOR DELTA2 that fit DELTA1 -----
# d1_default <- c("2%" = 1460, "4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "49%" = 60, "100%" = 1)
d1_default <- c(1, 60, 90, 120, 150, 180, 270, 360, 540, 730, 1460)

# d1_default <- c("4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "49%" = 45, "100%" = 1)
d_default <- sapply(d1_default, function(d){
  print(d)
  bsl <- sr(apap_2d(pars_fdf_slow, d, delay_default), f = "2d_v2") %>% 
    rescale_rcs(pop, merge = T)
  # bsl_vaccinated <- bsl[c(30, 60, 90, 180, 360),"cumV",1]
  
  # d_start <- floor(d*.25)
  d_start <- 1
  bsl_vaccinated <- bsl[d_start:min(d, 360),"cumV",1]
  
  sse <- sapply(ceiling(d/2):d, function(d2){
    fdf <- sr(apap_2d(pars_fdf_slow, d2, delay_fdf), f = "2d_v2") %>% 
      rescale_rcs(pop, merge = T)
    fdf_v <- fdf[d_start:min(d, 360),"cumV",1]
    # fdf_v <- fdf[c(30, 60, 90, 180, 360),"cumV",1]
    sum((fdf_v - bsl_vaccinated)^2)
  })
  sse_h <- sapply(ceiling(d/2):d, function(d2){
    fdf <- sr(apap_2d(pars_fdf_slow, c(rep(d2, 6), rep(d,3)), delay_hybrid), f = "2d_v2") %>% 
      # fdf <- sr(apap_2d(pars_fdf_slow, d2, delay_hybrid), f = "2d_v2") %>% 
      rescale_rcs(pop, merge = T)
    fdf_v <- fdf[d_start:min(d, 360),"cumV",1]
    sum((fdf_v - bsl_vaccinated)^2)
  })
  c(d,
    which.min(sse) + ceiling(d/2) - 1,
    which.min(sse_h) + ceiling(d/2) - 1)
})

# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# [1,]    1   60   90  120  150  180  270  360  540   730  1460
# [2,]    1   46   69   93  118  145  230  318  481   650  1300
# [3,]    1   52   77  102  128  155  239  325  496   680  1460
# > 

save(d1_default, d_default, file = "results/fdf-deltas.Rdata")
