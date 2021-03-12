# FINDING VALUES FOR DELTA2 that fit DELTA1 -----

# d1_default <- c("4%" = 730, "8%" = 360, "15%" = 180, "28%" = 90, "49%" = 45, "100%" = 1) fdf_speeds
fdf_deltas_raw <- lapply(as.list(fdf_speeds), function(d){
  print(d)
  bsl <- sr(apap_2d(pars_fdf_fast, d, delay_default), f = "2d_v2") %>% 
    rescale_rcs(pop, merge = T)
  # bsl_vaccinated <- bsl[c(30, 60, 90, 180, 360),"cumV",1]
  
  # d_start <- floor(d*.25)
  d_start <- 1
  bsl_vaccinated <- bsl[d_start:min(d, 360),"cumV",1]
  delta_vals <- ceiling(d/2):d
  
  sse <- sapply(delta_vals, function(d2){
    fdf <- sr(apap_2d(pars_fdf_fast, d2, delay_fdf), f = "2d_v2") %>% 
      rescale_rcs(pop, merge = T)
    fdf_v <- fdf[d_start:min(d, 360),"cumV",1]
    # fdf_v <- fdf[c(30, 60, 90, 180, 360),"cumV",1]
    sum((fdf_v - bsl_vaccinated)^2)
  })
  
  if (!all_k){
    sse_h <- sapply(delta_vals, function(d2){
      fdf <- sr(apap_2d(pars_fdf_fast, c(rep(d2, 6), rep(d,3)), delay_hybrid), f = "2d_v2") %>% 
        # fdf <- sr(apap_2d(pars_fdf_slow, d2, delay_hybrid), f = "2d_v2") %>% 
        rescale_rcs(pop, merge = T)
      fdf_v <- fdf[d_start:min(d, 360),"cumV",1]
      sum((fdf_v - bsl_vaccinated)^2)
    })
    
    data.frame(d, delta_vals, sse, sse_h)
  } else {
    
    for (k in c(1,2,3,4,5,6,7,8)){
      sse_h.tmp <- sapply(delta_vals, function(d2){
        fdf <- sr(apap_2d(pars_fdf_fast, c(rep(d2, k), rep(d,9-k)), delay_hybrid_k[,k]), f = "2d_v2") %>% 
          # fdf <- sr(apap_2d(pars_fdf_slow, d2, delay_hybrid), f = "2d_v2") %>% 
          rescale_rcs(pop, merge = T)
        fdf_v <- fdf[d_start:min(d, 360),"cumV",1]
        sum((fdf_v - bsl_vaccinated)^2)
      })
      sse_h.tmp <- data.frame(sse_h=sse_h.tmp)
      colnames(sse_h.tmp) <- paste0("sse_h",k)
      if (k==1){
        sse_h <- sse_h.tmp
      } else {
        sse_h <- cbind(sse_h,sse_h.tmp)
      }
    }
    data.frame(d, delta_vals, sse, sse_h)
  }
  
})

if (!all_k){
  fdf_deltas <- bind_rows(fdf_deltas_raw) %>%
    group_by(d) %>%
    summarise(d2 = delta_vals[which.min(sse)],
              d3 = delta_vals[which.min(sse_h)],
              sse = min(sse),
              sse_h = min(sse_h))
  
  save(fdf_speeds, fdf_deltas, file = "results/fdf-deltas.Rdata")
} else {
  fdf_deltas <- bind_rows(fdf_deltas_raw) %>%
    group_by(d) %>%
    summarise(d_fdf = delta_vals[which.min(sse)],
              d_sse_h1 = delta_vals[which.min(sse_h1)],
              d_sse_h2 = delta_vals[which.min(sse_h2)],
              d_sse_h3 = delta_vals[which.min(sse_h3)],
              d_sse_h4 = delta_vals[which.min(sse_h4)],
              d_sse_h5 = delta_vals[which.min(sse_h5)],
              d_sse_h6 = delta_vals[which.min(sse_h6)],
              d_sse_h7 = delta_vals[which.min(sse_h7)],
              d_sse_h8 = delta_vals[which.min(sse_h8)],
              sse = min(sse),
              sse_h1 = min(sse_h1),
              sse_h2 = min(sse_h2),
              sse_h3 = min(sse_h3),
              sse_h4 = min(sse_h4),
              sse_h5 = min(sse_h5),
              sse_h6 = min(sse_h6),
              sse_h7 = min(sse_h7),
              sse_h8 = min(sse_h8))
  
  save(fdf_speeds, fdf_deltas, file = "results/fdf-deltas-allk.Rdata")
}

# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# [1,]    1   60   90  120  150  180  270  360  540   730  1460
# [2,]    1   46   69   93  118  145  230  318  481   650  1300
# [3,]    1   52   77  102  128  155  239  325  496   680  1460
# > 
