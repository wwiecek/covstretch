rescale_rcs <- function(y, pop_sizes=rep(1, dim(y)[3]), merge = FALSE) {
  y_new <- y*rep(pop_sizes, each = dim(y)[1]*dim(y)[2])
  if(merge)
    y_new <- replicate(1, apply(y_new, c(1,2), sum))
  y_new
}

rescale_and_bind <- function(ll, pop) {
  lapply(ll, rescale_rcs, pop_sizes = pop, merge = T) %>% abind::abind(along=3)
}

plot_rcs <- function(y, compartment = "R", 
                     shade_weeks = c(0,0), 
                     start_date = as.Date("01-01-2021", format="%d-%m-%Y"), 
                     end_date = start_date + dim(y)[1] - 30,
                     lab_type = "", overlay_data = NULL,
                     long_names = ln,
                     ncol = 2) {
  gg_data <- as.data.frame(y[,compartment,]) %>%
    rownames_to_column("time") %>%
    mutate(time = as.numeric(time)) %>%
    gather(group_name, value, -time)
  
  if(length(compartment) > 1){
    if(dim(y)[3] > 1)
      gg_data <- separate(gg_data, group_name, c("compartment", "group_name"), sep = "\\.") 
    else
      gg_data$compartment <- gg_data$group_name
  }
  
  if(!is.null(start_date)){
    gg_data$time <- as.Date(gg_data$time, origin = start_date)
    gg_data <- gg_data[gg_data$time <= as.Date(end_date),]
  }
  
  if(max(gg_data$value, na.rm=T) > 1)
    xl <- "N individuals"
  else
    xl <- "Proportion"
  
  if(!is.null(long_names) && !is.null(gg_data$compartment))
    gg_data$compartment <- factor(gg_data$compartment, levels = names(long_names), labels = long_names)
    
  gg <- ggplot() + 
    # geom_rect(aes(xmin = as.Date(7*shade_weeks[1], origin = start_date), 
    # xmax = as.Date(7*shade_weeks[2], origin = start_date),
    # ymin=0, ymax=Inf), fill = "skyblue", alpha = .2) +
    geom_line(data = gg_data %>% mutate(group_name=factor(group_name,levels=rev(unique(gg_data$group_name))),
                                        compartment=factor(compartment,levels=rev(unique(gg_data$compartment)))),
              aes(x=time, y=value, group=group_name, color=group_name)) +
    {if(!is.null(start_date)) scale_x_date(limits = c(as.Date(start_date), as.Date(end_date)))} +
    # {if(length(compartment) == 1) labs(y = paste0(compartment_names_trt[compartment], lab_type))} +
    {if(dim(y)[3] == 1) theme(legend.position = "none")} +
    {if(length(compartment) > 1) facet_wrap( ~ compartment, ncol = ncol,
                                             scales = "free")} +
    scale_color_discrete(name = "") +
    ylab(xl)
  
  if(!is.null(overlay_data) && nrow(overlay_data) > 0){
    # overlay_data$time <- as.Date(overlay_data$time, origin = start_date)
    gg <- gg +
      geom_point(aes(x=time, y=value), data = overlay_data) +
      geom_line(aes(x=time, y=value), data = overlay_data)
  }
  gg
}


