########################
## Shape scientific data
########################

area_df <- loc_x %>%
  group_by(ICESNAME) %>%
  dplyr::summarise(area=n()) %>%
  mutate(perc_area = area/(sum(area))) %>%
  mutate(ICESNAME = factor(ICESNAME)) %>%
  arrange(ICESNAME)

## Sequential sampling (not that there is no preferential sampling)
# First select a statistical rectangle
sample_full_rect <- F
zonesize <- 0.2
n_zone <- 1 # c(1,3,5)
samp_df <- do.call(rbind,lapply(1:n_seq_com, function(j){
  
  ## Select a statistical rectangle
  rect_stat <- sample(x = area_df$ICESNAME,
                      size = 1,
                      replace=T,
                      prob =  area_df$perc_area)
  
  loc_x_rect <- loc_x %>%
    filter(ICESNAME == rect_stat)
  
  if(sample_full_rect == T){
    
    ## Sample within the statistical rectangle
    samp_cell <- sample(x = loc_x_rect$cell,
                        size = n_ping_per_seq,
                        replace=T)
    
  }else{
    
    samp_cell <- do.call(c,lapply(1:n_zone, function(j){
      ## Sample within the statistical rectangle
      center <- sample(x = loc_x_rect$cell,
                       size = 1,
                       replace=T)
      
      loc_x_zone <- loc_x_rect[which(loc_x_rect$x < loc_x_rect$x[which(loc_x_rect$cell == center)] + zonesize &
                                       loc_x_rect$x > loc_x_rect$x[which(loc_x_rect$cell == center)] - zonesize &
                                       loc_x_rect$y < loc_x_rect$y[which(loc_x_rect$cell == center)] + zonesize &
                                       loc_x_rect$y > loc_x_rect$y[which(loc_x_rect$cell == center)] - zonesize),]
      
      ## Sample within the statistical rectangle
      samp_cell <- sample(x = loc_x_zone$cell,
                          size = n_ping_per_seq/n_zone,
                          replace=T)
    }))
    
    
  }
  
  seq_id <- rep(j,n_ping_per_seq)
  df <- data.frame(cell=samp_cell,seq_id=seq_id,rect_stat)
  return(df)
  
}))

index_i <- samp_df$cell
boats_i <- samp_df$seq_id

c_com_x <- samp_df %>% 
  dplyr::count(cell) %>%
  full_join(loc_x) %>%
  arrange(cell) %>%
  data.frame() %>% 
  mutate(n = ifelse(is.na(n), 0,n)) %>%
  dplyr::select(n) %>%
  c() %>% unlist() %>%
  as.matrix()

## Simulate catch observations (zero-inflated lognormal distribution)
q1 <- -1
SD_obs <- exp(0.3)
y_i <- do.call(c,lapply(1:length(index_i), function(j){
  exp_catch <- S_x[index_i[j]]
  prob_encount <- 1-exp(- exp(q1) * exp_catch)
  abs_i <- rbinom(1,1,prob_encount)
  if(abs_i>0){
    y_i <- exp(rnorm(1,mean = log(exp_catch/prob_encount) - SD_obs^2/2, sd = SD_obs))
  }else{
    y_i <- 0
  }
}))

## Reallocation
k_sim <- 1
if(k_sim==1){ ## Reallocate data
  mean_y_i <- aggregate(x = y_i, 
                        by = list(unique.values = boats_i),
                        FUN = mean)[,2]
  
  sum_y_i <- aggregate(x = y_i, 
                       by = list(unique.values = boats_i),
                       FUN = sum)[,2]
  
  y_i2 <- rep(NA,length(y_i))
  y_i2 <- do.call(c,lapply(1:length(y_i), function(k){
    y_i2[k] <- mean_y_i[boats_i[k]]
  }))
}else{ ## Do not reallocate data
  y_i2 <- y_i
}

# # check:
# x11()
# par(mfrow = c(3,1))
# plot((S_x[index_i]),log(y_i))
# plot((y_i2),(y_i))
# plot((S_x[index_i]),log(y_i2))

#-----------------------------------------------------------------------------------
##-------------------------------- Scientific data ---------------------------------
#-----------------------------------------------------------------------------------

index_sci_i <- c()
n_cells <- nrow(loc_x)
n_samp_sci <- 75
nb_hauls_strata <- loc_x %>% 
  mutate(value=1) %>%
  group_by(strata) %>%
  dplyr::summarise(value = sum(value)) %>%
  mutate(hauls = round(value/n_cells*n_samp_sci))

# positions of hauls (straitified random sampling)
index_sci_i <- do.call(c,lapply(1:nrow(nb_hauls_strata),function(j){
  index_sci_i <- c(index_sci_i,sample(loc_x$cell[which(loc_x$strata == nb_hauls_strata$strata[j])], size=nb_hauls_strata$hauls[which(nb_hauls_strata$strata == nb_hauls_strata$strata[j])],replace=FALSE))
}))

# c_sci_x = ifelse(1:n_cells %in% index_sci_i, 1, 0) # shots for scientific data
# # check sampling of sci data
# verif_samp_sci <- cbind(loc_x,c_sci_x)
# simu_plot_1 <- ggplot()+
#   geom_point(data = verif_samp_sci,aes(x=x,y=y,col=strata)) +
#   geom_point(data = verif_samp_sci[which(verif_samp_sci$c_sci_x > 0),],aes(x=x,y=y)) +
#   theme_bw()

# scientific catch data (zero-inflated lognormal)
q1_sci <- 0
logSigma_sci <- -0.2
y_sci_i <- do.call(c,lapply(1:length(index_sci_i), function(j){
  exp_catch <- S_x[index_sci_i[j]]
  prob_encount <- 1-exp(- exp(q1_sci) * exp_catch)
  
  abs_sci_i <- rbinom(1,1,prob_encount)
  if(abs_sci_i>0){
    y_sci_i <-  exp(rnorm(1,mean = log(exp_catch)-exp(logSigma_sci)^2/2, sd = exp(logSigma_sci)))
  }else{
    y_sci_i <- 0 
  }
  return(y_sci_i)
  
}))
