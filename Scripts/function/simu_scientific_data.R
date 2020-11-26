###################################
## Simulate or load scientific data 
###################################

#' @title simu_scientific_data()
#' 
#' @param loc_x dataframe with grid coordinate and strata
#' @param grid_dim grid dimension
#' @param Strue_x spatial abundance data
#' @param scientific_data already simulated data for the scientific data
#' @param DataObs Observation model (1 : zero-inflated gamma, 2 : zero-inflated lognormal, 3 : lognormal)
#' @param zero.infl_model Zero-inflated model parametrisation
#' @param n_samp_sci number number of scientific samples
#' @param logSigma_sci observation error for the scientific data (log scale)
#' @param q1_sci parameter related to the zero-inflation for scientific data
#' @param q2_sci relative catchability of the scientific data
#' @param sci_sampling type of sampling plan for scientific data ("random_stratified";"fixed")
#' @param n_strate parameter for the size of the strata (n=9)
#' 
#' 
#' @return index_sci_i : cells sampled by the scientific survey
#' @return c_sci_x : vector filled with 0/1. If 1 : the cell has been sampled. If 0 :  the cell has not been sampled. 
#' @return y_sci_i : scientific observations 
#' 
#' @author B. Alglave


simu_scientific_data <- function(loc_x,
                                 grid_dim,
                                 Strue_x,
                                 zero.infl_model,
                                 n_samp_sci,
                                 logSigma_sci,
                                 q1_sci,
                                 q2_sci,
                                 n_strate,
                                 scientific_data){

  index_sci_i <- c()
  
  # ## Uniform sampling
  # index_sci_i <- sample(loc_x$cell,size=n_samp_sci)

  ## Random stratified sampling
  loc_x %>%
    tidyr::pivot_longer(cols = starts_with("strata"),names_to = "strata") %>%
    data.frame() %>% filter(value > 0) -> loc_x_2
  # nb of samples in each strata

  loc_x_2 %>% group_by(strata) %>%
    dplyr::summarise(value = sum(value)) %>%
    mutate(hauls = round(value/n_cells*n_samp_sci)) -> nb_hauls_strata

  index_sci_i <- do.call(c,lapply(1:nrow(nb_hauls_strata),function(j){
    index_sci_i <- c(index_sci_i,sample(loc_x_2$cell[which(loc_x_2$strata == nb_hauls_strata$strata[j])], size=nb_hauls_strata$hauls[which(nb_hauls_strata$strata == nb_hauls_strata$strata[j])],replace=FALSE))
  }))
  
  c_sci_x = ifelse(1:prod(grid_dim) %in% index_sci_i, 1, 0) # shots for scientific data
  
  # # check sampling of sci data
  # verif_samp_sci <- cbind(loc_x,c_sci_x)
  # verif_samp_sci %>%
  #   tidyr::pivot_longer(cols = c(starts_with("strata")),names_to = "strata") %>%
  #   data.frame() %>% filter(value > 0) -> verif_samp_sci
  # simu_plot_1 <- ggplot()+
  #   geom_point(data = verif_samp_sci,aes(x=x,y=y,col=strata)) +
  #   geom_point(data = verif_samp_sci[which(verif_samp_sci$c_sci_x > 0),],aes(x=x,y=y)) +
  #   theme_bw()
  
  # scientific data (zero-inflated)
  y_sci_i <- do.call(c,lapply(1:length(index_sci_i), function(j){
    exp_catch <- q2_sci * Strue_x[index_sci_i[j]] # expected catch
    # proba of encounter
    prob_encount <- 1-exp(- exp(q1_sci[1]) * exp_catch)
    
    abs_sci_i <- rbinom(1,1,prob_encount)
    if(abs_sci_i>0){
      # y_sci_i <- rlnorm(1,meanlog = log(exp_catch), sdlog = exp(logSigma_sci))
      # shape <- exp(logSigma_sci)^-2
      # y_sci_i <-  rgamma(1,shape=shape,scale= exp_catch * exp(logSigma_sci)^2)
      y_sci_i <-  exp(rnorm(1,mean = log(exp_catch)-exp(logSigma_sci)^2/2, sd = exp(logSigma_sci)))
    }else{
      y_sci_i <- 0 
    }
    return(y_sci_i)
    
  })) # q2_sci * Strue_x : expected catch     |     q1_sci : parameter linking number of zero and density for scientific data
  
  res <- list(index_sci_i = index_sci_i,
              c_sci_x = c_sci_x,
              y_sci_i=y_sci_i)
  
  return(res)
  
}


# scientific_data <- list()
# 
# # run the loop
# for(i in 1:100){
#   print(i)
#   set.seed( RandomSeed + i ) # i = 10
#   Strue_x <- latent_field[[i]]$Strue_x
#   scientific_data[[i]] <- simu_scientific_data(loc_x,grid_dim,sci_sampling,n_samp_sci,q1_sci,q2_sci,Strue_x,zero.infl_model,logSigma_sci)
# }
# 
# save(data = scientific_data,file = "C:/Users/test/Documents/GitHub/phd-baptiste_alglave/Results/Simu/simulated_data/scientific_data.RData")

