#####################################
## Simulate or load latent field data 
#####################################

#' @title simu_latent_field()
#' 
#' @param loc_x dataframe with grid coordinate and strata
#' @param latent_fields_simu if == "simulate" : simulate data, if == "load_from_RData" : load already simulated data, to have the same data between different simulation loops (here `latent_field` and `scientific_data`)
#' @param latent_field already simulated data for the abundance distribution equation
#' @param beta0 intercept of the abundance distribution equation
#' @param beta  fixed parameters of the abundance distribution equation
#' @param simu_tool simulation functions ("RandomFields" or "MaternCov")
#' @param SpatialScale spatial scale parameter. Must be provided for the `RandomFields` package
#' @param range range parameter. Must be provided for the "MaternCov" package.
#' @param nu nu parameter. Must be provided for the "MaternCov" package.
#' @param SD_x Standard deviation for the simulation of covariate related to the abundance distribution equation.
#' @param SD_delta Standard deviation for the random effect delta (abundance distribution).
#' @param SD_xfb Standard deviation for the simulation of covariate related to the sampling process.
#' @param SD_eta Standard devioation for the random effect eta (sampling process).
#' 
#' @return Cov_x : design matrix (matrix for covariates levles / values)
#' @return Strue_x Latent field values
#' @return beta : fixed parameters for the abundance distribution equation
#' @return delta_x : random effect for the abundance distribution equation
#' @return xfb_x : fixed parameters for the sampling process equation
#' @return eta_x : random effect for the sampling process equation
#' 
#' @author B. Alglave

# for example : run main_script and beginning of commercial_scientific_14

simu_latent_field <- function(loc_x,
                              latent_fields_simu,
                              latent_field,
                              scientific_data,
                              beta0,
                              beta,
                              simu_tool,
                              range,
                              nu,
                              SD_x,
                              SD_delta){
  
  # 1st param : continuous variable, 4 nexts : strata_1, strata_2, strata_3, strata_4
  beta[1] = runif(1,-.5,.5)

  # Simulate using Matérn covariance
  Cov_x = sim_GF_Matern(loc_x, nu, range, SD_x^2)[[1]]$y

  # Create design matrix for covariates
  Cov_x <- as.matrix(cbind(Cov_x,loc_x[,which(str_detect(colnames(loc_x),"strata"))]))
  
  # Create random noise
  delta_x <- rnorm(nrow(Cov_x),mean = 0, sd = SD_delta)
  
  # Total abundance
  Strue_x = exp(beta0 + as.matrix(Cov_x) %*% beta + delta_x)
  res <- list(Cov_x = Cov_x,
              Strue_x = Strue_x,
              beta = beta,
              delta_x = delta_x)
  return(res)
}

#-----------------------------------------------------------------------------------
# Create data to save latent fields in order to fit models on the same latent fields
#-----------------------------------------------------------------------------------

## Construct grid

# # create grid and cells + strata
# n_cells <- grid_dim['x'] * grid_dim['y']
# loc_x = expand.grid( "x"=1:grid_dim['x'], "y"=1:grid_dim['y'])
# loc_x = cbind(loc_x,cell = 1:n_cells)
# diff.dim_y_x <- grid_dim['y'] - grid_dim['x']
# loc_x %>%
#   mutate(strata_1 = ifelse(x <= n_strate & y > grid_dim['y'] - n_strate , 1, 0)) %>%
#   mutate(strata_4 = ifelse(x >= grid_dim['x'] - n_strate & y <= n_strate + diff.dim_y_x, 1, 0), 1, 0) %>%
#   mutate(strata_3 = ifelse((x < grid_dim['x'] - n_strate & y <= grid_dim['y'] - n_strate & x + y <= grid_dim['y'] + 1),1,0)) %>%
#   mutate(strata_2 = ifelse(strata_1 == 0 & strata_4 == 0 & strata_3 == 0, 1, 0)) %>%
#   dplyr::select(x,y,cell,strata_1,strata_2,strata_3,strata_4) -> loc_x
# 
# # build output list
# list_res <- list()
# 
# # run the loop
# for(i in 1:100){
#   print(i)
#   set.seed( RandomSeed + i ) # i = 10
# 
#   simu_latent_field_outputs <- simu_latent_field(beta,Spatial_sim_model,SD_delta,SpatialScale,loc_x)
#   list_res[[i]] <- simu_latent_field_outputs
# 
# }
# latent_field <- list_res
# save(data = latent_field,file = "C:/Users/test/Documents/GitHub/phd-baptiste_alglave/Results/Simu/simulated_data/latent_field.RData")
