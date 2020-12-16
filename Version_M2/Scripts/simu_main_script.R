####################
## Simulation script
####################
# B. Alglave based on Conn et al. (2017)

source("Scripts/function/load_packages.R")

#-----------------------------------------------------------
#-------------- Simulation/estimation settings -------------
#-----------------------------------------------------------

# Simulation name
simu_name = "test"

#---------------------
# Simulation scenarios
#---------------------

###
## Latent field
###

# Grid dimension
grid_side_x <- 25
grid_side_y <- 25
grid_dim = c("x"=grid_side_x, "y"=grid_side_y)
n_cells = grid_dim[1]*grid_dim[2]

# Intercepts and covariates of the latent field
beta0 = 2 # intercept
beta = c(0,0,0,0,0) # covariates
# parameters covariates vector : 1st param : continuous variable, 4 nexts : strata_1, strata_2, strata_3, strata_4
#la premiere covariable est continue et a une structure spatiale
#afin de simuler cette corrélation spatiale on a besoin de deux paramètres :
# - la portée de cette corrélation (nu)(+ la portée est grande et + la distance sur 
#   laquelle les points sont corrélées entre eux est importante, et donc plus la 
#   distance est grande, plus il faut "s'éloigner" pour avoir des points indépendants
# - la variance marginale de ce processus (écart type SD_x)(cad à quel point ce signal varie)


# Correlation structure of the continuous covariate 
simu_tool <- "MaternCov" # 2 way to simulate latent fields : "RandomFields", "MaternCov"
range = sqrt(prod(grid_dim))/(5)*2 # --> MaternCov 
nu = 1  #portée de la covariable continue
SD_x = 0.5 # ecart type/ variance marginale du processus de la var continue

#effets spatiaux aléatoires
SD_delta = SD_eta= 0.5


###
## Commercial sampling process
###

# intercept pour l'intensite lambda
beta0_fb = 2 

###
## Observation process
###

## Scientific data
# Number of scientific data
n_samp_sci = 50
# observation error 
logSigma_sci = log(1)
# zero-inflation parameter
q1_sci <- 1
# Relative catchability
q2_sci <- 1
# parameter controling size of the strata
n_strate <- 9 

## Commercial data
# Number of commercial samples
n_samp_com = 3000
# observation error 
logSigma_com = log(1)
# zero-inflation parameter
q1_com <- 1
# Relative catchability
q2_com <- 1
# Levels of preferential sampling
b_set = c(0,1,3) 



#--------------------
# Model configuration
#--------------------

# Data sources feeding the model
Data_source = c("scientific_commercial","scientific_only","commercial_only")

# b estimated or b fixed to 0 in estimation
EM_set = c("fix_b","est_b")
EM = EM_set[2]

# if 1, sampling process is accounted for in estimation
Samp_process = 1

# TMB model version
TmbFile = "Scripts/"

# if F, compute uncertainty for parameters estimates
ignore.uncertainty = T

## Loop indices
# Index for Results
counter <- 1 
# index of the first iteration
i0 <- 1
# Number of simulation
n_sim = 100

RandomSeed = 123456


#-------------------------------------------------------------
#-------------- Dataframes and model compilation -------------
#-------------------------------------------------------------

## Results dataframe
n_cov = 5
colnames_Results <- c("counter","sim","b_true","Data_source","type_b","Alpha","b_est","ObsMod","Sigma_com_true","Sigma_sci_true","Sigma_com","Sigma_sci","n_samp_com","n_samp_sci","N_true","N_est","SD_N","Convergence","LogLik","MSPE_S","k")

# add_cov_abs <- paste0("BiasBetaj_abs_",c(1:n_cov+1))
# add_cov_pos <- paste0("BiasBetaj_pos_",c(1:n_cov+1))
# colnames_Results <- c(colnames_Results,add_cov_abs,add_cov_pos)

Results = data.frame(matrix(NA,1,length(colnames_Results)))
colnames(Results)=colnames_Results

## list for simulated parameters and parameters estimates
List_param <- list()


####### Compile TMB model ########

# https://kaskr.github.io/adcomp/Introduction.html : for TMB details
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/troubleshooting.html 

TMB::compile(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple.cpp"),"-O1 -g",DLLFLAGS="")

# If problems with the PATH
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) # if "C:/Rtools/bin" is not in the PATH


#-------------------------------------------------------------
#------------------- Simulation/Estimation loop --------------
#-------------------------------------------------------------

################
# source Scripts
################
## Simulation loop function : simu_commercial_scientific()
source("Scripts/function/commercial_scientific_14.R")
## Simulate Matérn field : sim_GF_Matern()
source("Scripts/function/sim_GF_Matern.R")
## Fit model : fit_IM()
source("Scripts/function/fit_IM.R")
## Simulate latent field : simu_latent_field()
source("Scripts/function/simu_latent_field.R")
## Simulate scientific data : simu_scientific_data()
source("Scripts/function/simu_scientific_data.R")
## Simulate commercial data : simu_commercial_data()
source("Scripts/function/simu_commercial_data.R")

# file name for savinf outputs
Start_time.tot = Sys.time()
Start_time.tot_2 <- str_replace_all(Start_time.tot, " ", "_")
Start_time.tot_2 <- str_replace_all(Start_time.tot_2, ":", "_")
simu_file <- paste0("results/com_x_sci_data_14_scientific_commercial_simple-",Start_time.tot_2,"_",simu_name,"/")

# When simu crashes --> param to re-run the loop.
# Before re-runing the loop load the last Results_loop list (Results/Simu/)
restart_after_crash = F
if(restart_after_crash == T){
  counter <- max(Results$counter,na.rm=T)
  i0 <- max(Results$sim,na.rm=T) + 1
  n_sim <- n_sim
  Results <- Results
  List_param <- Results_loop$List_param
  simu_file <- "results/com_x_sci_data_14/com_x_sci_data_14-2020-08-29_11_59_26_Lognormal_Nsamp"
}

# load TMB model
dyn.load( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )
#dyn.unload( dynlib(paste0(TmbFile,"inst/executables/",Version,"_scientific_commercial") ) )
## loop
for(i in i0:n_sim){

    res <- simu_commercial_scientific(Results,
                                      simu_file,
                                      grid_dim,
                                      n_cells,
                                      beta0,
                                      beta,
                                      range,
                                      nu,
                                      SD_x,
                                      SD_delta,
                                      SD_eta,
                                      n_samp_sci,
                                      logSigma_sci,
                                      q1_sci,
                                      q2_sci,
                                      n_strate,
                                      n_samp_com,
                                      logSigma_com,
                                      q1_com,
                                      q2_com,
                                      b_set,
                                      Data_source,
                                      Samp_process,
                                      EM,
                                      RandomSeed,
                                      Version,
                                      TmbFile,
                                      ignore.uncertainty,
                                      counter,
                                      i,
                                      n_sim)
    
    Results <- res[[1]]
    List_param <- res[[2]]
    counter <- res[[3]]

}
dyn.unload( dynlib(paste0(TmbFile,"inst/executables/com_x_sci_data_14_scientific_commercial_simple") ) )

#-------------------------------------------------------------
#------------------------ Plot results -----------------------
#-------------------------------------------------------------

## Function for plotting elative bias of abundance, of b and MSPE : Plot_Results()
source("Scripts/function/plot_simu.R")

# Results <- Results_loop$Results
# Strue_df <- Results_loop$Strue_df
# PE_df <- Results_loop$PE_df
# SD.S_df <- Results_loop$SD.S_df

Plot_results_list <- Plot_Results(Results,b_set)
grid.arrange(Plot_results_list[[1]],Plot_results_list[[2]],Plot_results_list[[3]],ncol=3)


# ## (DEPRECATED) Plot parameters estimate (relative) bias
# ## See "in_progress" folder
# source("Scripts/rel_bias_param.R")
# # Results_loop file
# file <- ""
# load(file)
# RelBias_param(Results_loop)


