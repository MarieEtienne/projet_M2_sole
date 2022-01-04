
## Simulate observations
#-----------------------

# obs_model <- "gamma" # "gamma", "zerolognorm" (zero-inflated lognormal distribution)

par(mfrow = c(3,2))
for(set_sim in 1:3){
  
  sum_y_i <- c()
  log_pi_j <- c()
  sum_mu <- c()
  Var_D <- c()
  E_D <- c()
  Var_Dsup0 <- c()
  Sigma_D <- c()
  decl_j <- c()
  
  sum_shape <- c()
  
  mu <- runif(1,min=0,max=2)
  
  for(i in 1:n_sim){
    
    ## Data simulation
    if(obs_model == "zerolognorm"){
      y_i <- do.call(c,lapply(1:n_obs, function(j){
        exp_catch <- mu
        prob_encount <- 1-exp(- exp(q1) * exp_catch)
        abs_i <- rbinom(1,1,prob_encount)
        if(abs_i>0){
          y_i <- exp(rnorm(1,mean = log(exp_catch/prob_encount) - SD_obs^2/2, sd = SD_obs))
        }else{
          y_i <- 0
        }
      }))
    }else if(obs_model == "gamma"){
      y_i <- do.call(c,lapply(1:n_obs, function(j){
        exp_catch <- mu
        prob_encount <- 1-exp(- exp(q1) * exp_catch)
        abs_i <- rbinom(1,1,prob_encount)
        y_i <- rgamma(1, shape = 1/SD_obs^2, scale = exp_catch*SD_obs^2)
      }))
    }
    
    # y_i <- Data$y_i
    sum_y_i <- c(sum_y_i,sum(y_i))
    
    ## Moments
    if(obs_model == "zerolognorm"){
      log_pi_j <- c(log_pi_j,sum(rep(-1 * mu * exp(q1),n_obs)))
      
      pi_j <- exp(log_pi_j)
      p_i <- rep(exp(-1 * mu * exp(q1)),n_obs)
      
      if(sum(y_i) == 0){
        
        Var_D <- c(Var_D,0)
        sum_mu <- c(sum_mu,0)
        E_D <- c(E_D,0)
        Var_Dsup0 <- c(Var_Dsup0,0)
        Sigma_D <- c(Sigma_D,0)
        
      }else{
        
        sum_mu <- c(sum_mu,sum(rep(mu,n_obs)))
        Var_D <- c(Var_D,sum(mu^2 / (1 - p_i) * (exp(SD_obs^2) - (1 - p_i))))
        E_D = c(E_D, 1 / (1 - pi_j[i]) * sum_mu[i])
        Var_Dsup0 = c(Var_Dsup0 ,Var_D[i] / (1 -  pi_j[i]) - pi_j[i] * E_D[i]^2 / (1 - pi_j[i])^2)
        Sigma_D = c(Sigma_D,sqrt(log(Var_Dsup0[i] / E_D[i]^2 + 1)))
        
      }
      
      exp_decl <- E_D[i]
      prob_decl <- 1-pi_j[i]
      
      abs_i <- rbinom(1,1,prob_decl)
      if(abs_i>0){
        decl_j <- c(decl_j,exp(rnorm(1,mean = log(E_D[i]) - Sigma_D^2/2, sd = Sigma_D)))
      }else{
        decl_j <- c(decl_j,0)
      }
    }else if(obs_model == "gamma"){
      # shape = 1/CV^2;   scale = mean*CV^2
      
      sum_shape <- c(sum_shape,sum(rep(1/SD_obs^2,n_obs)))
      scale_par <- mu*SD_obs^2
      decl_j <- c(decl_j,rgamma(1, shape = sum_shape[i], scale = scale_par))
      
    }
      
  }
  
  hist(sum_y_i,main="",breaks=seq(min(sum_y_i,decl_j)-2*SD_obs,max(sum_y_i,decl_j)+2*SD_obs,length.out=20))
  if(set_sim == 1) title(main="Simulating punctual observations y_i and \n then summing over fishing sequence D_j",cex.main=1)
  
  # mtext(text = paste0("mu = ", round(mu,digits = 3)))
  hist(decl_j,main="",breaks=seq(min(sum_y_i,decl_j)-2*SD_obs,max(sum_y_i,decl_j)+2*SD_obs,length.out=20))
  if(set_sim == 1) title(main="Simulating declaration D_j with \n the aggregated observation model",cex.main=1)

}



# ## Simulated data
# #----------------
# ## Data simulation
# 
# q1 <- obj$par["q1"]
# SD_obs <- exp(obj$par["logSD_obs"])
# seq_i <- Data$seq_i
# sum_y_i <- c()
# log_pi_j <- c()
# sum_mu <- c()
# Var_D <- c()
# E_D <- c()
# Var_Dsup0 <- c()
# Sigma_D <- c()
# decl_j <- c()
# 
# for(i in (unique(Data$seq_i))){
#   
#   sum_y_i <- c(sum_y_i,sum(y_i[which(seq_i == i)]))
#   samp_cell <- index_i[which(seq_i == i)]
#   
#   # log_pi_j(seq_i(i)) += -1 * E_catch(i) * exp(q1)
#   log_pi_j <- c(log_pi_j,sum(-1 * S_x[samp_cell] * exp(q1)))
#   pi_j <- exp(log_pi_j)
#   p_i <- exp(-1 * S_x[samp_cell] * exp(q1))
#   
#   if(sum_y_i[i+1] == 0){
#     
#     Var_D <- c(Var_D,0)
#     sum_mu <- c(sum_mu,0)
#     E_D <- c(E_D,0)
#     Var_Dsup0 <- c(Var_Dsup0,0)
#     Sigma_D <- c(Sigma_D,0)
#     
#   }else{
#     
#     sum_mu <- c(sum_mu,sum(S_x[samp_cell]))
#     Var_D <- c(Var_D,sum(S_x[samp_cell]^2 / (1 - p_i) * (exp(SD_obs^2) - (1 - p_i))))
#     E_D = c(E_D, 1 / (1 - pi_j[i+1]) * sum_mu[i+1])
#     Var_Dsup0 = c(Var_Dsup0 ,Var_D[i+1] / (1 -  pi_j[i+1]) - pi_j[i+1] * E_D[i+1]^2 / (1 - pi_j[i+1])^2)
#     Sigma_D = c(Sigma_D,sqrt(log(Var_Dsup0[i+1] / E_D[i+1]^2 + 1)))
#     
#   }
#   
# }
# 
# 
# 
# 
# exp_decl <- E_D[i]
# prob_decl <- 1-pi_j[i]
# 
# abs_i <- rbinom(1,1,prob_decl)
# if(abs_i>0){
#   decl_j <- c(decl_j,exp(rnorm(1,mean = log(E_D[i]) - Sigma_D^2/2, sd = Sigma_D)))
# }else{
#   decl_j <- c(decl_j,0)
# }






# for(int i=0; i<n_ping_i; i++){
#   
#   E_catch(i) = S_x(index_i(i));
#   log_pi_j(seq_i(i)) += -1 * E_catch(i) * exp(q1);
#   
#   if(y_i(seq_i(i)) > 0){
#     
#     p_i(i) = 1 - exp(-1 * E_catch(i) * exp(q1));
#     sum_mu(seq_i(i)) += S_x(index_i(i));
#     Var_D(seq_i(i)) += square(E_catch(i)) / (1 - p_i(i)) * (exp(square(SD_obs)) - (1 - p_i(i)));
#     
#   }
#   
# }
# 
# pi_j = exp(log_pi_j);
# 
# 
# for(int j=0; j<n_i; j++){
# 
#   if( y_i(j) == 0 ){
# 
#     jnll_comp(0) -= log_pi_j(j);
# 
#   }
# 
#   if( y_i(j) > 0 ){
# 
#     E_D(j) = 1 / (1 - pi_j(j)) * sum_mu(j);
# 
#     Var_Dsup0(j) =  Report$Var_D / (1 - Report$pi_j) - Report$pi_j * Report$E_D^2 / (1 - Report$pi_j)^2 ;
# 
#     Sigma_D = sqrt(log(Var_Dsup0(j) / square(E_D(j)) + 1));
# 
#     jnll_comp(0) -= dlognorm(y_i(j),
#                              log(E_D(j)),
#                              Sigma_D,
#                              true);
# 
#   }
# 
# }