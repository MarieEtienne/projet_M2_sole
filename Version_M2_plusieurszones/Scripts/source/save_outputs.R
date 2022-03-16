###############################
## Save outputs of simumlations
###############################
## B. Alglave

simu.est.data <- list(S_x=S_x,
                  index_i=index_i,
                  y_sci_i=y_sci_i,
                  y_i=y_i,
                  boats_i=boats_i,
                  cov_x=cov_x,
                  index_sci_i=index_sci_i,
                  intercept=intercept,
                  beta1=beta1,
                  SD_delta=SD_delta,
                  range_delta=range_delta,
                  q1_com=q1,
                  Sigma_com=SD_obs,
                  q1_sci=q1_sci,
                  Sigma_sci=Sigma_sci,
                  nu = nu,
                  range_cov=range_cov,
                  SD_cov=SD_cov,
                  q1=q1,
                  n_seq_com=n_seq_com,
                  n_ping_per_seq=n_ping_per_seq,
                  zonesize=zonesize,
                  n_zone=n_zone,
                  sample_full_rect=sample_full_rect,
                  k_sim=k_sim,
                  unsampled_rect=unsampled_rect,
                  q1_sci=q1_sci,
                  n_samp_sci=n_samp_sci,
                  k_com=1,
                  aggreg_obs=aggreg_obs,
                  Estimation_model=Estimation_model_i,
                  counter=counter,
                  sim=i)



# index_sci_i=index_sci_i,
# c_sci_x=c_sci_x,
# y_sci_i=y_sci_i,
# b=b,
# index_com_i=index_com_i,
# y_com_i=y_com_i,
# # b_com_i=b_com_i,
# c_com_x=c_com_x,
# aggreg_obs=aggreg_obs,
# k=k,
# sequencesdepeche=sequencesdepeche,
# zonespersequence=zonespersequence,
# taillezone=taillezone,
# 
List_param <- list(data.info = simu.est.data,
                   Opt_par = Opt,
                   SD = SD,
                   Report = Report)

if(counter == 1){
  
  Results <- data.frame(counter=counter,
                        sim=i,
                        Estimation_model=Estimation_model_i,
                        aggreg_obs=aggreg_obs,
                        n_seq_com=n_seq_com,
                        n_ping_per_seq=n_ping_per_seq,
                        n_samp_sci=n_samp_sci,
                        N_true=sum(S_x),
                        N_est=sum(Report$S_x),
                        MSPE=sum((S_x-Report$S_x)^2)/n_cells,
                        beta1_true=beta1,
                        beta1_est=Report$beta_j[2])
  Results$fixed <- NA
  Results$random <- NA
  
}else{
  
  Results[counter,"counter"] <- counter
  Results[counter,"sim"] <- i
  Results[counter,"Estimation_model"] <- Estimation_model_i
  Results[counter,"aggreg_obs"] <- aggreg_obs
  Results[counter,"n_seq_com"] <- n_seq_com
  Results[counter,"n_ping_per_seq"] <- n_ping_per_seq
  Results[counter,"n_samp_sci"] <- n_samp_sci
  Results[counter,"N_true"] <- sum(S_x)
  if(fit_IM_res$Converge==0) Results[counter,"N_est"] <- SD$unbiased$value["total_abundance"]
  Results[counter,"MSPE"] <- sum((S_x-Report$S_x)^2)/n_cells
  Results[counter,"beta1_true"] <- beta1
  Results[counter,"beta1_est"] <- Report$beta_j[2]
  
  if(Estimation_model_i==2){
    
    Results[counter,"fixed"] <- fixed
    Results[counter,"random"] <- random
    
  }
  
}

List_param <- list(data.info=simu.est.data,
                   Report=fit_IM_res$Report,
                   SD=fit_IM_res$SD,
                   converge=fit_IM_res$Converge)

if(! dir.exists(simu_file)) dir.create(simu_file,recursive = T)
save(data=Results,file=paste0(simu_file,"Results.RData"))
save(data=List_param,file=paste0(simu_file,"List_param_",counter,".RData"))


