#----------------------
## Plot from List_param
#----------------------

folder_name <- "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/com_x_sci_data_14_scientific_commercial_simple-2022-01-04_18_55_00_SimuRef"
folder_name <- "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/com_x_sci_data_14_scientific_commercial_simple-2022-01-05_15_15_00_SimuRef"

library(HistogramTools)

#-------------
## If parallel
#-------------
for(i in 1:10){
  load(paste0(folder_name,"/",i,"/Results",i,".RData"))
  if(i==1) Results_2 <- Results
  if(i!=1) Results_2 <- rbind(Results,Results_2)
}
library(stringr)
for(i in 1:10){
  
  # Load Results
  load(paste0(folder_name,"/",i,"/Results",i,".RData"))
  
  # Load List_param
  list_file_i <- list.files(paste0(folder_name,"/",i))
  list_file_i <- list_file_i[which(str_detect(list_file_i,"List_param"))]
  
  Results$aggreg_obs <- NA
  Results$sequencesdepeche <- NA
  Results$zonespersequence <- NA
  Results$n_samp_com <- NA
  Results$rho_S <- NA
  Results$beta1_true <- NA
  Results$beta1_est <- NA
  
  # Run over List_param
  for(counter in Results$counter){
    print(paste0("i: ",i," | counter: ",counter))
    load(paste0(folder_name,"/",i,"/List_param_",counter,".RData"))
    Results$aggreg_obs[which(Results$counter == counter)] <- List_param$data.info$aggreg_obs
    Results$sequencesdepeche[which(Results$counter == counter)] <- List_param$data.info$sequencesdepeche
    Results$zonespersequence[which(Results$counter == counter)] <- List_param$data.info$zonespersequence
    Results$n_samp_com[which(Results$counter == counter)] <- List_param$data.info$n_samp_com
    Results$rho_S[which(Results$counter == counter)] <- cor(List_param$data.info$Strue_x,List_param$Report$S_x)
    Results$beta1_true[which(Results$counter == counter)] <- List_param$data.info$beta[1]
    Results$beta1_est[which(Results$counter == counter)] <- List_param$Report$beta_j[2]
  }
  
  # Build full Results
  if(i==1) Results_2 <- Results
  if(i!=1) Results_2 <- rbind(Results,Results_2)
  
}

#----------------
## If no parallel
#----------------
load(paste0(folder_name,"/Results.RData"))
Results_2 <- Results

# Load List_param
list_file_i <- list.files(paste0(folder_name,"/"))
list_file_i <- list_file_i[which(str_detect(list_file_i,"List_param"))]
  
Results$aggreg_obs <- NA
Results$sequencesdepeche <- NA
Results$zonespersequence <- NA
Results$n_samp_com <- NA
Results$rho_S <- NA
Results$beta1_true <- NA
Results$beta1_est <- NA
Results$alpha <- NA
Results$beta <- NA
Results$gamma <- NA
Results$mspe <- NA
Results$mspe <- NA
Results$sigma_com_true <- NA
Results$sigma_com_est <- NA
Results$q1_com_est <- NA
Results$q1_com_true <- NA
Results$intercept_est <- NA
Results$intercept_true <- NA

# Run over List_param
for(counter in Results$counter){
  
  print(paste0("counter: ",counter))
  load(paste0(folder_name,"/List_param_",counter,".RData"))
  
  ## Config simu
  Results$aggreg_obs[which(Results$counter == counter)] <- List_param$data.info$aggreg_obs
  Results$aggreg_obs[which(Results$counter == counter)] <- List_param$data.info$aggreg_obs
  Results$sequencesdepeche[which(Results$counter == counter)] <- List_param$data.info$sequencesdepeche
  Results$zonespersequence[which(Results$counter == counter)] <- List_param$data.info$zonespersequence
  Results$n_samp_com[which(Results$counter == counter)] <- List_param$data.info$n_samp_com
  
  ## Coefficient de correlation
  Results$rho_S[which(Results$counter == counter)] <- cor(List_param$data.info$Strue_x,List_param$Report$S_x)

  ## MSPE
  Results$mspe[which(Results$counter == counter)] <- sum((log(List_param$data.info$Strue_x) - log(List_param$Report$S_x))^2) / n_cells
  
  ## Observation parameter
  Results$sigma_com_true[which(Results$counter == counter)] <- exp(List_param$data.info$logSigma_com)
  if(List_param$Opt_par$convergence == 0) Results$sigma_com_est[which(Results$counter == counter)] <- exp(List_param$SD$par.fixed["logSigma_com"])
  Results$q1_com_true[which(Results$counter == counter)] <- List_param$data.info$q1_com
  if(List_param$Opt_par$convergence == 0) Results$q1_com_est[which(Results$counter == counter)] <- List_param$SD$par.fixed["q1_com"]
  
  ## Species-habitat relationship
  Results$beta1_true[which(Results$counter == counter)] <- List_param$data.info$beta[1]
  Results$beta1_est[which(Results$counter == counter)] <- List_param$Report$beta_j[2]
  
  ## Intercept
  Results$intercept_true[which(Results$counter == counter)] <- List_param$data.info$beta0
  Results$intercept_est[which(Results$counter == counter)] <- List_param$Report$beta_j[1]
  
  
  # SPAEF
  S_sim <- List_param$data.info$Strue_x
  S_est <- List_param$Report$S_x
  alpha <- cor(S_sim,S_est)
  beta <- (sd(S_sim)/mean(S_sim)) / (sd(S_est)/mean(S_est))
  
  max_S <- max(c(S_sim,S_est))
  seq_S_levels <- seq(from = 0, to = max_S, length.out = 15)
  hist_S_sim <- hist(S_sim,breaks = seq_S_levels,plot = F)
  hist_S_est <- hist(S_est,breaks = seq_S_levels,plot = F)
  gamma <- intersect.dist(hist_S_sim,hist_S_est)
  
  Results$alpha[which(Results$counter == counter)] <- alpha
  Results$beta[which(Results$counter == counter)] <- beta
  Results$gamma[which(Results$counter == counter)] <- gamma
  
}


#-----------------------------------------------------------------------------------------------------------

library(doBy)
library(gt)
Results_conv <- Results
Results_conv$one <- 1
Results_conv <- Results_conv %>% 
  filter(b_true == 0)
summaryBy(Convergence+one~
            sequencesdepeche+
            # zonespersequence+
            aggreg_obs+
            # b_true+
            n_samp_com+
            reallocation,
          data=Results_conv,
          FUN=sum) %>%
  dplyr::select(n_samp_com,
                sequencesdepeche,
                # zonespersequence,
                reallocation,
                aggreg_obs,
                Convergence.sum,
                one.sum) %>%
  dplyr::rename(fishing_sequence = sequencesdepeche,
                data_realloc = reallocation,
                lkl_level = aggreg_obs) %>%
  mutate(lkl_level = ifelse(lkl_level == T,"Dj","Yi"),
         data_realloc = ifelse(data_realloc == 1,"Yes","No"),
         perc_convergence = round((1 - Convergence.sum / one.sum)*100,digits = 3)) %>%
  dplyr::select(-Convergence.sum,-one.sum) %>%
  gt() %>%
  tab_options(
    summary_row.background.color = "#ACEACE80",
    grand_summary_row.background.color = "#990000",
    row_group.background.color = "#FFEFDB80",
    heading.background.color = "#EFFBFC",
    column_labels.background.color = "#EFFBFC",
    stub.background.color = "#EFFBFC",
    table.font.color = "#323232",
    table_body.hlines.color = "#989898",
    table_body.border.top.color = "#989898",
    heading.border.bottom.color = "#989898",
    row_group.border.top.color = "#989898",
    row_group.border.bottom.style = "none",
    stub.border.style = "dashed",
    stub.border.color = "#989898",
    stub.border.width = "1px",
    summary_row.border.color = "#989898",
    table.width = "80%",
    column_labels.font.weight = "bold"
  )


Results <- Results %>%
  filter(b_true == 0 & Convergence == 0)
Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
Results[,"RelBias_beta"]=(Results[,"beta1_est"]-Results[,"beta1_true"])/Results[,"beta1_true"]
Results[,"Bias_b"]=(Results[,"b_est"]-Results[,"b_true"]) / ifelse(Results[,"b_true"] != 0,Results[,"b_true"],1)
Results[,"Bias_sigma_com"]=(Results[,"sigma_com_est"]-Results[,"sigma_com_true"]) / Results[,"sigma_com_true"]
Results[,"Bias_q1_com"]=(Results[,"q1_com_est"]-Results[,"q1_com_true"]) / Results[,"q1_com_true"]
Results[,"Bias_intercept"]=(Results[,"intercept_est"]-Results[,"intercept_true"]) / Results[,"intercept_true"]

Results <- Results %>%
  mutate(relloc_aggreg = paste0("realloc.sim:",reallocation,"_realloc.est:",ifelse(aggreg_obs==T,1,0)))

Results$zonespersequence <- as.factor(Results$zonespersequence)
Results$b_true <- as.factor(Results$b_true)

#----------------------------------------------------------------------
df <- Results %>%
  filter(n_samp_com == 1000)
ggplot()+
  geom_boxplot(data = df,
               aes(x = relloc_aggreg,
                   y = Bias_q1_com,
                   fill = relloc_aggreg))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)

ggplot()+
  geom_boxplot(data = df,
               aes(x = relloc_aggreg,
                   y = Bias_sigma_com,
                   fill = relloc_aggreg))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)

ggplot()+
  geom_boxplot(data = df,
               aes(x = relloc_aggreg,
                   y = Bias_intercept,
                   fill = relloc_aggreg))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)

#----------------------------------------------------------------------

RelBias_N_plot <- ggplot()+
  geom_boxplot(data = Results,
               aes(x = zonespersequence,
                   y = RelBias_N,
                   fill = relloc_aggreg))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  ylab("Relative bias of biomass")+
  facet_wrap(.~factor(n_samp_com))+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)+
  ylim(-1,5)

RelBias_beta_plot <- ggplot()+
  geom_boxplot(data = Results,
               aes(x = zonespersequence,
                   y = RelBias_beta,
                   fill = relloc_aggreg))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  ylab("Relative bias of beta1")+
  facet_wrap(.~factor(n_samp_com))+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)


MSPE_S_plot <- ggplot()+
  geom_boxplot(data = Results,
               aes(x = zonespersequence,
                   y = log(mspe),
                   fill = relloc_aggreg))+
  theme_bw()+
  # ylab("MSPE")+
  facet_wrap(.~factor(n_samp_com))+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)



Bias_SPAEF_plot <- ggplot()+
  geom_boxplot(data = Results,
               aes(x = zonespersequence,
                   y = 1 - sqrt( (alpha-1)^2 + (beta-1)^2 + (gamma-1)^2 ),
                   fill = relloc_aggreg))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  ylab("SPAEF")+
  facet_wrap(.~factor(n_samp_com))+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)+
  ylim(0,1)


library(cowplot)

plot.full <- plot_grid(
  RelBias_N_plot,
  RelBias_beta_plot,
  MSPE_S_plot,ncol = 1)


legend <- ggpubr::as_ggplot(cowplot::get_legend(RelBias_N_plot+theme(legend.position = "bottom", legend.title = element_blank())))

plot.full <- plot_grid(plot.full,legend,ncol=1,rel_heights = c(1,0.1))

ggsave(filename = "plot_full.png",width = 9,height = 9)
