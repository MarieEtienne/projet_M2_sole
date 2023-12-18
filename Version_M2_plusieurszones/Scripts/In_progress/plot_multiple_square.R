#----------------------------------
## Load results for multiple square
#----------------------------------
library(HistogramTools)
library(gt)


folder_name <- "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/severa_rect-2022-01-11_17_16_00_SimuTest"
folder_name <- "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/severa_rect-2022-01-12_13_03_00_SimuUnsampledRect"
folder_name <- "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/severa_rect-2022-01-13_11_14_00_SimuFullArea"

folder_c <- c("C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/severa_rect-2022-01-12_13_03_00_SimuUnsampledRect",
              "C:/R_projects/projet_M2_sole/Version_M2_plusieurszones/results/severa_rect-2022-01-13_11_14_00_SimuFullArea")


for(folder_i in folder_c){
  
  load(paste0(folder_i,"/Results.RData"))
  
  # Load List_param
  list_file_i <- list.files(paste0(folder_i,"/"))
  list_file_i <- list_file_i[which(str_detect(list_file_i,"List_param"))]
  
  Results$rho_S <- NA
  Results$MAPE <- NA
  Results$mspe <- NA
  Results$alpha <- NA
  Results$beta <- NA
  Results$gamma <- NA
  Results$converge <- NA
  Results$N_est.2 <- NA
  Results$N_true.2 <- NA
  Results$Range <- NA
  Results$SD_Range <- NA
  Results$MargeSD <- NA
  Results$SD_MargeSD <- NA
  Results$sigma_com_true <- NA
  Results$sigma_com_est <- NA
  Results$q1_com_est <- NA
  Results$q1_com_true <- NA
  Results$intercept_est <- NA
  Results$intercept_true <- NA
  
  
  if(str_detect(folder_i,"UnsampledRect")) simu_type <- "Unsampled Rectangles"
  if(str_detect(folder_i,"FullArea")) simu_type <- "All area sampled"
  
  # Run over List_param
  for(counter in Results$counter){
    
    print(paste0("simu_type: ", simu_type," | counter: ",counter))
    load(paste0(folder_i,"/List_param_",counter,".RData"))
    Results$aggreg_obs[which(Results$counter == counter)] <- List_param$data.info$aggreg_obs
    Results$converge[which(Results$counter == counter)] <- List_param$converge
    Results$rho_S[which(Results$counter == counter)] <- cor(log(List_param$data.info$S_x),log(List_param$Report$S_x))
    Results$MAPE[which(Results$counter == counter)] <- sum(abs(log(List_param$data.info$S_x) - log(List_param$Report$S_x))) / n_cells
    Results$mspe[which(Results$counter == counter)] <- sum((log(List_param$data.info$S_x) - log(List_param$Report$S_x))^2) / n_cells
    Results$N_est.2[which(Results$counter == counter)] <- sum(List_param$Report$S_x / exp(List_param$Report$beta_j[1]))
    Results$N_true.2[which(Results$counter == counter)] <- sum(List_param$data.info$S_x / exp(List_param$data.info$intercept))
    
    if(List_param$converge == 0) Results$MargeSD[which(Results$counter == counter)] <- 1 / sqrt(4*pi) / exp(List_param$SD$par.fixed["logtau"]) / exp(List_param$SD$par.fixed["logkappa"])
    # if(List_param$converge == 0) Results$SD_MargeSD[which(Results$counter == counter)] <- List_param$SD$sd[which(names(List_param$SD$value) %in% "MargSD")]
    
    if(List_param$converge == 0) Results$Range[which(Results$counter == counter)] <- sqrt(8)/exp(List_param$SD$par.fixed["logkappa"])
    # if(List_param$converge == 0) Results$SD_Range[which(Results$counter == counter)] <- List_param$SD$sd[which(names(List_param$SD$value) %in% "Range")]
    
    ## Observation parameter
    Results$sigma_com_true[which(Results$counter == counter)] <- List_param$data.info$Sigma_com
    if(List_param$converge == 0) Results$sigma_com_est[which(Results$counter == counter)] <- exp(List_param$SD$par.fixed["logSigma_com"])
    Results$q1_com_true[which(Results$counter == counter)] <- List_param$data.info$q1_com
    if(List_param$converge == 0) Results$q1_com_est[which(Results$counter == counter)] <- List_param$SD$par.fixed["q1_com"]
    
    ## Species-habitat relationship
    Results$beta1_true[which(Results$counter == counter)] <- List_param$data.info$beta[1]
    Results$beta1_est[which(Results$counter == counter)] <- List_param$Report$beta_j[2]
    
    ## Intercept
    Results$intercept_true[which(Results$counter == counter)] <- List_param$data.info$intercept
    Results$intercept_est[which(Results$counter == counter)] <- List_param$Report$beta_j[1]
    
    
    # SPAEF
    S_sim <- log(List_param$data.info$S_x)
    S_est <- log(List_param$Report$S_x)
    alpha <- cor(S_sim,S_est)
    beta <- (sd(S_sim)/mean(S_sim)) / (sd(S_est)/mean(S_est))
    max_S <- max(c(S_sim,S_est))
    min_S <- min(c(S_sim,S_est))
    
    seq_S_levels <- seq(from = min_S, to = max_S, length.out = 15)
    hist_S_sim <- hist(S_sim,breaks = seq_S_levels,plot = F)
    hist_S_est <- hist(S_est,breaks = seq_S_levels,plot = F)
    gamma <- intersect.dist(hist_S_sim,hist_S_est)
    
    Results$alpha[which(Results$counter == counter)] <- alpha
    Results$beta[which(Results$counter == counter)] <- beta
    Results$gamma[which(Results$counter == counter)] <- gamma
    
  }
  
  Results <- Results %>% mutate(simu_type = simu_type)
  
  if(folder_i == folder_c[1]){
    Results_2 <- Results
  }else{
    Results_2 <- rbind(Results_2,Results)
  }
  
}

Results_2$Model <- NA
Results_2$Model[which(Results_2$Estimation_model==1)] <- "Integrated model"
Results_2$Model[which(Results_2$Estimation_model==2)] <- "Scientific model"
Results_2$Model[which(Results_2$Estimation_model==3)] <- "Commercial model"


#---------------------
## Convergence results
#---------------------
Results_conv <- Results_2
Results_conv$one <- 1
summaryBy(converge+one~
            aggreg_obs+
            Model,
          data=Results_conv,
          FUN=sum) %>%
  dplyr::rename(lkl_level = aggreg_obs) %>%
  mutate(lkl_level = ifelse(lkl_level == T,"Dj","Yi"),
         perc_convergence = round((1 - converge.sum / one.sum)*100,digits = 3)) %>%
  filter(! (lkl_level == "Dj" & Model == "Scientific model") ) %>%
  dplyr::select(-converge.sum,-one.sum) %>%
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

#----------------------------
## Boxplot of all simulations 
#----------------------------
Results_plot <- Results_2 %>%
  filter(converge == 0)
Results_plot[,"RelBias_N"]=(Results_plot[,"N_est"]-Results_plot[,"N_true"])/Results_plot[,"N_true"]
Results_plot[,"RelBias_N.2"]=(Results_plot[,"N_est.2"]-Results_plot[,"N_true.2"])/Results_plot[,"N_true.2"]
Results_plot[,"RelBias_beta"]=(Results_plot[,"beta1_est"]-Results_plot[,"beta1_true"])/Results_plot[,"beta1_true"]
Results_plot[,"Bias_sigma_com"]=(Results_plot[,"sigma_com_est"]-Results_plot[,"sigma_com_true"]) / Results_plot[,"sigma_com_true"]
Results_plot[,"Bias_q1_com"]=(Results_plot[,"q1_com_est"]-Results_plot[,"q1_com_true"]) / Results_plot[,"q1_com_true"]
Results_plot[,"Bias_intercept"]=(Results_plot[,"intercept_est"]-Results_plot[,"intercept_true"]) / Results_plot[,"intercept_true"]

Results_plot$obs <- NA
Results_plot$obs[which(Results_plot$aggreg_obs == T)] <- "Dj"
Results_plot$obs[which(Results_plot$aggreg_obs == F)] <- "Yi"
Results_plot$obs[which(Results_plot$Estimation_model == 2)] <- "Scientific model"
Results_plot$obs <- factor(Results_plot$obs,levels = c("Yi","Dj","Scientific model"))

#-----------------------------------------------------------
ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = Range,
                   fill = Model))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)+
  geom_hline(yintercept = 0.6)+
  ylim(0,2)

ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = MargeSD,
                   fill = Model))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)+
  geom_hline(yintercept = 1)+
  ylim(0,2)


ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = Bias_intercept,
                   fill = Model))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)
  
ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = Bias_q1_com,
                   fill = Model))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)

ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = sigma_com_est,
                   fill = Model))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)

#-----------------------------------------------------------

# Results_plot_fullArea <- Results_plot
# Results_plot_UnsampArea <- Results_plot

beta_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = beta1_est,
                   fill = Model))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.title = element_blank(),
        # legend.position = "none",
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)+
  ylab(TeX("$\\beta_S$")) +
  xlab("")+
  geom_hline(yintercept = 2,col="darkgrey",linetype="dashed",size=1)

MSPE_S_plot <- ggplot()+
  geom_boxplot(data = Results_plot,
               aes(x = obs,
                   y = mspe,
                   fill = Model))+
  theme_bw()+
  theme(legend.title = element_blank(),
        # legend.position = "none",
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)+
  ylab("MSPE")+xlab("")

plot_grid(beta_plot,MSPE_S_plot,ncol=1)

ggplot(Results_plot)+
  # geom_point(aes(x = obs,
  #                y = Range + 1.96 * SD_Range,
  #                fill = Model),col="blue")+
  geom_boxplot(aes(x = obs,
                 y = Range + 1.96 * SD_Range,
                 fill = Model),col="blue")+
  # geom_point(aes(x = obs,
  #                y = Range,
  #                fill = Model))+
  theme_bw()+
  theme(legend.title = element_blank(),
        # legend.position = "none",
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)+
  ylab("")+xlab("")+
  ylim(0,2.5)+
  geom_hline(yintercept = 0.6,linetype="dashed")+
  ggtitle("Range")



ggplot(Results_plot)+
  
  # geom_point(aes(x = obs,
  #                y = MargeSD + 1.96 * SD_MargeSD,
  #                fill = Model),col="blue")+
  geom_boxplot(aes(x = obs,
                 y = MargeSD  + 1.96 * SD_MargeSD,
                 fill = Model),col="blue")+
  
  # geom_point(aes(x = obs,
  #                  y = MargeSD,
  #                  fill = Model))+
  theme_bw()+
  theme(legend.title = element_blank(),
        # legend.position = "none",
        aspect.ratio = 1)+
  facet_wrap(.~simu_type)+
  ylab("")+xlab("")+
  # ylim(0,2)+
  geom_hline(yintercept = 1,linetype="dashed")+
  ggtitle("MargeSD")

# MAPE_plot <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = obs,
#                    y = MAPE,
#                    fill = Model))+
#   theme_bw()+
#   theme(legend.title = element_blank(),
#         # legend.position = "none",
#         aspect.ratio = 1)+
#   facet_wrap(.~simu_type)
# 
# Corr_plot <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = obs,
#                    y = rho_S,
#                    fill = Model))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   theme(legend.title = element_blank(),
#         # legend.position = "none",
#         aspect.ratio = 1)+
#   facet_wrap(.~simu_type)
# 
# RelBias_N_plot <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = Model,
#                    y = RelBias_N.2,
#                    fill = obs))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   ylab("Relative bias of biomass")+
#   theme(legend.title = element_blank(),
#         # legend.position = "none",
#         aspect.ratio = 1)+
#   facet_wrap(.~simu_type)
# 
# consistency_fixed_plot <- ggplot()+
#   geom_boxplot(data = Results_2,
#                aes(x = aggreg_obs,
#                    y = fixed,
#                    fill = aggreg_obs))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   theme(legend.title = element_blank(),
#         # legend.position = "none",
#         aspect.ratio = 1)+
#   facet_wrap(.~simu_type)
# 
# consistency_random_plot <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = aggreg_obs,
#                    y = random,
#                    fill = aggreg_obs))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   ylab("SPAEF")+
#   theme(legend.title = element_blank(),
#         # legend.position = "none",
#         aspect.ratio = 1)
# 
# Bias_SPAEF_plot <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = Model,
#                    y = 1 - sqrt( (alpha-1)^2 + (beta-1)^2 + (gamma-1)^2 ),
#                    fill = aggreg_obs))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   ylab("SPAEF")+
#   theme(legend.title = element_blank(),
#         # legend.position = "none",
#         aspect.ratio = 1)




