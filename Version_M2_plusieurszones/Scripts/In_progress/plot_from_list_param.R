#----------------------
## Plot from List_param
#----------------------

folder_name <- "C:/Users/balglave/Desktop/com_x_sci_data_14_scientific_commercial_simple-2021-09-19_15_07_00_SimuTest"

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

library(doBy)
library(gt)
Results_conv <- Results_2
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


Results <- Results_2 %>%
  filter(b_true == 0 & Convergence == 0)
Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
Results[,"RelBias_beta"]=(Results[,"beta1_est"]-Results[,"beta1_true"])/Results[,"beta1_true"]
Results[,"Bias_b"]=(Results[,"b_est"]-Results[,"b_true"]) / ifelse(Results[,"b_true"] != 0,Results[,"b_true"],1)

Results <- Results %>%
  mutate(relloc_aggreg = paste0("realloc.sim:",reallocation,"_realloc.est:",ifelse(aggreg_obs==T,1,0)))

Results$zonespersequence <- as.factor(Results$zonespersequence)
Results$b_true <- as.factor(Results$b_true)

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
        aspect.ratio = 1)

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

Bias_b_plot <- ggplot()+
  geom_boxplot(data = Results,
               aes(x = zonespersequence,
                   y = Bias_b,
                   fill = relloc_aggreg))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  ylab("Bias of b")+
  facet_wrap(.~factor(n_samp_com))+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)

MSPE_S_plot <- ggplot()+
  geom_boxplot(data = Results,
               aes(x = zonespersequence,
                   y = log(MSPE_S),
                   fill = relloc_aggreg))+
  theme_bw()+
  # ylab("MSPE")+
  facet_wrap(.~factor(n_samp_com))+
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)


library(cowplot)

plot.full <- plot_grid(
  RelBias_N_plot,
  RelBias_beta_plot,
  MSPE_S_plot,ncol = 1)


legend <- ggpubr::as_ggplot(cowplot::get_legend(RelBias_N_plot+theme(legend.position = "bottom", legend.title = element_blank())))

plot.full <- plot_grid(plot.full,legend,ncol=1,rel_heights = c(1,0.1))

ggsave(filename = "plot_full.png",width = 9,height = 9)
