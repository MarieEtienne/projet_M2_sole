


# Results$n_samp_com
# Results$nsci_ncom <- 0
# Results$nsci_ncom[which(Results$n_samp_sci == 0 | Results$n_samp_com == 200)] <- "0_200"
# Results$nsci_ncom[which(Results$n_samp_sci == 50 | Results$n_samp_com == 150)] <- "50_150"
# Results$nsci_ncom[which(Results$n_samp_sci == 100 | Results$n_samp_com == 100)] <- "100_100"
# Results$nsci_ncom[which(Results$n_samp_sci == 150 | Results$n_samp_com == 50)] <- "150_50"
# Results$nsci_ncom[which(Results$n_samp_sci == 200 | Results$n_samp_com == 0)] <- "200_0"

# Results$nsci_ncom <- factor(Results$nsci_ncom,levels = c("200_0","150_50","100_100","50_150","0_200"))
Results$Data_source <- factor(Results$Data_source,levels = c("scientific_only","commercial_only","scientific_commercial"))

Results[,"OnBoundary"]=rep(0,nrow(Results))
Results[which(abs(Results[,"b_est"])>10),"OnBoundary"]=1
Results[,"Converge"]=Results[,"OnBoundary"]+Results[,"Convergence"]
#For figures, limit to results for which |b estimate| < 10
Results[which(Results[,"Converge"]==2),"Converge"]=1

Converge_table = summaryBy(Converge~sci_sampling+Data_source,data=Results,FUN=sum)

# length(Results$Convergence[which(Results$Convergence == 1 & Results$Data_source == "scientific_only")])/length(Results$Convergence[which(Results$Data_source == "scientific_only")])

# Filter simulation which converged
WhichFailed <- which(Results[,"Convergence"]!=0)
Results_cvg.failed <- Results[WhichFailed,]
WhichDone = which(Results[,"Convergence"]==0)
Results=Results[WhichDone,]

# add relative bias for abundance
Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
Results[,"Bias_b"]=Results[,"b_est"]-Results[,"b_true"]

# Handle data to have pretty graph 
Results %>%
  filter(Data_source == "scientific_only") -> Results_sc_only
Results_sc_only <- Results_sc_only[which(is.na(Results_sc_only$b_true)),]
Results %>%
  filter(Data_source != "scientific_only") -> Results_com

Results_sc_only$b_true <- b_set[1]
df_1 <- Results_sc_only
Results_sc_only <- do.call(rbind.data.frame, lapply(2:length(b_set),function(i){
  df_1$b_true <- b_set[i]
  Results_sc_only <- rbind(Results_sc_only,df_1)
}))

Results_2 <- rbind(Results_sc_only,Results_com)
# Results_sc_only %>% mutate(b_true = b_set[1]) -> Results_sc_only_0
# Results_sc_only %>% mutate(b_true = 1) -> Results_sc_only_0
# 
# Results <- rbind(Results_sc_only,Results_com)

Results_2$b_est[which(is.na(Results_2[,"b_est"]))] <- 0
Results_plot <- Results_2[which(Results_2[,"b_est"]<10),]


Bias_b.N <- ggplot()+
  geom_point(data = Results_plot, aes(x = Bias_b, y = RelBias_N, col = as.character(b_true)))+
  xlim(-2,2)+ylim(-1,1)+
  facet_wrap(.~nsci_ncom)

Bias_b.MSPE <- ggplot()+
  geom_point(data = Results_plot, aes(x = Bias_b, y = MSPE_S, col = as.character(b_true)))+
  xlim(-2,2)+ylim(0,2000)+
  facet_wrap(.~nsci_ncom + Data_source)

RelBias_N <- ggplot()+
  geom_boxplot(data = Results_plot, # [which(Results_plot$Data_source == "scientific_only"),]
               aes(x = as.factor(n_samp_com),
                   y = RelBias_N,
                   fill= Data_source))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_fill_manual(breaks = c("scientific_only", "commercial_only", "scientific_commercial"),
                    values=c("#56B4E9","#999999","#E69F00"))+
  ylab("Biais relatif de l'abondance")+
  xlab("")+ylim(-1,2.5)+
  facet_wrap(.~factor(b_true))

Bias_b <- ggplot()+
  geom_boxplot(data = Results_plot,  # [which(Results_plot$Data_source == "scientific_commercial"),]
               aes(x = as.factor(n_samp_com),
                   y = Bias_b,
                   fill = Data_source))+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_bw()+
  theme(legend.position='none')+
  scale_fill_manual(breaks = c("commercial_only", "scientific_commercial"), 
                    values=c("#999999","#E69F00"))+
  ylab("Biais de b")+
  xlab("")+
  facet_wrap(.~as.factor(b_true))

MSPE <- ggplot()+
  geom_boxplot(data = Results_plot,  # [which(Results_plot$Data_source == "scientific_only"),]
               aes(x = as.factor(n_samp_com),
                   y = MSPE_S,
                   fill = Data_source))+
  theme_bw()+
  theme(legend.position='none')+
  scale_fill_manual(breaks = c("scientific_only","commercial_only", "scientific_commercial"), 
                    values=c("#56B4E9","#999999","#E69F00"))+
  ylab("MSPE")+
  xlab("")+
  ylim(0,2500)+
  facet_wrap(.~as.factor(b_true))

grid.arrange(RelBias_N,Bias_b,MSPE)

# save(data = Results_plot,file = "C:/R_projects/phd-baptiste_alglave/Rmd/Paper_1/res/Results_plot_q1_effect_commercial_sample.RData")


ggplot()+
  geom_boxplot(data = Results_plot,  # [which(Results_plot$Data_source == "scientific_only"),]
               aes(x = as.factor(b_true),
                   y = RelBias_N,
                   fill = Data_source))+
  theme_bw()+
  theme(legend.position='none')+
  scale_fill_manual(breaks = c("scientific_only","commercial_only", "scientific_commercial"), 
                    values=c("#56B4E9","#999999","#E69F00"))+
  ylab("Biais relatif de l'abondance")+
  xlab("")+
  ylim(-1.5,1.5)+
  facet_wrap(.~as.factor(nsci_ncom))

ggplot()+
  geom_boxplot(data = Results_plot,  # [which(Results_plot$Data_source == "scientific_only"),]
               aes(x = as.factor(b_true),
                   y = Bias_b,
                   fill = Data_source))+
  theme_bw()+
  theme(legend.position='none')+
  scale_fill_manual(breaks = c("scientific_only","commercial_only", "scientific_commercial"), 
                    values=c("#56B4E9","#999999","#E69F00"))+
  ylab("Biais de b")+
  xlab("")+
  ylim(-5,5)+
  facet_wrap(.~as.factor(nsci_ncom))


# Bias Sigma_sci
Sigma_sci <- ggplot()+
  geom_boxplot(data = Results_plot[which(Results_plot$Data_source!="commercial_only"),],
               aes(x = nsci_ncom,
                   y = Bias_Sigma_sci,
                   fill = Data_source))+
  scale_fill_manual(breaks = c("scientific_only", "commercial_only", "scientific_commercial"),
                    values=c("#56B4E9","#999999","#E69F00"))+
  theme_bw()+
  theme(legend.position='none')+
  ylab("Biais de Sigma_sci")+
  xlab("")+
  facet_wrap(.~as.factor(b_true))


grid.arrange(RelBias_N,Bias_b,MSPE_S,nrow=3)

save_plot_eval <- list(RelBias_N=RelBias_N,
                       Bias_b=Bias_b,
                       MSPE_S=MSPE_S,
                       Sigma_sci=Sigma_sci)
save(file = paste0("images/Simu/",Version,"/RData/",simu_name,"_Alpha",Alpha,"_plot.RData"),save_plot_eval)
