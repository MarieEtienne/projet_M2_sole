###############
## Plot_Results
###############
#' @title Plot_Results()
#' 
#' @param Results summary dataframe of the simulation/estimation loop
#' @param b_set levels of preferential sampling
#' 
#' @return RelBias_N : Plot of the relative bias of abundance
#' @return Bias_b : Plot of the relative bias of b
#' @return MSPE_S : Plot of the mean squared prediction error
#' @return Converge_table : convergence table
#' 

Plot_Results <- function(Results,b_set){
  
  # Identify simulations/estimations that have failed
  Results[,"OnBoundary"]=rep(0,nrow(Results))
  Results[which(abs(Results[,"b_est"])>10),"OnBoundary"]=1
  Results[,"Converge"]=Results[,"OnBoundary"]+Results[,"Convergence"]
  
  # For figures, limit to results for which |b estimate| < 10
  Results[which(Results[,"Converge"]==2),"Converge"]=1
  
  # Convergence table
  Converge_table = summaryBy(Converge~b_true+Data_source,data=Results,FUN=sum)
  
  # Filter simulation which converged : on ne conserve dans les résultats que les
  #simulations ou l'algorithme a convergé
  WhichFailed <- which(Results[,"Convergence"]!=0)
  Results_cvg.failed <- Results[WhichFailed,]
  WhichDone = which(Results[,"Convergence"]==0)
  Results=Results[WhichDone,]
  
  # add relative bias for abundance
  Results[,"RelBias_N"]=(Results[,"N_est"]-Results[,"N_true"])/Results[,"N_true"]
  Results[,"Bias_b"]=(Results[,"b_est"]-Results[,"b_true"]) / ifelse(Results[,"b_true"] != 0,Results[,"b_true"],1)
  # plot estimates of 'b'
  Which_plot = which(Results[,"type_b"] == "est_b")
  
  # Handle data to have pretty graph (a bit dirty...)
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
  Results_2$b_est[which(is.na(Results_2[,"b_est"]))] <- 0
  Results_plot <- Results_2[which(Results_2[,"b_est"]<10),]
  
  # Plots
  Results_plot$Data_source <- factor(Results_plot$Data_source, levels = c("scientific_only", "commercial_only", "scientific_commercial"))
  
  #construction du boxplot pour le biais de l'abondance
  RelBias_N <- ggplot()+
    geom_boxplot(data = Results_plot,
                 aes(x = as.factor(b_true),
                     y = RelBias_N,
                     fill = Data_source))+
    geom_hline(yintercept = 0,linetype="dashed")+
    theme_bw()+
    theme(legend.position='none')+
    scale_fill_manual(breaks = c("scientific_only", "commercial_only", "scientific_commercial"), 
                      values=c("#56B4E9","#999999","#E69F00"))+
    ylab("Relative bias of abundance")+
    xlab(" ")
  
  #construction du graphique pour le biais de b 
  Bias_b <- ggplot()+
    geom_boxplot(data = Results_plot,
                 aes(x = as.factor(b_true),
                     y = Bias_b,
                     fill = Data_source))+
    geom_hline(yintercept = 0,linetype="dashed")+
    theme_bw()+
    theme(legend.position='none')+
    scale_fill_manual(breaks = c("commercial_only", "scientific_commercial"), 
                      values=c("#999999","#E69F00"))+
    ylab("Relative bias of b")+
    xlab(" ")
  
  #construction du graphique pour le MSPE
  MSPE_S <- ggplot()+
    geom_boxplot(data = Results_plot,
                 aes(x = as.factor(b_true),
                     y = MSPE_S,
                     fill = Data_source))+
    theme_bw()+
    theme(legend.position='none')+
    scale_fill_manual(breaks = c("scientific_only", "commercial_only", "scientific_commercial"), 
                      values=c("#56B4E9","#999999","#E69F00"))+
    ylab("MSPE")+
    xlab(" ")
  
  res <- list(RelBias_N=RelBias_N,
              Bias_b=Bias_b,
              MSPE_S=MSPE_S,
              Converge_table=Converge_table)
  
  return(res)
}





# save(file = paste0("images/Simu/",Version,"/RData/",simu_name,"_Alpha",Alpha,"_plot.RData"),save_plot_eval)

#################
# All other plots
#################

## Plot Strue_df vs PE_df and SD.S_df
# Strue_df <- Results_loop$Strue_df
# PE_df <- Results_loop$PE_df
# SD.S_df <- Results_loop$SD.S_df
# 
# 
# PE_df <- PE_df[,c(1:5)]
# colnames(PE_df)[1:625] <- paste0("cell_",c(1:n_cells))
# 
# SD.S_df <- SD.S_df[,-c(1:5)]
# colnames(SD.S_df)[1:625] <- paste0("cell_",c(1:n_cells))
# 

# Strue_df %>%
#   pivot_longer("cell_1":paste0("cell_",n_cells)) %>%
#   dplyr::rename(cell=name, Strue=value) %>%
#   dplyr::select(-metrics) -> Strue_df.2
# 
# PE_df %>%
#   pivot_longer("cell_1":paste0("cell_",n_cells)) %>%
#   dplyr::rename(cell=name, PE=value) %>%
#   dplyr::select(-metrics) -> PE_df.2
# 
# SD.S_df %>%
#   pivot_longer("cell_1":paste0("cell_",n_cells)) %>%
#   dplyr::rename(cell=name, SD.S=value) %>%
#   dplyr::select(-metrics) -> SD.S_df.2
# 
# Strue_df.2 <- data.table(na.omit(Strue_df.2))
# PE_df.2 <- data.table(na.omit(PE_df.2))
# SD.S_df.2 <- data.table(na.omit(SD.S_df.2))
# Strue_vs_PE <- merge(Strue_df.2, PE_df.2, by =c("counter","sim","b_true","Data_source","cell"))
# Strue_vs_SD.S <- merge(Strue_df.2, SD.S_df.2, by =c("counter","sim","b_true","Data_source","cell"))
# 
# Strue_vs_PE_plot <- ggplot(Strue_vs_PE)+
#   geom_point(aes(x=Strue,y=PE,col=Data_source))+
#   theme_bw()+
#   facet_wrap(.~b_true)
# 
# # x11();plot(Strue_vs_PE_plot)
# 
# Strue_vs_SD.S_plot <- ggplot(Strue_vs_SD.S)+
#   geom_point(aes(x=Strue,y=SD.S,col=Data_source))+
#   theme_bw()+
#   facet_wrap(.~b_true)
# 
# # to assess effect of commercial sample size on estimates and predictions
# library(RColorBrewer)
# n_color <- 3
# palette <- "RdBu"
# color_palette <- brewer.pal(n = n_color, name = palette)
# 
# Bias_N_fw <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = as.factor(n_samp_com),
#                    y = RelBias_N,
#                    fill = Data_source))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 90,size=10))+
#   # scale_fill_manual(" ",breaks = c("50", "150", "300"),
#   #                   values=c(color_palette[1],color_palette[2],color_palette[3]))+
#   ylab("Biais relatif de l'abondance")+
#   xlab("")+
#   facet_wrap(.~b_true)+ylim(-2,2)
# 
# Bias_b_fw <- ggplot()+
#   geom_boxplot(data = Results_plot[Results_plot[,"Data_source"]!="scientific_only",],
#                aes(x = as.factor(n_samp_com),
#                    y = Bias_b,
#                    fill = Data_source))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 90,size=10))+
#   # scale_fill_manual(" ",breaks = c("50", "150", "300"),
#   #                   values=c(color_palette[1],color_palette[2],color_palette[3]))+
#   ylab("Biais de l'échantillonnage préférentiel (b)")+
#   xlab(" ")+
#   facet_wrap(.~b_true)+ylim(-5,5)
# 
# MSPE_S_fw <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = as.factor(n_samp_com),
#                    y = MSPE_S,
#                    fill = Data_source))+
#   theme_bw()+
#   ylab("Carré moyen des erreurs de prédiction")+
#   xlab(" ")+
#   theme(axis.text.x=element_text(angle = 90,size=10))+
#   # scale_fill_manual(" ",breaks = c("50", "150", "300"),
#   #                   values=c(color_palette[1],color_palette[2],color_palette[3]))+
#   facet_wrap(.~b_true)+ylim(0,1500)
# 
# #plot bias as function of b, estimation method
# #first, rearrange relative bias in form for ggplot
# Bias.df=Results_plot
# Bias.df[,"B"]=Bias.df[,"b_true"]
# Bias.df[,"Bias"]=Bias.df[,"RelBias_N"]
# 
# 
# #plot proportion relative bias
# bias.plot.1 = ggplot(Bias.df,aes(factor(Data_source),Bias))+geom_boxplot()+facet_grid(~b_true) #,scales="free")
# bias.plot.1 = bias.plot.1 + theme(axis.text.x=element_text(angle = 90))
# bias.plot.1 = bias.plot.1 + labs(x = "Estimation model", y="Proportion relative bias")
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# 
# 
# #bias of b parameter
# BiasB.df = Bias.df[which(!is.na(Bias.df[,"Bias_b"])),]
# bias.plot = ggplot(BiasB.df[BiasB.df[,"Data_source"]!="scientific_only",],aes(factor(b_true),Bias_b))+geom_boxplot() +facet_grid(~Data_source)#,scales="free")
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + ylim(-1,1.5) + labs(x = "True b parameter", y="Absolute bias")
# bias.plot
# 
# ## produce plot of bias of species-habitat relationship parameters
# bias.plot = ggplot(Bias.df,aes(factor(b_true),BiasBetaj1))+geom_boxplot()+facet_grid(~Data_source) #,scales="free")
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot
# 
# bias.plot = ggplot(Bias.df,aes(factor(b_true),BiasBetaj2))+geom_boxplot()+facet_grid(~Data_source) #,scales="free")
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot
# 
# ## produce plot of bias of sampling intensity-covariates relationship parameters
# bias.plot = ggplot(Bias.df,aes(factor(b_true),BiasBetak1))+geom_boxplot()+facet_grid(~Data_source) #,scales="free")
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot
# 
# bias.plot = ggplot(Bias.df,aes(factor(b_true),BiasBetak2))+geom_boxplot()+facet_grid(~Data_source) #,scales="free")
# #bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
# #bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
# bias.plot=bias.plot + labs(x = "Estimation model", y="Absolute bias")
# bias.plot


#   # Bias Sigma_sci
# Sigma_sci <- ggplot()+
#   geom_boxplot(data = Results_plot,
#                aes(x = as.factor(b_true),
#                    y = Bias_Sigma_sci,
#                    fill = Data_source))+
#   scale_fill_manual(breaks = c("scientific_only", "commercial_only", "scientific_commercial"),
#                     values=c("#56B4E9","#999999","#E69F00"))+
#   theme_bw()+
#   theme(legend.position='none')+
#   ylab("Biais Sigma_sci")+
#   xlab(" ")
# 
# Results_plot %>%
#   mutate(RelBias_exp.b = ((exp(b_est)-exp(b_true))/exp(b_true))) -> Results_plot
# 
# RelBias_exp.b <- ggplot()+
#   geom_boxplot(data = Results_plot[which(Results_plot$Data_source != "scientific_only"),],
#                aes(x = as.factor(b_true),
#                    y = RelBias_exp.b,
#                    fill = Data_source))+
#   geom_hline(yintercept = 0,linetype="dashed")+
#   theme_bw()+
#   theme(legend.position='none')+
#   scale_fill_manual(breaks = c("commercial_only", "scientific_commercial"), 
#                     values=c("#999999","#E69F00"))+
#   ylab("Biais relatif b")+
#   xlab(" ")
#
# Bias_b.N <- ggplot()+
#   geom_point(data = Results_plot[Results_plot[,"Data_source"]!="scientific_only",], aes(x = Bias_b, y = RelBias_N, col = as.character(b_true)))+
#   facet_wrap(.~Data_source)
# # Bias table
# sd <- function(x)sqrt(var(x))
# Bias_table = summaryBy(Bias_b+RelBias_N~b_true+Data_source,data=Results_plot,FUN=c(median,mean,sd))


