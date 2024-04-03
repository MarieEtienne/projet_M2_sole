
consistency_df <- Results %>% 
  mutate(Model = ifelse(aggreg_obs==T,"Declaration model","Reallocated model"))
consistency_df$Model <- factor(consistency_df$Model,levels = c("Reallocated model","Declaration model"))
consistency_plot <- ggplot(consistency_df)+
  geom_boxplot(aes(x=Model,y=`p-val-fixed`,fill=Model))+
  theme_classic()+
  theme(aspect.ratio = 2.5/3,
        legend.position = "bottom",
        legend.title = element_blank())+
  xlab("")+ylab("p-value")+
  scale_fill_manual(breaks = c("Reallocated model","Declaration model"), 
                    values=c("#F8766D", "#00BA38"))
  

ggsave(consistency_plot,filename="consistency_plot.png",width = 4.5,height = 4.5)
