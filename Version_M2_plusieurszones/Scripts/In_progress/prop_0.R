it <- 1
results_df <- data.frame(it=1,q1=q1,prop=0)
for(j in 1:100){
  
  for(q1 in c(-5,-4,-3,-1)){
    
    print(paste0(q1," | ",j))  
    source("Scripts/source/simulate_real.R")
    
    results_df[it,"it"] <- it
    results_df[it,"q1"] <- q1
    results_df[it,"prop"] <- length(which(sum_y_i == 0))/length(sum_y_i)
    it <- it + 1
  }
  
}

ggplot(results_df)+
  geom_boxplot(aes(x=factor(q1),y=prop))

