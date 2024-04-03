
test <- tacsatEflalo[,colnames(tacsatEflalo_2)[which(str_detect(colnames(tacsatEflalo_2),"LE_KG"))]]
for(i in 1:ncol(test)){
  
  vec_i <- test[,i]
  
  if(i==1) df <- data.frame(sp=colnames(test)[i],prop=length(which(vec_i==0))/length(vec_i))
  if(i>1){
    
    toto <- data.frame(sp=colnames(test)[i],prop=length(which(vec_i==0))/length(vec_i))
    df <- rbind(toto,df)
      
  } 

}

df %>% arrange(desc(prop))
