
library(tidyverse);library(statforbiology):library(viridis)

FULL_Out_NAT <- data.frame()
for (repli in 1:length(out)) {
  FULL_Out_NAT <- bind_rows(FULL_Out_NAT,out[[repli]]$Out_NAT)
}

#4. Plot and analysis ####
#FULL_Out_NAT %>%
fake_output_Total %>% 
  group_by(Scenario,Type_NATS,n_sites)  %>% 
  summarise(mean_edge_dens=mean(edge_dens),mean_mean_grStre=mean(mean_grStre)) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=mean_edge_dens,colour=Type_NATS))+
  geom_line(aes(x=n_sites,y=mean_edge_dens,colour=Type_NATS,linetype=as.factor(Type_NATS)))+
  theme_classic()+facet_wrap(Scenario~.,scales="free")

scenario <- unique(FULL_Out_NAT$Scenario)
type <- unique(FULL_Out_NAT$Type_NATS)
output_slope <- data.frame()
for (repli in 1:4) {
  Repli_full_ouput_temp <- FULL_Out_NAT %>% filter(Replicates==repli)  
  for (sceni in 1:length(scenario)) {
    Sceni_full_ouput_temp <-  Repli_full_ouput_temp %>% filter(Scenario==scenario[sceni])
    for (ind_type in 1:length(type)) {
      # We create a "temporary" file filtered according to the TypeNat selected. 
      full_ouput_temp <- Sceni_full_ouput_temp %>% filter(Type_NATS==type[ind_type])
      random <- Sceni_full_ouput_temp %>% filter(Type_NATS=="Rand") %>% 
        group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      
      model_Random <- lm(mean_grStre~(as.numeric(n_sites)), data=random)
      
      Model_nats <- predict(model_Random)
      
      #lm <- summary(lm(mean_grStre~(as.numeric(n_sites)), data=full_output_temp_temp))
      #slope <- lm$coefficients[2,1]
      full_output_temp_temp <- full_ouput_temp %>%
        group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      slope <- abs(Model_nats-(full_output_temp_temp$mean_grStre))
      #fltrem cada any
      if(sceni==1){
        full_output_temp_temp <- full_ouput_temp %>%
          group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.)) %>% 
          filter(n_sites<30) 
      }else{
        full_output_temp_temp <- full_ouput_temp %>%
          group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      }
      
      Y <-  full_output_temp_temp$edge_dens
      X <- full_output_temp_temp$n_sites
      control1 <- nls.control(maxiter= 1000,tol=1e-02, warnOnly=TRUE)
      model_S <- nls(Y~NLS.asymReg(X, init, m, plateau),control=control1)
      output_slope2 <- data.frame("Replicates"=repli,
                                  "Scenario"=scenario[sceni],
                                  "Type_NATS"=type[ind_type],
                                  #"site"=full_output_temp_temp$n_sites,
                                  "Slope"=(slope),
                                  "m"=coefficients(model_S)[2],
                                  "plat"=coefficients(model_S)[3])
      
      output_slope <- rbind(output_slope,output_slope2)
    }#Type Nats
  }# sceni
}# Repli


output_slope%>% 
  group_by(Replicates,Scenario,Type_NATS) %>% summarise(m=mean(m),plat=mean(plat)) %>% 
  ggplot(aes(x = Type_NATS, y=m))+
  geom_jitter(shape=21,aes(fill=as.factor(Replicates)))+
  geom_boxplot(aes(colour=Type_NATS),alpha=0.5)+
  scale_fill_viridis(discrete = T)+
  facet_wrap(.~Scenario)+labs(fill="Replicates")+
  theme_classic()

output_slope%>% 
  ggplot() +
  aes(x = Type_NATS, y =Slope, color = Type_NATS) +
  geom_boxplot()+
  facet_wrap(.~Scenario,scales="free")+
  theme_classic()




