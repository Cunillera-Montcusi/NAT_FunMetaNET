
library(tidyverse);library(viridis)
library(drc);library(nlme);library(statforbiology)

load("ClusterNATs/Sim_ClusterNATs.RData")

Final_O_put$Out_NAT %>% 
  group_by(Year,Type_NATS,n_sites)  %>% 
  summarise(mean_edge_dens=mean(edge_dens),mean_mean_grStre=mean(mean_grStre)) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=mean_edge_dens,colour=Type_NATS))+
  geom_line(aes(x=n_sites,y=mean_edge_dens,colour=Type_NATS,linetype=as.factor(Type_NATS)))+
  theme_classic()+facet_wrap(Year~.,scales="free")


year <- unique(Final_O_put$Out_NAT$Year)
type <- unique(Final_O_put$Out_NAT$Type_NATS)
iter <- unique(Final_O_put$Out_NAT$iter)
output_slope <- data.frame()

for (sceni in 1:length(year)) {
    Sceni_full_ouput_temp <-  Final_O_put$Out_NAT %>% filter(Year==year[sceni])
    for (ind_type in 1:length(type)) {
      for (itera in 1:length(iter)) {
      # We create a "temporary" file filtered according to the TypeNat selected. 
      full_ouput_temp <- Sceni_full_ouput_temp %>% filter(Type_NATS==type[ind_type],iter==iter[itera]) # En cas de voler filtrar per iteració
      random <- Sceni_full_ouput_temp %>% filter(Type_NATS=="Random",iter==iter[itera]) %>% 
        group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      
      model_Random <- lm(mean_grStre~(as.numeric(n_sites)), data=random)
      Model_nats <- predict(model_Random)
      
      #lm <- summary(lm(mean_grStre~(as.numeric(n_sites)), data=full_output_temp_temp))
      #slope <- lm$coefficients[2,1]
      full_output_temp_temp <- full_ouput_temp %>%group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      slope <- abs(Model_nats-(full_output_temp_temp$mean_grStre))
      
      # Abans feiem la mitjana per tot però ara el model ja considera el patró mitjà 
      #if(sceni==1){
      #  full_output_temp_temp <- full_ouput_temp %>%
      #    group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.)) %>% 
      #    filter(n_sites<30) 
      #}else{
      #  full_output_temp_temp <- full_ouput_temp %>%
      #    group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      #}
      full_output_temp_temp <- full_ouput_temp
      Y <-  full_output_temp_temp$edge_dens
      X <- full_output_temp_temp$n_sites
      control1 <- nls.control(maxiter= 1000,tol=1e-02, warnOnly=TRUE)
      model_S <- nls(Y~NLS.asymReg(X, init, m, plateau),control=control1)
      output_slope2 <- data.frame("Year"=year[sceni],
                                  "Type_NATS"=type[ind_type],
                                  "iter"=iter[itera],
                                  #"site"=full_output_temp_temp$n_sites,
                                  "GrStr_Obs_vs_Rand"=(slope),
                                  "ED_Curve_accel"=coefficients(model_S)[2],
                                  "ED_Plateau"=coefficients(model_S)[3])
      
      output_slope <- rbind(output_slope,output_slope2)
    }#itera
  }#Type Nats
}# sceni

gridExtra::grid.arrange(

output_slope%>% 
  group_by(Year,Type_NATS) %>% 
  summarise(Std_Dev=sd(ED_Curve_accel),ED_Curve_accel=median(ED_Curve_accel)) %>% 
  ggplot(aes(x = Year, y=ED_Curve_accel))+
  geom_errorbar(aes(x = Year, y=ED_Curve_accel,ymin=ED_Curve_accel-Std_Dev ,ymax=ED_Curve_accel+Std_Dev,
                    colour=as.factor(Type_NATS)),alpha=0.5)+
  geom_point(shape=21,aes(fill=as.factor(Type_NATS)))+
  geom_line(aes(colour=Type_NATS))+
  scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
  labs(fill="Type NATS")+ guides(colour="none")+
  #facet_wrap(.~Type_NATS)+labs(fill="Year")+
  theme_classic(),

output_slope%>% 
  group_by(Year,Type_NATS) %>% 
  summarise(Std_Dev=sd(ED_Curve_accel),ED_Plateau=mean(ED_Plateau)) %>% 
  ggplot(aes(x = Year, y=ED_Plateau))+
  geom_point(shape=21,aes(fill=as.factor(Type_NATS)))+
  geom_errorbar(aes(x = Year, y=ED_Plateau,ymin=ED_Plateau-Std_Dev ,ymax=ED_Plateau+Std_Dev,
                    colour=as.factor(Type_NATS)),alpha=0.5)+
  geom_line(aes(colour=Type_NATS))+
  scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
  labs(fill="Type NATS")+ guides(colour="none")+
  #facet_wrap(.~Type_NATS)+labs(fill="Year")+
  scale_y_continuous(limits=c(0.25,0.75))+
  theme_classic(),

output_slope%>% 
  group_by(Year,Type_NATS) %>% 
  summarise(Std_Dev=sd(GrStr_Obs_vs_Rand),GrStr_Obs_vs_Rand=median(GrStr_Obs_vs_Rand)) %>% 
  ggplot(aes(x = Year, y=GrStr_Obs_vs_Rand ))+
  geom_point(shape=21,aes(fill=as.factor(Type_NATS)))+
  geom_errorbar(aes(x = Year, y=GrStr_Obs_vs_Rand ,ymin=GrStr_Obs_vs_Rand -Std_Dev ,ymax=GrStr_Obs_vs_Rand +Std_Dev,
                    colour=as.factor(Type_NATS)),alpha=0.5)+
  geom_line(aes(colour=Type_NATS))+
  scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
  labs(fill="Type NATS")+ guides(colour="none")+
  #facet_wrap(.~Type_NATS)+labs(fill="Year")+
  theme_classic(),

nrow=3)


