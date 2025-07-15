
library(tidyverse);library(viridis)
library(drc);library(nlme);library(statforbiology)

load("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/12. FunMetaNet/NAT/NAT_FunMetaNET/FakeMetas/OUT_fake_NATs.RData")

FULL_Out_NAT <- data.frame()
for (repli in 1:length(out)) {
  FULL_Out_NAT <- bind_rows(FULL_Out_NAT,out[[repli]]$Out_NAT)
}

#4. Plot and analysis ####
FULL_Out_NAT %>%
  group_by(Scenario,Type_NATS,n_sites)  %>% 
  summarise(mean_edge_dens=mean(edge_dens),mean_mean_grStre=mean(mean_grStre)) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=mean_edge_dens,colour=Type_NATS))+
  geom_line(aes(x=n_sites,y=mean_edge_dens,colour=Type_NATS,linetype=as.factor(Type_NATS)))+
  theme_classic()+
  facet_wrap(Scenario~.,scales="free")

scenario <- unique(FULL_Out_NAT$Scenario)
type <- unique(FULL_Out_NAT$Type_NATS)
iter <- unique(FULL_Out_NAT$iter)
output_slope <- data.frame()
for (repli in 1:length(out)) {
  Repli_full_ouput_temp <- FULL_Out_NAT %>% filter(Replicates==repli)  
  for (sceni in 1:length(scenario)) {
    Sceni_full_ouput_temp <-  Repli_full_ouput_temp %>% filter(Scenario==scenario[sceni])
    for (ind_type in 1:length(type)) {
      #for (itera in 1:length(iter)) {
      # We create a "temporary" file filtered according to the TypeNat selected. 
      full_ouput_temp <- Sceni_full_ouput_temp %>% filter(Type_NATS==type[ind_type])#,iter==iter[itera]) # En cas de voler filtrar per iteració
      random <- Sceni_full_ouput_temp %>% filter(Type_NATS=="Rand") %>%#,iter==iter[itera]) %>% 
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
      #  full_output_temp_temp <- full_ouput_temp %>%group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      #}
      full_output_temp_temp <- full_ouput_temp
      Y <-  full_output_temp_temp$edge_dens
      X <- full_output_temp_temp$n_sites
      control1 <- nls.control(maxiter= 1000,tol=1e-02, warnOnly=TRUE)
      model_S <- nls(Y~NLS.asymReg(X, init, m, plateau),control=control1)
      output_slope2 <- data.frame("Replicates"=repli,
                                  "Scenario"=scenario[sceni],
                                  "Type_NATS"=type[ind_type],
                                  #"iter"=iter[itera],
                                  #"site"=full_output_temp_temp$n_sites,
                                  "GrStr_Obs_vs_Rand"=(slope),
                                  "ED_Curve_accel"=coefficients(model_S)[2],
                                  "ED_Plateau"=coefficients(model_S)[3])
      
      output_slope <- rbind(output_slope,output_slope2)
     # }#itera
    }#Type Nats
  }# sceni
}# Repli


out_Model_List <- list()
for (ScenariosS in 1:length(unique(output_slope$Scenario))) {
Comparison_data <- output_slope %>% group_by(Replicates,Scenario,Type_NATS) %>% 
                                    summarise(ED_Curve_accel =mean(ED_Curve_accel ),.groups = "drop" )%>%
                                    filter(Scenario==unique(output_slope$Scenario)[ScenariosS]) %>% 
                                    filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  

Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("Rand","Env","Dist"))
Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))

model <- lmerTest::lmer(ED_Curve_accel  ~Type_NATS + (1|Replicates) , data=Comparison_data)
out_Model_List[[ScenariosS]] <- summary(model)
cat("For", unique(output_slope$Scenario)[ScenariosS], 
    "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
    "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")
}

output_slope%>%
  mutate(Scenario=
           case_when(str_detect(Scenario ,"Env") ~ "Environmental-driven assembly",
                     str_detect(Scenario ,"Spa") ~ "Spatial-driven assembly",
                     str_detect(Scenario ,"Both") ~ "Both drivers assembly",
                     str_detect(Scenario ,"Null") ~ "Null assembly")) %>% 
  mutate(Scenario=factor(Scenario,levels = c("Environmental-driven assembly",
                                             "Spatial-driven assembly",
                                             "Both drivers assembly",
                                             "Null assembly"))) %>%
  group_by(Replicates,Scenario,Type_NATS) %>% 
  summarise(ED_Curve_accel =mean(ED_Curve_accel),.groups = "drop")%>%
  filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  %>%
  mutate(Type_NATS=
           case_when(str_detect(Type_NATS ,"Dist") ~ "Distance",
                     str_detect(Type_NATS ,"Env") ~ "Environment",
                     str_detect(Type_NATS ,"Rand") ~ "Random")) %>% 
  ggplot(aes(x = Type_NATS, y=ED_Curve_accel  ))+
  geom_violin(aes(fill=as.factor(Type_NATS),color=as.factor(Type_NATS)),alpha=0.1,linewidth=2)+
  geom_jitter(shape=21,aes(fill=as.factor(Type_NATS)),size=3,alpha=0.6,width = 0.1)+
  scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
  labs(fill="Additive approach",x="",y="Edge density slope")+ 
  guides(colour="none")+
  facet_wrap(.~Scenario,scale="free")+
  theme_classic()




gridExtra::grid.arrange(
  
  output_slope %>% 
    group_by(Scenario,Type_NATS) %>% 
    summarise(Std_Dev=sd(ED_Curve_accel),ED_Curve_accel=(ED_Curve_accel)) %>% 
    ggplot(aes(x = Type_NATS, y=ED_Curve_accel))+
    #geom_errorbar(aes(x = Type_NATS, y=ED_Curve_accel,ymin=ED_Curve_accel-Std_Dev ,ymax=ED_Curve_accel+Std_Dev,
    #                  colour=as.factor(Type_NATS)),alpha=0.5)+
    geom_point(shape=21,aes(fill=as.factor(Type_NATS)))+
    #geom_line(aes(colour=Type_NATS))+
    scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
    labs(fill="Type NATS")+ guides(colour="none")+
    facet_wrap(.~Scenario)+
    theme_classic(),
  
  output_slope%>% 
    group_by(Replicates,Scenario,Type_NATS) %>% 
    summarise(Std_Dev=sd(ED_Curve_accel),ED_Plateau=mean(ED_Plateau)) %>% 
    ggplot(aes(x = as.factor(Replicates), y=ED_Plateau))+
    geom_point(shape=21,aes(fill=as.factor(Type_NATS)))+
    geom_errorbar(aes(x = as.factor(Replicates), y=ED_Plateau,ymin=ED_Plateau-Std_Dev ,ymax=ED_Plateau+Std_Dev,
                      colour=as.factor(Type_NATS)),alpha=0.5)+
    #geom_line(aes(colour=Type_NATS))+
    scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
    labs(fill="Type NATS")+ guides(colour="none")+
    facet_wrap(.~Scenario)+
    theme_classic(),
  
  output_slope%>% 
    group_by(Scenario,Type_NATS) %>% 
    summarise(Std_Dev=sd(GrStr_Obs_vs_Rand),GrStr_Obs_vs_Rand=mean(GrStr_Obs_vs_Rand)) %>% 
    ggplot(aes(x = Type_NATS, y=GrStr_Obs_vs_Rand ))+
    geom_point(shape=21,aes(fill=as.factor(Type_NATS)))+
    geom_errorbar(aes(x = Type_NATS, y=GrStr_Obs_vs_Rand ,ymin=GrStr_Obs_vs_Rand -Std_Dev ,ymax=GrStr_Obs_vs_Rand +Std_Dev,
                      colour=as.factor(Type_NATS)),alpha=0.5)+
    geom_line(aes(colour=Type_NATS))+
    scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
    labs(fill="Type NATS")+ guides(colour="none")+
    facet_wrap(.~Scenario, scales="free")+
    theme_classic(),
  
  nrow=3)



