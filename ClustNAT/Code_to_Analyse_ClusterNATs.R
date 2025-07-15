
library(tidyverse);library(viridis)
library(drc);library(nlme);library(statforbiology)

load("ClustNAT/ClusterNATs.RData")

Merged_Final_Outputs <- bind_rows(
  Final_Output$lake$Out_NAT %>% mutate(Ecosy="Lake"),
  Final_Output$river$Out_NAT %>% mutate(Ecosy="River")
  )

output_slope <- data.frame()
for (Ecosys in 1:2) {
Final_O_put <- Merged_Final_Outputs %>% filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys])
  
year <- unique(Final_O_put$Year)
type <- unique(Final_O_put$Type_NATS)
iter <- unique(Final_O_put$iter)
for (sceni in 1:length(year)) {
    Sceni_full_ouput_temp <-  Final_O_put %>% filter(Year==year[sceni])
    for (ind_type in 1:length(type)) {
      #for (itera in 1:length(iter)) {
      # We create a "temporary" file filtered according to the TypeNat selected. 
      full_ouput_temp <- Sceni_full_ouput_temp %>% 
                          filter(Type_NATS==type[ind_type])#,iter==iter[itera]) # En cas de voler filtrar per iteració
      random <- Sceni_full_ouput_temp %>% filter(Type_NATS=="Random") %>% # ,iter==iter[itera]
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
      output_slope2 <- data.frame("Ecosy"=unique(Merged_Final_Outputs$Ecosy)[Ecosys],
                                  "Year"=year[sceni],
                                  "Type_NATS"=type[ind_type],
                                  #"iter"=iter[itera],
                                  #"site"=full_output_temp_temp$n_sites,
                                  "GrStr_Obs_vs_Rand"=(slope),
                                  "ED_Curve_accel"=coefficients(model_S)[2],
                                  "ED_Plateau"=coefficients(model_S)[3])
      
      output_slope <- rbind(output_slope,output_slope2)
    #}#itera
  }#Type Nats
}# sceni

} # Ecosy
  
# Anàlisi significació global
out_Model_List <- list()
out_Plot_List <- list()
#std.error <- function(x) sd(x)/sqrt(length(x))  
for (Ecosys in 1:2) {
Comparison_data <- output_slope %>%
  filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
  group_by(Year,Type_NATS) %>% 
  summarise(ED_Curve_accel =mean(ED_Curve_accel ),.groups = "drop" )%>%
  filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  

Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("Random","Environment","Distance"))
Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))

model <- lmerTest::lmer(ED_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
WholeModel_out <- summary(model)
cat("For the", unique(Merged_Final_Outputs$Ecosy)[Ecosys], 
    "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
    "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")

# Anàlisi corrent cada 5 anys (comparem com son les diferències en packs de 5 anys)
output_slope_Running_Mean <- data.frame()
Model_YearPack <- list()
for (Year_beg in 1:(length(year)-5)) {
Posit_beg <- Year_beg
Posit_end <- Posit_beg+5

Comparison_data <- output_slope %>% 
  filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
  group_by(Year,Type_NATS) %>% 
  filter(Year%in%year[Posit_beg:Posit_end]) %>% 
  summarise(ED_Curve_accel =mean(ED_Curve_accel ),.groups = "drop" )%>%
  filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  



Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("Random","Environment","Distance"))
Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))

model <- lmerTest::lmer(ED_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
Model_YearPack[[Year_beg]] <- summary(model)
cat("For the", unique(Merged_Final_Outputs$Ecosy)[Ecosys], 
    "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
    "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")

Sign_df <- data.frame(Type_NATS=c("Random","Environment","Distance"), 
                      Sign=c("X.NoSign",
                             ifelse(summary(model)$coefficients[2,5]<0.05,"Sign","X.NoSign"),
                             ifelse(summary(model)$coefficients[3,5]<0.05,"Sign","X.NoSign")))
output_slope_Running_Mean <- bind_rows(output_slope_Running_Mean,
                             Comparison_data%>% group_by(Type_NATS) %>%
                                                summarise(Sd_ED_Curve_accel =sd(ED_Curve_accel),
                                                          ED_Curve_accel =mean(ED_Curve_accel),
                                                          .groups = "drop") %>%
                                                mutate(Year=Posit_beg) %>%
                                                left_join(Sign_df,by="Type_NATS"))
}
  
out_Model_List[[Ecosys]] <-list("Full_Model"=WholeModel_out,
                                "YearPakc_Model"=Model_YearPack) 

out_Plot_List[[Ecosys]] <- gridExtra::arrangeGrob(
    output_slope%>%
      filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
      group_by(Year,Type_NATS) %>% 
      summarise(ED_Curve_accel =mean(ED_Curve_accel),.groups = "drop")%>%
      filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  %>% 
      #group_by(Scenario,Type_NATS) %>% 
      #summarise(Std_Dev=sd(ED_Curve_accel),ED_Curve_accel=ED_Curve_accel) %>% 
      ggplot(aes(x = Type_NATS, y=ED_Curve_accel  ))+
      #geom_errorbar(aes(x = Type_NATS, y=ED_Curve_accel,ymin=ED_Curve_accel-Std_Dev ,ymax=ED_Curve_accel+Std_Dev,
      #                  colour=as.factor(Type_NATS)),alpha=0.5)+
      geom_violin(aes(fill=as.factor(Type_NATS),color=as.factor(Type_NATS)),alpha=0.2,linewidth=2)+
      geom_jitter(shape=21,aes(fill=as.factor(Type_NATS)),size=2,alpha=0.6,width = 0.1)+
      scale_y_continuous(limits=c(0.08,0.25))+
      scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
      labs(title=unique(Merged_Final_Outputs$Ecosy)[Ecosys],
           subtitle = "A)",
           x="",y="Edge density slope")+ 
      guides(colour="none",fill="none")+
      theme_classic(),
    
    output_slope_Running_Mean%>%
      ggplot(aes(x = Year, y=ED_Curve_accel))+
      geom_errorbar(aes(x = Year, y=ED_Curve_accel,ymin=ED_Curve_accel-Sd_ED_Curve_accel,ymax=ED_Curve_accel+Sd_ED_Curve_accel ,
                        colour=as.factor(Type_NATS)),alpha=0.5)+
      geom_line(aes(colour=Type_NATS),linewidth=1.5, alpha=0.2)+
      geom_point(shape=21,aes(fill=as.factor(Type_NATS),colour=as.factor(Type_NATS)),size=4,alpha=0.6,stroke=1)+
      geom_jitter(inherit.aes = F,shape=8,size=4,height = -0.05,width = 0.25,stroke=1.5,
                 aes(y=rep(0.21,(length(year)-5)*3),x=Year,
                     colour=as.factor(Type_NATS),alpha=Sign))+
      scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
      scale_alpha_manual(values=c(1,0.01))+
      scale_y_continuous(limits=c(0.08,0.25))+
      labs(title="",subtitle ="B)",
          fill="Additive approach",y="Edge density slope",x="Time")+ 
      guides(colour="none",alpha="none")+
      #facet_wrap(.~Type_NATS)+labs(fill="Year")+
      theme_classic()+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()),
    
ncol=2 ,widths=c(1,2))

} # Ecosy

gridExtra::grid.arrange(out_Plot_List[[1]],out_Plot_List[[2]],nrow=2)


















gridExtra::grid.arrange(

output_slope%>% 
  filter(Type_NATS%in%c("Random","Distance","Environment")) %>% 
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


