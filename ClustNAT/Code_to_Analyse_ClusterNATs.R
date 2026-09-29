
library(tidyverse);library(viridis)
library(drc);library(nlme);library(statforbiology)

load("ClustNAT/ClusterNATs.RData")

Merged_Final_Outputs <- bind_rows(
  Final_Output$lake$Out_NAT %>% mutate(Ecosy="Lake"),
  Final_Output$river$Out_NAT %>% mutate(Ecosy="Stream")
  )

Merged_Final_Outputs <- Merged_Final_Outputs %>% mutate(Type_NATS=ifelse(Type_NATS=="Random","RO",
                                          ifelse(Type_NATS=="Distance","SD",
                                          ifelse(Type_NATS=="Environment","ES","Error")))) %>% 
                         mutate(Type_NATS=factor(Type_NATS,levels=c("SD","ES","RO")))

## NATs dataset ####

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
      random <- Sceni_full_ouput_temp %>% filter(Type_NATS=="RO") %>% # ,iter==iter[itera]
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

Letters_Plot <- list(paste(LETTERS[1:2],")",sep=""),paste(LETTERS[3:4],")",sep=""))
  
# Anàlisi significació global
out_Model_List <- list()
out_Plot_List <- list()
size_effects <- data.frame()
#std.error <- function(x) sd(x)/sqrt(length(x))  
for (Ecosys in 1:2) {
Comparison_data <- output_slope %>%
  filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
  group_by(Year,Type_NATS) %>% 
  summarise(ED_Curve_accel =mean(ED_Curve_accel ),.groups = "drop" )%>%
  filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  

year <- unique(Comparison_data$Year)

Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("RO","ES","SD"))
Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))

model <- lmerTest::lmer(ED_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
WholeModel_out <- summary(model)
cat("For the", unique(Merged_Final_Outputs$Ecosy)[Ecosys], 
    "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
    "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")

out_size_effects <- data.frame(
  "Ecosy"=unique(Merged_Final_Outputs$Ecosy)[Ecosys],
  "RO"=summary(model)$coefficients[1,1:2],
  "ES"=summary(model)$coefficients[2,1:2],
  "SD"=summary(model)$coefficients[3,1:2]) %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(3:5)

size_effects <- size_effects %>% bind_rows(out_size_effects)

# Anàlisi corrent cada 5 anys (comparem com son les diferències en packs de 5 anys)
output_slope_Running_Mean <- data.frame()
Model_YearPack <- list()
for (Year_beg in 1:(length(year)-4)) {
Posit_beg <- Year_beg
Posit_end <- Posit_beg+5

Comparison_data <- output_slope %>% 
  filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
  group_by(Year,Type_NATS) %>% 
  filter(Year%in%year[Posit_beg:Posit_end]) %>% 
  summarise(ED_Curve_accel =mean(ED_Curve_accel ),.groups = "drop" )%>%
  filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  



Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("RO","ES","SD"))
Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))

model <- lmerTest::lmer(ED_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
Model_YearPack[[Year_beg]] <- summary(model)
cat("For the", unique(Merged_Final_Outputs$Ecosy)[Ecosys], 
    "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
    "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")

Sign_df <- data.frame(Type_NATS=c("RO","ES","SD"), 
                      Sign=c("X.NoSign",
                             ifelse(summary(model)$coefficients[2,5]<0.05,"Sign","X.NoSign"),
                             ifelse(summary(model)$coefficients[3,5]<0.05,"Sign","X.NoSign")))
output_slope_Running_Mean <- bind_rows(output_slope_Running_Mean,
                             Comparison_data%>% group_by(Type_NATS) %>%
                                                summarise(Sd_ED_Curve_accel =sd(ED_Curve_accel),
                                                          ED_Curve_accel =mean(ED_Curve_accel),
                                                          .groups = "drop") %>%
                                                mutate(Year=paste(year[Posit_beg:Posit_end][1],year[Posit_beg:Posit_end][5],sep="-")) %>%
                                                left_join(Sign_df,by="Type_NATS"))

}
  
out_Model_List[[Ecosys]] <-list("Full_Model"=WholeModel_out,
                                "YearPakc_Model"=Model_YearPack) 

Test_eff_size <- out_size_effects%>%mutate(name=factor(name,levels=c("SD","ES","RO"))) %>% 
                          pivot_wider(names_from = rowname, values_from  = value) 
  
Ref_Eff_Size <- Test_eff_size %>% filter(name=="RO")

out_Plot_List[[Ecosys]] <- gridExtra::arrangeGrob(
    # output_slope%>%
    #   filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
    #   group_by(Year,Type_NATS) %>% 
    #   summarise(ED_Curve_accel =mean(ED_Curve_accel),.groups = "drop")%>%
    #   filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  %>% 
    #   #group_by(Scenario,Type_NATS) %>% 
    #   #summarise(Std_Dev=sd(ED_Curve_accel),ED_Curve_accel=ED_Curve_accel) %>% 
    #   ggplot(aes(x = Type_NATS, y=ED_Curve_accel  ))+
    #   #geom_errorbar(aes(x = Type_NATS, y=ED_Curve_accel,ymin=ED_Curve_accel-Std_Dev ,ymax=ED_Curve_accel+Std_Dev,
    #   #                  colour=as.factor(Type_NATS)),alpha=0.5)+
    #   geom_violin(aes(fill=as.factor(Type_NATS),color=as.factor(Type_NATS)),alpha=0.2,linewidth=2)+
    #   geom_jitter(shape=21,aes(fill=as.factor(Type_NATS)),size=2,alpha=0.6,width = 0.1)+
    #   scale_y_continuous(limits=c(0.08,0.25))+
    #   scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
    #   labs(title=unique(Merged_Final_Outputs$Ecosy)[Ecosys],
    #        subtitle = "A)",
    #        x="",y="Edge density slope")+ 
    #   guides(colour="none",fill="none")+
    #   theme_classic(),

  Test_eff_size %>% #filter(name!="RB") %>% 
    left_join(Ref_Eff_Size, by=c("Ecosy")) %>% 
    mutate(Size_Effect=Estimate.y+Estimate.x) %>% 
    mutate(Size_Effect=ifelse(name.x=="RO",Estimate.y,Size_Effect)) %>% 
    
    ggplot()+
    geom_rect(data = Ref_Eff_Size, aes(xmax =Estimate+`Std. Error`,xmin = Estimate-`Std. Error`,ymin=0,ymax=Inf),
              alpha=0.1,fill="grey30")+#viridis(n = 3,option = "D",direction = 1)[3])+
    geom_linerange(aes(x=Size_Effect,y=name.x, xmin=Size_Effect-`Std. Error.x`,xmax=Size_Effect+`Std. Error.x`),
                   linewidth=1.2)+
    geom_point(aes(y=name.x,x=Size_Effect,fill=name.x,colour=name.x),size=4,shape=21)+
    scale_fill_viridis(option = "D",direction = 1,discrete = T)+
    scale_colour_viridis(option = "D",direction = 1,discrete = T)+
    scale_x_continuous(limits=c(0.13,0.20))+
    labs(subtitle= paste(Letters_Plot[[Ecosys]][1],unique(Merged_Final_Outputs$Ecosy)[Ecosys],"metacommunity"),y="",
         colour="Additive approach",fill="Additive approach",x="Effect size")+
    theme_classic()+
    theme(legend.position = "none",
        strip.text.y = element_text(angle = 0), 
          panel.background = element_rect(colour="black")),
    
    output_slope_Running_Mean%>%
      mutate(Type_NATS=factor(Type_NATS,levels=c("SD","ES","RO"))) %>% 
      ggplot(aes(x = Year, y=ED_Curve_accel))+
      geom_errorbar(aes(x = Year, y=ED_Curve_accel,ymin=ED_Curve_accel-Sd_ED_Curve_accel,ymax=ED_Curve_accel+Sd_ED_Curve_accel ,
                        colour=as.factor(Type_NATS)),alpha=0.5)+
      geom_line(aes(colour=Type_NATS),linewidth=1.5, alpha=0.2)+
      geom_point(shape=21,aes(fill=as.factor(Type_NATS),colour=as.factor(Type_NATS)),size=4,alpha=0.6,stroke=1)+
      geom_jitter(inherit.aes = F,shape=8,size=2,height = -0.05,width = 0.25,stroke=1.5,
                 aes(y=rep(0.21,(length(year)-4)*3),x=Year,
                     colour=as.factor(Type_NATS),alpha=Sign))+
      scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
      scale_alpha_manual(values=c(1,0.01))+
      scale_y_continuous(limits=c(0.08,0.25))+
      labs(subtitle= paste(Letters_Plot[[Ecosys]][2],unique(Merged_Final_Outputs$Ecosy)[Ecosys],"metacommunity"),
          fill="Additive approach",y="Edge density accumulation rate",x="Period")+ 
      guides(colour="none",alpha="none")+
      #facet_wrap(.~Type_NATS)+labs(fill="Year")+
      theme_classic()+
      theme(axis.text.x = element_text(angle=45,hjust = 1)),
    
ncol=2 ,widths=c(1,2))

} # Ecosy

png(filename = "ClustNAT/NATs_results.png", width = 1000*3.5,height = 700*3.5,units = "px",res = 300)
gridExtra::grid.arrange(out_Plot_List[[1]],out_Plot_List[[2]],nrow=2)
dev.off()

## FUNCTIONAL dataset ####
# Treatment and analysis of the functional traits
DF_FunIndices <- read.csv2(file="ClustNAT/FunctionalIndices_Tot.csv") %>% 
                  mutate(Type_NATs=ifelse(Type_NATs=="Random","RO",
                                   ifelse(Type_NATs=="Distance","SD",
                                   ifelse(Type_NATs=="Environment","ES","Error")))) %>% 
                  mutate(Type_NATs=factor(Type_NATs,levels=c("SD","ES","RO"))) %>% 
                  rename("Type_NATS"="Type_NATs") %>%
                  separate(Combi_ID, c("Year","n_sites","iter")) %>% 
                  mutate(System=ifelse(System=="lake","Lake",
                                ifelse(System=="river","Stream","ERROR"))) %>% 
                  rename("Ecosy"="System") %>% 
                  mutate(Year=as.numeric(Year),n_sites=as.numeric(n_sites))

output_slope <- data.frame()
for (Ecosys in 1:2) {
  Final_O_put <- DF_FunIndices %>% filter(Ecosy==unique(DF_FunIndices$Ecosy)[Ecosys]) %>% 
                                   dplyr::select(Ecosy,Type_NATS,Year,n_sites,iter,FRic)
  
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
      random <- Sceni_full_ouput_temp %>% filter(Type_NATS=="RO") %>% # ,iter==iter[itera]
                group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      
      # model_Random <- lm(mean_grStre~(as.numeric(n_sites)), data=random)
      # Model_nats <- predict(model_Random)
      
      #lm <- summary(lm(mean_grStre~(as.numeric(n_sites)), data=full_output_temp_temp))
      #slope <- lm$coefficients[2,1]
      # full_output_temp_temp <- full_ouput_temp %>%group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      # slope <- abs(Model_nats-(full_output_temp_temp$mean_grStre))
      
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
      Y <-  full_output_temp_temp$FRic
      X <- full_output_temp_temp$n_sites
      control1 <- nls.control(maxiter= 1000,tol=1e-02, warnOnly=TRUE)
      model_S <- nls(Y~NLS.asymReg(X, init, m, plateau),control=control1)
      output_slope2 <- data.frame("Ecosy"=unique(Merged_Final_Outputs$Ecosy)[Ecosys],
                                  "Year"=year[sceni],
                                  "Type_NATS"=type[ind_type],
                                  #"iter"=iter[itera],
                                  #"site"=full_output_temp_temp$n_sites,
                                  #"GrStr_Obs_vs_Rand"=(slope),
                                  "FRic_Curve_accel"=coefficients(model_S)[2],
                                  "FRic_Plateau"=coefficients(model_S)[3])
      
      output_slope <- rbind(output_slope,output_slope2)
      #}#itera
    }#Type Nats
  }# sceni
  
} # Ecosy

Letters_Plot <- list(paste(LETTERS[1:2],")",sep=""),paste(LETTERS[3:4],")",sep=""))

# Anàlisi significació global
out_Model_List <- list()
out_Plot_List <- list()
size_effects <- data.frame()
#std.error <- function(x) sd(x)/sqrt(length(x))  
for (Ecosys in 1:2) {
  Comparison_data <- output_slope %>%
    filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
    group_by(Year,Type_NATS) %>% 
    summarise(FRic_Curve_accel =mean(FRic_Curve_accel ),.groups = "drop" )%>%
    filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  
  
  year <- unique(Comparison_data$Year)
  
  Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("RO","ES","SD"))
  Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))
  
  model <- lmerTest::lmer(FRic_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
  WholeModel_out <- summary(model)
  cat("For the", unique(Merged_Final_Outputs$Ecosy)[Ecosys], 
      "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
      "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")
  
  out_size_effects <- data.frame(
    "Ecosy"=unique(Merged_Final_Outputs$Ecosy)[Ecosys],
    "RO"=summary(model)$coefficients[1,1:2],
    "ES"=summary(model)$coefficients[2,1:2],
    "SD"=summary(model)$coefficients[3,1:2]) %>% 
    tibble::rownames_to_column() %>% 
    pivot_longer(3:5)
  
  size_effects <- size_effects %>% bind_rows(out_size_effects)
  
  # Anàlisi corrent cada 5 anys (comparem com son les diferències en packs de 5 anys)
  output_slope_Running_Mean <- data.frame()
  Model_YearPack <- list()
  for (Year_beg in 1:(length(year)-4)) {
    Posit_beg <- Year_beg
    Posit_end <- Posit_beg+5
    
    Comparison_data <- output_slope %>% 
      filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
      group_by(Year,Type_NATS) %>% 
      filter(Year%in%year[Posit_beg:Posit_end]) %>% 
      summarise(FRic_Curve_accel =mean(FRic_Curve_accel ),.groups = "drop" )%>%
      filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  
    
    
    
    Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("RO","ES","SD"))
    Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))
    
    model <- lmerTest::lmer(FRic_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
    Model_YearPack[[Year_beg]] <- summary(model)
    cat("For the", unique(Merged_Final_Outputs$Ecosy)[Ecosys], 
        "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
        "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")
    
    Sign_df <- data.frame(Type_NATS=c("RO","ES","SD"), 
                          Sign=c("X.NoSign",
                                 ifelse(summary(model)$coefficients[2,5]<0.05,"Sign","X.NoSign"),
                                 ifelse(summary(model)$coefficients[3,5]<0.05,"Sign","X.NoSign")))
    output_slope_Running_Mean <- bind_rows(output_slope_Running_Mean,
                                           Comparison_data%>% group_by(Type_NATS) %>%
                                             summarise(Sd_FRic_Curve_accel =sd(FRic_Curve_accel),
                                                       FRic_Curve_accel =mean(FRic_Curve_accel),
                                                       .groups = "drop") %>%
                                             mutate(Year=paste(year[Posit_beg:Posit_end][1],year[Posit_beg:Posit_end][2],sep="-")) %>%
                                             left_join(Sign_df,by="Type_NATS"))
  }
  
  out_Model_List[[Ecosys]] <-list("Full_Model"=WholeModel_out,
                                  "YearPakc_Model"=Model_YearPack) 
  
  Test_eff_size <- out_size_effects%>%mutate(name=factor(name,levels=c("SD","ES","RO"))) %>% 
    pivot_wider(names_from = rowname, values_from  = value) 
  
  Ref_Eff_Size <- Test_eff_size %>% filter(name=="RO")
  
  out_Plot_List[[Ecosys]] <- gridExtra::arrangeGrob(
    # output_slope%>%
    #   filter(Ecosy==unique(Merged_Final_Outputs$Ecosy)[Ecosys]) %>% 
    #   group_by(Year,Type_NATS) %>% 
    #   summarise(FRic_Curve_accel =mean(FRic_Curve_accel),.groups = "drop")%>%
    #   filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  %>% 
    #   #group_by(Scenario,Type_NATS) %>% 
    #   #summarise(Std_Dev=sd(FRic_Curve_accel),FRic_Curve_accel=FRic_Curve_accel) %>% 
    #   ggplot(aes(x = Type_NATS, y=FRic_Curve_accel  ))+
    #   #geom_errorbar(aes(x = Type_NATS, y=FRic_Curve_accel,ymin=FRic_Curve_accel-Std_Dev ,ymax=FRic_Curve_accel+Std_Dev,
    #   #                  colour=as.factor(Type_NATS)),alpha=0.5)+
    #   geom_violin(aes(fill=as.factor(Type_NATS),color=as.factor(Type_NATS)),alpha=0.2,linewidth=2)+
    #   geom_jitter(shape=21,aes(fill=as.factor(Type_NATS)),size=2,alpha=0.6,width = 0.1)+
    #   #scale_y_continuous(limits=c(0.08,0.25))+
    #   scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
    #   labs(title=unique(Merged_Final_Outputs$Ecosy)[Ecosys],
    #        subtitle = "A)",
    #        x="",y="Edge density slope")+ 
    #   guides(colour="none",fill="none")+
    #   theme_classic(),
    Test_eff_size %>% #filter(name!="RB") %>% 
      left_join(Ref_Eff_Size, by=c("Ecosy")) %>% 
      mutate(Size_Effect=Estimate.y+Estimate.x) %>% 
      mutate(Size_Effect=ifelse(name.x=="RO",Estimate.y,Size_Effect)) %>% 
      
      ggplot()+
      geom_rect(data = Ref_Eff_Size, aes(xmax =Estimate+`Std. Error`,xmin = Estimate-`Std. Error`,ymin=0,ymax=Inf),
                alpha=0.1,fill="grey30")+#viridis(n = 3,option = "D",direction = 1)[3])+
      geom_linerange(aes(x=Size_Effect,y=name.x, xmin=Size_Effect-`Std. Error.x`,xmax=Size_Effect+`Std. Error.x`),
                     linewidth=1.2)+
      geom_point(aes(y=name.x,x=Size_Effect,fill=name.x,colour=name.x),size=4,shape=21)+
      scale_fill_viridis(option = "D",direction = 1,discrete = T)+
      scale_colour_viridis(option = "D",direction = 1,discrete = T)+
      scale_x_continuous(limits=c(0.13,0.20))+
      labs(subtitle= paste(Letters_Plot[[Ecosys]][1],unique(Merged_Final_Outputs$Ecosy)[Ecosys],"metacommunity"),
           y="",x="Effect size",
           colour="Additive approach",fill="Additive approach")+
      theme_classic()+
      theme(legend.position = "none",
            strip.text.y = element_text(angle = 0), 
            panel.background = element_rect(colour="black")),
    
    output_slope_Running_Mean%>%
      mutate(Type_NATS=factor(Type_NATS,levels=c("SD","ES","RO"))) %>% 
      ggplot(aes(x = Year, y=FRic_Curve_accel))+
      geom_errorbar(aes(x = Year, y=FRic_Curve_accel,ymin=FRic_Curve_accel-Sd_FRic_Curve_accel,ymax=FRic_Curve_accel+Sd_FRic_Curve_accel ,
                        colour=as.factor(Type_NATS)),alpha=0.5)+
      geom_line(aes(colour=Type_NATS),linewidth=1.5, alpha=0.2)+
      geom_point(shape=21,aes(fill=as.factor(Type_NATS),colour=as.factor(Type_NATS)),size=4,alpha=0.6,stroke=1)+
      geom_jitter(inherit.aes = F,shape=8,size=2,height = -0.05,width = 0.25,stroke=1.5,
                  aes(y=rep(0.21,(length(year)-4)*3),x=Year,
                      colour=as.factor(Type_NATS),alpha=Sign))+
      scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
      scale_alpha_manual(values=c(1,0.01))+
      #scale_y_continuous(limits=c(0.08,0.25))+
      labs(subtitle= paste(Letters_Plot[[Ecosys]][2],unique(Merged_Final_Outputs$Ecosy)[Ecosys],"metacommunity"),
           fill="Additive approach",y="Edge density accumulation rate",x="Period")+ 
      guides(colour="none",alpha="none")+
      #facet_wrap(.~Type_NATS)+labs(fill="Year")+
      theme_classic()+
      theme(axis.text.x = element_text(angle=45,hjust = 1)),
    
    ncol=2 ,widths=c(1,2))
  
} # Ecosy

png(filename = "ClustNAT/FRic_results.png", width = 1000*3.5,height = 700*3.5,units = "px",res = 300)
gridExtra::grid.arrange(out_Plot_List[[1]],out_Plot_List[[2]],nrow=2,top="FRic values for the three additive approaches")
dev.off()


## Accumulation curves ####
png(filename = "ClustNAT/All_Curves_NATs_Lakes.png",width =2000*3,height = 2000*3,res = 300) 
Merged_Final_Outputs %>%
  filter(Ecosy=="Lake") %>% 
  #filter(Year==2000) %>% 
  group_by(Year,Type_NATS,n_sites)  %>% 
  mutate(mean_edge_dens=mean(edge_dens),mean_mean_grStre=mean(mean_grStre)) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=edge_dens,colour=Type_NATS),alpha=0.1)+
  geom_line(aes(x=n_sites,y=edge_dens,colour=Type_NATS,group=iter),alpha=0.2)+
  geom_line(aes(x=n_sites,y=mean_edge_dens,group=iter),colour="black", linewidth=1)+
  scale_color_viridis(discrete=T)+
  labs(title="Lakes edge density cummulative curves",
       colour="Additive approach")+
  theme_classic()+
  facet_wrap(Type_NATS~Year,ncol=7,strip.position = 'top')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0)
  )
dev.off()

png(filename = "ClustNAT/All_Curves_NATs_Rivers.png",width =2000*3,height = 2000*3,res = 300) 
Merged_Final_Outputs %>%
  filter(Ecosy=="Stream") %>% 
  #filter(Year==2000) %>% 
  group_by(Year,Type_NATS,n_sites)  %>% 
  mutate(mean_edge_dens=mean(edge_dens),mean_mean_grStre=mean(mean_grStre)) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=edge_dens,colour=Type_NATS),alpha=0.1)+
  geom_line(aes(x=n_sites,y=edge_dens,colour=Type_NATS,group=iter),alpha=0.2)+
  geom_line(aes(x=n_sites,y=mean_edge_dens,group=iter),colour="black", linewidth=1)+
  scale_color_viridis(discrete=T)+
  labs(title="Stream edge density cummulative curves",
       colour="Additive approach")+
  theme_classic()+
  facet_wrap(Type_NATS~Year,ncol=6,strip.position = 'top')+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(angle = 0)
  )
dev.off()
