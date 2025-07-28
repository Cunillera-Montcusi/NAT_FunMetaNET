
#install.packages("devtools")
#devtools::install_github("onofriAndreaPG/aomisc")
library(tidyverse);library(viridis)
library(drc);library(nlme);library(statforbiology)

#calcul de pendents LLACS
load("NATs_lakes_rand_def.RData")
rand_ecosyst <- output
load("NATs_lakes_env_def.RData")
env_ecosyst <- output
load("NATs_lakes_Dist_def.RData")
dist_ecosyst <- output
load("NATs_lakes_DisEnv_def.RData")
DisEnv_ecosyst <- output
load("NATs_lakes_DisDist_def.RData")
DisDist_ecosyst <- output

#calcul de pendents RIUS
load("NATs_riv_rand_def.RData")
rand_ecosyst <- output
load("NATs_riv_env_def.RData")
env_ecosyst <- output
load("NATs_riv_Dist_def.RData")
dist_ecosyst <- output
load("NATs_riv_disEnv_def.RData")
DisEnv_ecosyst <- output
load("NATs_riv_DisDits_def.RData")
DisDist_ecosyst <- output


full_output <- bind_rows(
  rand_ecosyst%>% mutate(Type_NATS="Random"),
  env_ecosyst%>% mutate(Type_NATS="Environment"),
  dist_ecosyst%>% mutate(Type_NATS="Distance"),
  DisEnv_ecosyst%>% mutate(Type_NATS="Dis_Environment"),
  DisDist_ecosyst%>% mutate(Type_NATS="Dis_Distance"))

full_output <- full_output %>% filter(Type_NATS%in%c("Random","Environment","Distance"))

year <- unique(full_output$Year)
type <- unique(full_output$Type_NATS)
iter <- unique(full_output$iter)
output_slope <- data.frame()

for (sceni in 1:length(year)) {
  Sceni_full_ouput_temp <-  full_output %>% filter(Year==year[sceni])
  for (ind_type in 1:length(type)) {
  #for (itera in 1:length(iter)) {
  # We create a "temporary" file filtered according to the TypeNat selected. 
  full_ouput_temp <- Sceni_full_ouput_temp %>% filter(Type_NATS==type[ind_type]) # En cas de voler filtrar per iteració #,iter==iter[itera]
  random <- Sceni_full_ouput_temp %>% filter(Type_NATS=="Random") %>% # ,iter==iter[itera]
                                      group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
  
  model_Random <- lm(mean_grStre~(as.numeric(n_sites)), data=random)
  Model_nats <- predict(model_Random)  
  
  
  full_output_temp_temp <- full_ouput_temp %>%group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
  slope <- abs(Model_nats-(full_output_temp_temp$mean_grStre))
  
  full_output_temp_temp <- full_ouput_temp
  Y <-  full_output_temp_temp$edge_dens
  X <- full_output_temp_temp$n_sites
  control1 <- nls.control(maxiter= 1000,tol=1e-02, warnOnly=TRUE)
  model_S <- nls(Y~NLS.asymReg(X, init, m, plateau),control=control1)
  output_slope2 <- data.frame("Year"=year[sceni],
                              "Type_NATS"=type[ind_type],
                              #"iter"=iter[itera],
                              #"site"=full_output_temp_temp$n_sites,
                              "GrStr_Obs_vs_Rand"=(slope),
                              "ED_Curve_accel"=coefficients(model_S)[2],
                              "ED_Plateau"=coefficients(model_S)[3])
    
    output_slope <- rbind(output_slope,output_slope2)
  }
}

# Anàlisi significació global
out_Model_List <- list()
Comparison_data <- output_slope %>% group_by(Year,Type_NATS) %>% 
    summarise(ED_Curve_accel =mean(ED_Curve_accel ),.groups = "drop" )%>%
    filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  
  
Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("Random","Environment","Distance"))
Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))
  
model <- lmerTest::lmer(ED_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
summary(model)
cat("For the Lakes", 
    "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
    "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")

# Anàlisi corrent cada 5 anys (comparem com son les diferències en packs de 5 anys)
std.error <- function(x) sd(x)/sqrt(length(x))
output_slope_Running_Mean <- data.frame()
for (Year_beg in 1:(length(year)-5)) {
  Posit_beg <- Year_beg
  Posit_end <- Posit_beg+5
  
  Comparison_data <- output_slope %>% group_by(Year,Type_NATS) %>% 
    filter(Year%in%year[Posit_beg:Posit_end]) %>% 
    summarise(ED_Curve_accel =mean(ED_Curve_accel ),.groups = "drop" )%>%
    filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  
  
  output_slope_Running_Mean <- bind_rows(output_slope_Running_Mean,
                                         Comparison_data%>% group_by(Type_NATS) %>% 
                                           summarise(Sd_ED_Curve_accel =std.error(ED_Curve_accel),
                                                     ED_Curve_accel =mean(ED_Curve_accel),
                                                     .groups = "drop" ) %>% 
                                           mutate(Year=Posit_beg))
  
  Comparison_data$Type_NATS <- factor(Comparison_data$Type_NATS,levels=c("Random","Environment","Distance"))
  Comparison_data <- within(Comparison_data, Type_NATS <- relevel(Type_NATS, ref = 1))
  
  model <- lmerTest::lmer(ED_Curve_accel  ~Type_NATS + (1|Year) , data=Comparison_data)
  summary(model)
  cat("For the Rivers", 
      "Env is", ifelse(summary(model)$coefficients[2,5]<0.05,"Different","NOdifferent"),  
      "and Disp is",ifelse(summary(model)$coefficients[3,5]<0.05,"Different","NOdifferent"), "\n")
}


gridExtra::grid.arrange(
  gridExtra::arrangeGrob(

output_slope%>% group_by(Year,Type_NATS) %>% 
  summarise(ED_Curve_accel =mean(ED_Curve_accel ))%>%
  filter(!Type_NATS%in%c("Dis_Env","Dis_Dist"))  %>% 
  #group_by(Scenario,Type_NATS) %>% 
  #summarise(Std_Dev=sd(ED_Curve_accel),ED_Curve_accel=ED_Curve_accel) %>% 
  ggplot(aes(x = Type_NATS, y=ED_Curve_accel  ))+
  #geom_errorbar(aes(x = Type_NATS, y=ED_Curve_accel,ymin=ED_Curve_accel-Std_Dev ,ymax=ED_Curve_accel+Std_Dev,
  #                  colour=as.factor(Type_NATS)),alpha=0.5)+
  geom_violin(aes(fill=as.factor(Type_NATS),color=as.factor(Type_NATS)),alpha=0.2,linewidth=2)+
  geom_jitter(shape=21,aes(fill=as.factor(Type_NATS)),size=2,alpha=0.6,width = 0.1)+
  scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
  labs(title="Rivers", fill="Type NATS",x="")+ guides(colour="none",fill="none")+
  theme_classic(),

output_slope_Running_Mean%>%
  ggplot(aes(x = Year, y=ED_Curve_accel))+
  geom_errorbar(aes(x = Year, y=ED_Curve_accel,ymin=ED_Curve_accel-Sd_ED_Curve_accel,ymax=ED_Curve_accel+Sd_ED_Curve_accel ,
                    colour=as.factor(Type_NATS)),alpha=0.5)+
  geom_line(aes(colour=Type_NATS),linewidth=1.5, alpha=0.2)+
  geom_point(shape=21,aes(fill=as.factor(Type_NATS),colour=as.factor(Type_NATS)),size=4,alpha=0.6,stroke=1)+
  scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
  labs(fill="Type NATS")+ guides(colour="none")+
  #facet_wrap(.~Type_NATS)+labs(fill="Year")+
  theme_classic(),

ncol=2 ,widths=c(1,2))
)






output_slope%>% group_by(Year,Type_NATS) %>% 
  summarise(ED_Curve_accel =mean(ED_Curve_accel ))%>%
  ggplot(aes(x = Year, y=ED_Curve_accel))+
  #geom_errorbar(aes(x = Year, y=ED_Curve_accel,ymin=ED_Curve_accel-Std_Dev ,ymax=ED_Curve_accel+Std_Dev,
  #                  colour=as.factor(Type_NATS)),alpha=0.5)+
  geom_line(aes(colour=Type_NATS),linewidth=1.5, alpha=0.2)+
  geom_point(shape=21,aes(fill=as.factor(Type_NATS),colour=as.factor(Type_NATS)),size=4,alpha=0.6,stroke=1)+
  scale_fill_viridis(discrete = T)+  scale_colour_viridis(discrete = T)+
  labs(fill="Type NATS")+ guides(colour="none")+
  #facet_wrap(.~Type_NATS)+labs(fill="Year")+
  theme_classic()








# Merge with richness
unique(full_output$Type_NATS)

output_slope

Match_Rich_DisEnv <- full_output_riv %>% filter(TypeNAT=="Dis_Env") %>% dplyr::select(c(Year,n_sites,iter,edge_dens)) %>% 
  mutate(rowname=paste(Year,n_sites,iter,sep = "_")) %>% 
  left_join(data.frame(important.indices) %>% tibble::rownames_to_column() %>% dplyr::select(c(rowname,NumbSpecies)),
            by="rowname") 

Y <-  Match_Rich_DisEnv$edge_dens
X <- Match_Rich_DisEnv$NumbSpecies
modelDisENv <- nls(Y~NLS.asymReg(X, init, m, plateau))
Res <- residuals(modelDisENv)

Match_Rich_DisEnv %>% mutate(Res_Rich_ED=Res) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=abs(Res_Rich_ED)))

Match_Rich_DisEnv <- Match_Rich_DisEnv%>% filter(Year==2015)
Match_Rich_Rand <- Match_Rich_Rand%>% filter(Year==2015)

#relacio n_sites amb rich i esge_dens
Match_Rich_DisEnv %>% ggplot()+
  geom_point(aes(x=n_sites,y=edge_dens))
Y <-  Match_Rich_DisEnv$edge_dens
X <- Match_Rich_DisEnv$n_sites
model1 <- nls(Y~NLS.asymReg(X, init, m, plateau))
Y <-  Match_Rich_DisEnv$NumbSpecies
X <- Match_Rich_DisEnv$n_sites
model2 <- nls(Y~NLS.asymReg(X, init, m, plateau))
predic1 <- predict(model1)
predic2 <- predict(model2)
diff_disenv <- abs(predic1-predic2)
Match_Rich_DisEnv %>% mutate(Res_diff_DisEnv=diff_disenv) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=Res_diff_DisEnv))

Match_Rich %>% ggplot()+
  geom_point(aes(x=n_sites,y=edge_dens))
Y <-  Match_Rich$edge_dens
X <- Match_Rich$n_sites
model1 <- nls(Y~NLS.asymReg(X, init, m, plateau))
Y <-  Match_Rich$NumbSpecies
X <- Match_Rich$n_sites
model2 <- nls(Y~NLS.asymReg(X, init, m, plateau))
predic1 <- predict(model1)
predic2 <- predict(model2)
diff_rand <- abs(predic1-predic2)
Match_Rich_Rand%>% mutate(Res_diff_rand=diff) %>% mutate(Res_diff_DisEnv=diff_disenv)%>%
  pivot_longer(cols=7:8) %>% group_by(Year,n_sites,name) %>%
  ggplot()+
  geom_line(aes(x=n_sites,y=value, colour=name))

#RIVERS####
load("NATs_riv_rand_def.RData")
rand_riv <- output
load("NATs_riv_env_def.RData")
env_riv <- output
load("NATs_riv_dist_def.RData")
dist_riv <- output
load("NATs_riv_disEnv_def.RData")
DisEnv_riv <- output
load("NATs_riv_DisDits_def.RData")
DisDist_riv <- output
 
full_output_riv <- bind_rows(
  rand_riv%>% mutate(TypeNAT="Rand"),
  env_riv%>% mutate(TypeNAT="Env"),
  dist_riv%>% mutate(TypeNAT="Dist"),
  DisEnv_riv%>% mutate(TypeNAT="Dis_Env"),
 DisDist_riv%>% mutate(TypeNAT="Dis_Dist"))

years <- unique(full_output_riv$Year)
type <- unique(full_output_riv$TypeNAT)
output_slope_riv <- data.frame()

### FIRST LOOP - TypeNAT
for (ind_type in 1:length(type)) {
  # We create a "temporary" file filtered according to the TypeNat selected. 
  full_ouput_temp <- full_output_riv %>% filter(TypeNAT==type[ind_type])
  # All years of that typeNAT
  year_type <- unique(full_ouput_temp$Year)
  
  ### 2ND LOOP - YEAR
  for (year in 1:length(year_type)) {
    
    random <- full_output_riv %>% filter(TypeNAT=="Rand",Year==years[year]) %>% 
                        group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
      
    model_Random <- lm(mean_grStre~(as.numeric(n_sites)), data=random)
    
    Model_nats <- predict(model_Random)
    
    #fltrem cada any
    full_output_temp_temp <- full_ouput_temp %>% filter(Year==years[year])%>%
      group_by(n_sites)%>%summarise_if(is.numeric, ~mean(.))
    #lm <- summary(lm(mean_grStre~(as.numeric(n_sites)), data=full_output_temp_temp))
    #slope <- lm$coefficients[2,1]
    
    slope <- abs(Model_nats-(full_output_temp_temp$mean_grStre))
    
    Y <-  full_output_temp_temp$edge_dens
    X <- full_output_temp_temp$n_sites
    model_S <- nls(Y~NLS.asymReg(X, init, m, plateau))
    output_slope2 <- data.frame("TypeNAT"=type[ind_type],
                                "site"=full_output_temp_temp$n_sites,
                                "Year"=year_type[year], "Slope"=(slope),
                                "m"=coefficients(model_S)[2],
                                "plat"=coefficients(model_S)[3])
  
    output_slope_riv <- rbind(output_slope_riv,output_slope2)
  }
}

#grafics####
library(viridis)
output_sviridisoutput_slope%>% group_by(TypeNAT)%>%summarise(slope=mean(m))%>%
  ggplot()+geom_point(aes(x=TypeNAT, y=slope))
output_slope %>% filter(!TypeNAT%in%"Rand") %>% 
  group_by(TypeNAT,Year)%>%summarise(LM_Residuals=mean(Slope))%>%
               ggplot()+ 
                geom_boxplot(aes(x=TypeNAT, y=LM_Residuals,fill=as.factor(TypeNAT)))+ theme_classic()+ theme(legend.position = "none")+
  labs(x="Type of addition", y="Residual difference")+ scale_fill_manual(values=c("#453781ff","#238A8DFF","#55C667FF","#FDE725FF","grey"))+
  ylim(c(0,10000))

output_slope %>% filter(TypeNAT%in%c("Env","Dist","Dis_Env","Dis_Dist")) %>%group_by(TypeNAT,Year) %>% 
  summarise(Slope_m=mean(Slope)) %>% 
  ggplot()+geom_point(aes(x=Year,y=Slope_m,colour=TypeNAT),size=4)+scale_colour_manual(values=c("#453781ff","#238A8DFF","#55C667FF","#FDE725FF","grey"))+
  geom_line(aes(x=Year,y=Slope_m,colour=TypeNAT),size=0.5)+
  theme_classic()+theme(legend.position = "none")+ labs(y="LM Residuals")
  


output_slope%>% filter(TypeNAT%in%c("Env","Dist")) %>% #filteSlopeoutput_slope%>% filter(TypeNAT%in%c("Env","Dist")) %>% #filter(Year==2010) %>% 
  ggplot()+ geom_boxplot(aes(y=Slope,colour=as.factor(Year)))+
  facet_wrap(TypeNAT~.,ncol = 1)+coord_flip()


output_slope%>% filter(Year==2009) %>% 
ggplot() +
  aes(x = TypeNAT, y = slope, color = TypeNAT) +
  geom_jitter() +
  theme(legend.position = "none")
# 1st method:

model <- oneway.test(m ~ TypeNAT,
            data = output_slope,
            var.equal = TRUE # assuming equal variances
)
model <- aov(output_slope_riv$m~output_slope_riv$TypeNAT)

TukeyHSD(model)

#plotejant els residus de lm
data <- full_output_riv %>%
  pivot_longer(cols = 4:11) %>% filter(name=="mean_grStre") %>%
  group_by(TypeNAT,n_sites, Year) %>% 
  summarise(Mean_val=mean(value),sd_val=sd(value))
data2 <- data %>% subset(TypeNAT=="Dis_Env")
data2 <- data2 %>% subset(Year==2015)
fit <- lm(Mean_val~(as.numeric(n_sites)), data=data2) # fit the model
data2$predicted <- Model_nats   # Save the predicted values
data2$residuals <- slope 
ggplot(data2, aes(x = n_sites, y = Mean_val)) + geom_line(color="#55C667FF", size=1)+
  geom_line(aes(x=n_sites,y=predicted),colour="#DCE319FF", size=1)+
  #geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +     # regression line  
  geom_segment(aes(xend = n_sites, yend = predicted), colour="grey50",alpha = 1) + # draw line from point to line
  #geom_point(aes(color = abs(data2$residuals), size = abs(data2$residuals))) +  # size of the points
  #scale_color_continuous(low = "green", high = "red") +             # colour of the points mapped to residual size - green smaller, red larger
  guides(color = FALSE, size = FALSE) +                             # Size legend removed
  geom_point(aes(y = predicted), shape = 1) +
  theme_classic()+labs(x="Type of addition", y="Mean strength")
#_____________________________
