library(tidyverse)
#install.packages("devtools")
#devtools::install_github("onofriAndreaPG/aomisc")
library(aomisc)
#calcul de pendents LLACS
load("NATs_lakes_rand_def.RData")
rand_lakes <- output
load("NATs_lakes_env_def.RData")
env_lakes <- output
load("NATs_lakes_Dist_def.RData")
dist_lakes <- output
load("NATs_lakes_DisEnv_def.RData")
DisEnv_lakes <- output
load("NATs_lakes_DisDist_def.RData")
DisDist_lakes <- output


full_output <- bind_rows(
  rand_lakes%>% mutate(TypeNAT="Rand"),
  env_lakes%>% mutate(TypeNAT="Env"),
  dist_lakes%>% mutate(TypeNAT="Dist"),
  DisEnv_lakes%>% mutate(TypeNAT="Dis_Env"),
  DisDist_lakes%>% mutate(TypeNAT="Dis_Dist"))

years <- unique(full_output$Year)
type <- unique(full_output$TypeNAT)
output_slope <- data.frame()

### FIRST LOOP - TypeNAT
for (ind_type in 1:length(type)) {
  # We create a "temporary" file filtered according to the TypeNat selected. 
  full_ouput_temp <- full_output %>% filter(TypeNAT==type[ind_type])
  # All years of that typeNAT
  year_type <- unique(full_ouput_temp$Year)

### 2ND LOOP - YEAR
  for (year in 1:length(year_type)) {
    
    random <- full_output %>% filter(TypeNAT=="Rand",Year==years[year]) %>% 
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
    
    output_slope <- rbind(output_slope,output_slope2)
  }
}


# Merge with richness
unique(full_output$TypeNAT)

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
