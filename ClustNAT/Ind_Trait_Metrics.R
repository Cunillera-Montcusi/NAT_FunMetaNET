
library(tidyverse);library(viridis)
library(drc);library(nlme);library(statforbiology)

load("ClustNAT/ClusterNATs.RData")

Merged_Traits <- bind_rows(
  Final_Output$lake$Out_NAT_Trait %>% mutate(Ecosy="Lake"),
  Final_Output$river$Out_NAT_Trait %>% mutate(Ecosy="River")
  ) %>% mutate(Trait_Cat = 
           case_when(
             # Size category 
             str_detect(Tra ,"X..0.25.cm") ~ "Size",
             str_detect(Tra ,"X..0.25.0.5.cm") ~ "Size",
             str_detect(Tra ,"X..0.5.1.cm") ~ "Size",
             str_detect(Tra ,"X..1.2.cm") ~ "Size",
             str_detect(Tra ,"X..2.4.cm") ~ "Size",
             str_detect(Tra ,"X..4.8.cm") ~ "Size",
             str_detect(Tra ,"X..8.cm") ~ "Size",
             #Lifespan 
             str_detect(Tra ,"X.1.year") ~ "Life_Span",
             str_detect(Tra ,"X..1.year") ~ "Life_Span",
             # SOMETHING THAT I DO NOT KNOW!!!
             str_detect(Tra ,"X..1") ~ "Reproductive_cycles", # WARNING !!!!
             str_detect(Tra ,"X1") ~ "Reproductive_cycle",# WARNING !!!!
             str_detect(Tra ,"X..1.1") ~ "Reproductive_cycle",# WARNING !!!!
             # Aquatic lifestage 
             str_detect(Tra ,"egg") ~ "Aquatic_stage",
             str_detect(Tra ,"larva") ~ "Aquatic_stage",
             str_detect(Tra ,"nymph") ~ "Aquatic_stage",
             str_detect(Tra ,"adult") ~ "Aquatic_stage",
             # Dispersal strategy
             str_detect(Tra ,"aquatic.passive") ~ "Dispersal",
             str_detect(Tra ,"aquatic.active") ~ "Dispersal",
             str_detect(Tra ,"aerial.passive") ~ "Dispersal",
             str_detect(Tra ,"aerial.active") ~ "Dispersal",
             # Reproduction
             str_detect(Tra ,"ovoviviparity") ~ "Reproduction",
             str_detect(Tra ,"isolated.eggs..free") ~ "Reproduction",
             str_detect(Tra ,"isolated.eggs..cemented") ~ "Reproduction",
             str_detect(Tra ,"clutches..cemented.or.fixed") ~ "Reproduction",
             str_detect(Tra ,"clutches..free") ~ "Reproduction",
             str_detect(Tra ,"clutches..in.vegetation") ~ "Reproduction",
             str_detect(Tra ,"clutches..terrestrial") ~ "Reproduction",
             str_detect(Tra ,"asexual.reproduction") ~ "Reproduction",
             # Resistance 
             str_detect(Tra ,"eggs..statoblasts") ~ "Resistance",
             str_detect(Tra ,"cocoons") ~ "Resistance",
             str_detect(Tra ,"housings.against.desiccation") ~ "Resistance",
             str_detect(Tra ,"diapause.or.dormancy") ~ "Resistance",
             str_detect(Tra ,"none") ~ "Resistance",
             # Respiration
             str_detect(Tra ,"interstitial") ~ "Respiration",
             str_detect(Tra ,"tegument") ~ "Respiration",
             str_detect(Tra ,"gill") ~ "Respiration",
             str_detect(Tra ,"plastron") ~ "Respiration",
             str_detect(Tra ,"spiracle") ~ "Respiration",
             str_detect(Tra ,"hydrostatic.vesicle") ~ "Respiration",
             # Movement capacity
             str_detect(Tra ,"flier") ~ "Movement",
             str_detect(Tra ,"surface.swimmer") ~ "Movement",
             str_detect(Tra ,"full.water.swimmer") ~ "Movement",
             str_detect(Tra ,"crawler") ~ "Movement",
             str_detect(Tra ,"burrower") ~ "Movement",
             str_detect(Tra ,"temporarily.attached") ~ "Movement",
             str_detect(Tra ,"permanently.attached") ~ "Movement",
             # Food source
             str_detect(Tra ,"microorganisms") ~ "Food_source",
             str_detect(Tra ,"detritus...1mm") ~ "Food_source",
             str_detect(Tra ,"dead.plant....1mm") ~ "Food_source",
             str_detect(Tra ,"living.microphytes") ~ "Food_source",
             str_detect(Tra ,"living.macrophytes") ~ "Food_source",
             str_detect(Tra ,"dead.animal....1mm") ~ "Food_source",
             str_detect(Tra ,"living.microinvertebrates") ~ "Food_source",
             str_detect(Tra ,"living.macroinvertebrates") ~ "Food_source",
             str_detect(Tra ,"vertebrates") ~ "Food_source",
             # Way of feeding
             str_detect(Tra ,"absorber") ~ "Way_of_feeding",
             str_detect(Tra ,"deposit.feeder") ~ "Way_of_feeding",
             str_detect(Tra ,"shredder") ~ "Way_of_feeding",
             str_detect(Tra ,"scraper") ~ "Way_of_feeding",
             str_detect(Tra ,"filter.feeder") ~ "Way_of_feeding",
             str_detect(Tra ,"piercer") ~ "Way_of_feeding",
             str_detect(Tra ,"predator") ~ "Way_of_feeding",
             str_detect(Tra ,"parasite") ~ "Way_of_feeding",
             # Habitat
             str_detect(Tra ,"river.channel") ~ "Habitat",
             str_detect(Tra ,"banks..connected.side.arms") ~ "Habitat",
             str_detect(Tra ,"ponds..pools..disconnected.side.arms") ~ "Habitat",
             str_detect(Tra ,"marshes..peat.bogs") ~ "Habitat",
             str_detect(Tra ,"temporary.waters") ~ "Habitat",
             str_detect(Tra ,"lakes") ~ "Habitat",
             str_detect(Tra ,"groundwaters") ~ "Habitat",
             # SOMETHING THAT I DO NOT KNOW!!!
             str_detect(Tra ,"crenon") ~ "Longitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"epirithron") ~ "Longitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"metarithron") ~ "Longitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"hyporithron") ~ "Longitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"epipotamon") ~ "Longitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"metapotamon") ~ "Longitudinal_distribution", # WARNING !!!!
             # SOMETHING THAT I DO NOT KNOW!!!
             str_detect(Tra ,"estuary") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"outside.river.system") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"lowlands") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"piedmont.level") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"alpine.level") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"X2...Pyrenees") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"X4...Alps") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"X8...Vosges..Jura..Massif.Central") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"X13a...lowlands..oceanic.") ~ "Latitudinal_distribution", # WARNING !!!!
             str_detect(Tra ,"X13b...lowlands..mediterranean.") ~ "Latitudinal_distribution", # WARNING !!!!
             # SOMETHING THAT I DO NOT KNOW!!! SUbstrate? 
             str_detect(Tra ,"flags.boulders.cobbles.pebbles") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"gravel") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"sand") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"silt") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"macrophytes") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"microphytes") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"twigs.roots") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"organic.detritus.litter") ~ "Substrate_preference", # WARNING !!!!
             str_detect(Tra ,"mud") ~ "Substrate_preference", # WARNING !!!!
             # SOMETHING THAT I DO NOT KNOW!!! 
             str_detect(Tra ,"null") ~ "Velocity", # WARNING !!!!
             str_detect(Tra ,"slow") ~ "Velocity", # WARNING !!!!
             str_detect(Tra ,"medium") ~ "Velocity", # WARNING !!!!
             str_detect(Tra ,"fast") ~ "Velocity", # WARNING !!!!
             # Habitat nutrients 
             str_detect(Tra ,"oligotrophic") ~ "Habitat_Quality",
             str_detect(Tra ,"mesotrophic") ~ "Habitat_Quality",
             str_detect(Tra ,"eutrophic") ~ "Habitat_Quality",
             # Habitat preference
             str_detect(Tra ,"fresh.water") ~ "Habitat_preference",  
             str_detect(Tra ,"brackish.water") ~ "Habitat_preference",  
             str_detect(Tra ,"psychrophilic") ~ "Habitat_preference",  
             str_detect(Tra ,"thermophilic") ~ "Habitat_preference",  
             str_detect(Tra ,"eurythermic") ~ "Habitat_preference",  
             str_detect(Tra ,"xenosaprobic") ~ "Habitat_preference",  
             str_detect(Tra ,"oligosaprobic") ~ "Habitat_preference",  
             str_detect(Tra ,"b.mesosaprobic") ~ "Habitat_preference",  
             str_detect(Tra ,"a.mesosaprobic") ~ "Habitat_preference",  
             str_detect(Tra ,"polysaprobic") ~ "Habitat_preference",  
             # SOMETHING THAT I DO NOT KNOW!!! SUbstrate? 
             str_detect(Tra ,"X..4") ~ "pH", # WARNING !!!!
             str_detect(Tra ,"X..4.4.5") ~ "pH", # WARNING !!!!
             str_detect(Tra ,"X..4.5.5") ~ "pH", # WARNING !!!!
             str_detect(Tra ,"X..5.5.5") ~ "pH", # WARNING !!!!
             str_detect(Tra ,"X..5.5.6") ~ "pH", # WARNING !!!!
             str_detect(Tra ,"X..6") ~ "pH" # WARNING !!!!
           ))
  
 
Big_Drop <- Merged_Traits %>% 
  filter(n_sites%in%c(1,39,56)) %>% 
  filter(Type_NATS=="Random") %>% 
  group_by(Ecosy,Type_NATS,n_sites,Trait_Cat,Tra) %>% 
  summarise(Mean_Gr=mean(log(gr_Stre+1)),
            Mean_Bet=mean(log(betw_Stre+1))
            ) %>%
  dplyr::select(-c(Mean_Gr)) %>% 
  mutate(n_sites=ifelse(n_sites<10,"One","All")) %>% 
  pivot_wider(names_from = n_sites, values_from = Mean_Bet) %>% 
  mutate(penis=All-One) %>%
  mutate(Drop=ifelse(penis<(-2),"Yes","No")) %>% 
  dplyr::select(c(Ecosy,Type_NATS,Trait_Cat,Tra,Drop))

write.csv2(x = Big_Drop %>% filter(Drop=="Yes"), file = "C:/Users/David CM/Desktop/Traits_Drop.csv")
  

Merged_Traits %>% 
  #filter(n_sites%in%c(1)) %>% 
  filter(Type_NATS=="Random") %>% 
  filter(Tra%in%c("aquatic.passive","aquatic.active","aerial.passive","aerial.active")) %>% 
  group_by(Ecosy,Type_NATS,n_sites,Tra) %>% 
  summarise(Mean_Gr=mean(log(gr_Stre+1)),
            Mean_Bet=mean(log(betw_Stre+1)),
            SD_Gr=sd(log(gr_Stre+1)),
            SD_Bet=sd(log(betw_Stre+1))) %>% 
  ggplot(aes(x=Mean_Gr,y=Mean_Bet, shape=Tra,colour=as.factor(n_sites)))+
  geom_linerange(aes(x=Mean_Gr, xmin=(Mean_Gr-SD_Gr),xmax=(Mean_Gr+SD_Gr)))+
  geom_linerange(aes(y=Mean_Bet,ymin=(Mean_Bet-SD_Bet),ymax=(Mean_Bet+SD_Bet)))+
  geom_point(size=4)+
  scale_color_viridis(discrete = T,direction = -1)+
  facet_wrap(Ecosy~Type_NATS)+
  theme_classic()

gridExtra::grid.arrange(
Merged_Traits %>% 
  filter(Type_NATS=="Random") %>% 
  filter(n_sites%in%c(1,39,56)) %>% 
  group_by(Ecosy,Type_NATS,n_sites,Trait_Cat,Tra) %>% 
  summarise(Mean_Gr=mean(log(gr_Stre+1)),
            Mean_Bet=mean(log(betw_Stre+1)),
            SD_Gr=sd(log(gr_Stre+1)),
            SD_Bet=sd(log(betw_Stre+1))) %>% 
  mutate(n_sites=ifelse(n_sites<10,"One","All")) %>% 
  left_join(Big_Drop,by=c("Ecosy","Type_NATS","Trait_Cat","Tra")) %>% 
  ggplot()+
  geom_vline(xintercept = 5,size=2,colour="grey50")+
  geom_hline(yintercept = 40,size=2,colour="grey50")+
  geom_linerange(aes(y=Mean_Bet, x=Mean_Gr, xmin=(Mean_Gr-SD_Gr),
                     xmax=(Mean_Gr+SD_Gr),colour=as.factor(n_sites)),alpha=0.2)+
  geom_linerange(aes(x=Mean_Gr, y=Mean_Bet,ymin=(Mean_Bet-SD_Bet),
                     ymax=(Mean_Bet+SD_Bet),colour=as.factor(n_sites)),alpha=0.2)+
  geom_line(aes(x=Mean_Gr,y=Mean_Bet,group=Tra,linewidth=Drop,alpha=Drop))+
  geom_point(aes(x=Mean_Gr,y=Mean_Bet,color=as.factor(n_sites),shape=as.factor(n_sites)),
             size=4,alpha=0.5)+
  #geom_smooth(aes(group=Tra),linewidth=0.2,method="lm",colour="black",se=F)+
  scale_color_manual(values = c(viridis(n = 2,option = "C",direction = 1)))+
  scale_linewidth_manual(values=c(0.5,2))+
  scale_alpha_manual(values=c(0.2,0.6))+
  scale_shape_manual(values = c(17,16))+
  facet_wrap(.~Ecosy, ncol=2)+
  labs(title="Cross-traits network role",
       y="Trait Betweenness (log-scale)",x="Trait Strength (log-scale)")+
  theme_classic()+
  theme(legend.position = "none"),

# `Plot per traits specivis`
Merged_Traits %>% 
  filter(Type_NATS=="Random") %>% 
  filter(Trait_Cat%in%c("Dispersal","Food_source","Reproduction")) %>% 
  filter(n_sites%in%c(1,39,56)) %>% 
  group_by(Ecosy,Type_NATS,n_sites,Trait_Cat,Tra) %>% 
  summarise(Mean_Gr=mean(log(gr_Stre+1)),
            Mean_Bet=mean(log(betw_Stre+1)),
            SD_Gr=sd(log(gr_Stre+1)),
            SD_Bet=sd(log(betw_Stre+1))) %>%
  mutate(n_sites=ifelse(n_sites<10,"One","All")) %>% 
  left_join(Big_Drop,by=c("Ecosy","Type_NATS","Trait_Cat","Tra")) %>% 
  ggplot()+
  geom_linerange(aes(y=Mean_Bet, x=Mean_Gr, xmin=(Mean_Gr-SD_Gr),
                     xmax=(Mean_Gr+SD_Gr),colour=as.factor(n_sites)))+
  geom_linerange(aes(x=Mean_Gr, y=Mean_Bet,ymin=(Mean_Bet-SD_Bet),
                     ymax=(Mean_Bet+SD_Bet),colour=as.factor(n_sites)))+
  geom_line(aes(x=Mean_Gr,y=Mean_Bet,group=Tra,linewidth=Drop,alpha=Drop))+
  geom_point(aes(x=Mean_Gr,y=Mean_Bet,colour=as.factor(n_sites),shape=as.factor(n_sites)),size=4,alpha=0.5)+
  #geom_smooth(aes(group=Tra),linewidth=0.2,method="lm",colour="black",se=F)+
  scale_color_manual(values = c(viridis(n = 2,option = "C",direction = 1)))+
  scale_linewidth_manual(values=c(0.5,2))+
  scale_alpha_manual(values=c(0.2,0.6))+
  scale_shape_manual(values = c(17,16))+
  facet_wrap(Trait_Cat~Ecosy, ncol=2)+
  guides(shape="none")+
  labs(y="Trait Betweenness (log-scale)",x="Trait Strength (log-scale)",
       colour="Nº added communities")+
  theme_classic(),
ncol=2)

Merge_Traits <- bind_rows(
  Final_Output$lake$Out_NAT_Trait %>% mutate(Ecosy="Lake"),
  Final_Output$river$Out_NAT_Trait %>% mutate(Ecosy="River")
) %>% 
  group_by(Ecosy,Type_NATS,n_sites,Tra) %>% 
  summarise(Mean_Gr=mean(log(gr_Stre+1)),
            Mean_Bet=mean(log(betw_Stre+1)),
            SD_Gr=sd(log(gr_Stre+1)),
            SD_Bet=sd(log(betw_Stre+1))) %>% 
  filter(Type_NATS=="Random")

Out <- data.frame()
for (Ecosystem in unique(Merge_Traits$Ecosy)) {
  Eco_Merge_Traits <- Merge_Traits %>% filter(Ecosy==Ecosystem)
  for (Traits in unique(Eco_Merge_Traits$Tra)) {
  Tra_Eco_Merge_Traits <- Eco_Merge_Traits %>% filter(Tra==Traits)
test_model <- lm(Mean_Bet~Mean_Gr,data =Tra_Eco_Merge_Traits)
Out_temp <- data.frame(
  "Ecosy"=Ecosystem,
  "Tra"=Traits,
  "Slope"=test_model$coefficients[2])
Out <- bind_rows(Out,Out_temp)
  }
}

Out %>%
  #filter(Slope<(-5)) %>% 
  filter(Tra%in%c("aquatic.passive","aquatic.active","aerial.passive","aerial.active",
                  "deposit.feeder","shredder",
                  "scraper","filter.feeder","piercer","predator","parasite")) %>% 
  #ggplot(aes(y=Tra,x=log(abs(Slope)+1),colour=Ecosy))+
  ggplot(aes(y=Tra,x=Slope,colour=Ecosy))+
  geom_vline(xintercept = 0)+
  geom_point()+
  theme(axis.text.y = element_text(size=5))

