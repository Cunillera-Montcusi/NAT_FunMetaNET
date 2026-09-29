

library(tidyverse)

tEst <- Final_Output$lake$Out_NAT_Trait %>% 
  mutate(Trait_Cat = 
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
             
           )) # End of the paragraph

unique(tEst$Trait_Cat)


