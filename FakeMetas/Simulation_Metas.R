

source("FakeMetas/H2020_Lattice_expKernel_Jenv_TempMeta.R")

library(tidyverse); library(cooccur)

Lakes_coord <- readxl::read_excel("data/Lakes/lake_coord_surface.xlsx")

J.freshwater<-rep(200,nrow(Lakes_coord)) # J is the size of each community. Constant in this case
id_NOmodule <- rep(1,nrow(Lakes_coord)) # Modules just mean if we want some sites to belong to the same module. 
pool_200 <- rep(1,210) # Distribution of the species pool #rlnorm(n = 200,5,1) 
Meta_t0 <- matrix(nrow = length(pool_200), ncol =nrow(Lakes_coord), 1) #Previous Metacommunity (for considering time relevance)

Dist_Matr <- geosphere::distm(Lakes_coord[,2:3])/1000

# Filter assignation and ellaboration
# We generate a filter that benefits some packs of 15 species and punishes the others. 
# A filter that benefits is a "high" value -- 0.99 #Never write 1!
# A filter that punishes is a "low" value -- 0.1 
FilterS_to_Assign <- list()
for (filter_sp in seq(from=1,to=210,by=15)) {
  filter_Sp_end <- filter_sp+14 # Set positions to sequence from begining to end
  Filter<- rep(0.1,210) # Create a vector with "low" filters
  Filter[filter_sp:filter_Sp_end] <- 0.99 # Add to the "good" species the good filters
  
  Pos_To_Store <-which(seq(from=1,to=210,by=15)==filter_sp) # Look for which position you are in the loop
  FilterS_to_Assign[[Pos_To_Store]] <- Filter # Store the filter values for each species
}

# We create a "non-filter" scenario where everywhere and everyone is 0.99 
Filter_Scen <- matrix(ncol=nrow(Lakes_coord),nrow =length(rep(1,210)) ,data=0.99)
Filter_NO_Filter <- Filter_Scen # We store this filter for the "NO filter scenario"
# Here we select the filters that we want that can be selected randomly or manually to set specific 
## environmental selection (e.g. only 1 type of filter is present everywhere for example... or whatevere)
Rand_Filt <- sample(rep(1:14,56),size=56,replace = F) # We randomly select different tipes of filters to assign 

# Based onthe previous selection we assign the corresponding filters in a table that will be later used 
## in the coalescent model to set environmental impact. 
for (Assign_Filter in 1:length(Rand_Filt)) {
  Filter_Scen[,Assign_Filter] <- FilterS_to_Assign[[Rand_Filt[Assign_Filter]]]
}

#Traits matrix for later. We just create a trait matrix that will simulate changes in the 
## traits of our species and that will correspond to the filters that we will be using. 
Trait_Matrix <- matrix(nrow = 210,ncol=118,data = 0)
for (Sp_Tr in seq(1,210,3)) {
  Sp_Tr_End <- Sp_Tr+2
  Col_Beg <- Sp_Tr_End/3
  Col_End <- Col_Beg+24
  Trait_Matrix[Sp_Tr:Sp_Tr_End,Col_Beg:Col_End] <- 1
}
Trait_Matrix <- as.data.frame(Trait_Matrix) %>% mutate(Sp_Name=colnames(as.data.frame(t(Trait_Matrix))),.before=V1)
colnames(Trait_Matrix) <- c("Sp_Name",paste("Tra",seq(1:(ncol(Trait_Matrix)-1)),sep = "_"))

Scenarios <- c("Sc_Spa")#,"Sc_Env","Sc_SpaEnv")
Community_Scenarios <- list()

# D50 corresponding to the distance at which probability is 50%
dispersal_test <- c(10)#,4000,300)
Applied_Filter <- list(Filter_NO_Filter)#,Filter_Scen,Filter_Scen)

for (Scenari in 1:length(Scenarios)) {
a <- NULL # We create an output object for each iteration
b <- list()
for (it in 1:20) { # We repeat 10 times the same process
  output <- H2020_Coalescent.and.lottery.exp.Kernel.J_TempMtcom_tempIT(
    Meta.pool = pool_200, # Species pool
    m.pool = 0.001, # Regional dispersal which is always constant 
    Js = J.freshwater, # Size of the communities (AKA: number of individuals/population contained in each community)
    id.module = id_NOmodule, # id of modules if there are some - NOT used for us
    filter.env = Applied_Filter[[Scenari]], # Pollution scenarios (created at 2. Pollution assignation.R)
    M.dist =Dist_Matr, # Distance matrix which corresponds to the STconmat (created at 1. OCnet - STconmat.R)
    D50 = dispersal_test[Scenari], # Dispersal distance scenario 
    m.max = 1, # Maximum migration
    tempo_imp = 0, # Relvance of "temporal" effect
    temp_Metacom = Meta_t0, # Metacommunity at time 0 (all species are equally favored)
    temp_it = 0, # Number of temporal iterations
    id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=pool_200, # If there are some communit. that should be fixed
    Lottery=F, 
    it=300, 
    prop.dead.by.it=0.05, # Lottery parameters, nº iterations and proportion of dead organisms  
    id.obs=1:nrow(Lakes_coord)) # Information if we would like to keep specific results only
  a <- rbind(a,output[[1]])
  b[[it]] <- output[[2]]
}# it

Out_Community <- matrix(nrow=nrow(b[[1]]),ncol=ncol(b[[1]]),data=0)
for (iterat in 1:length(b)) {Out_Community <- Out_Community+b[[iterat]]}
Out_Community <- Out_Community/length(b)

Community_Scenarios[[Scenari]] <- Out_Community
}# Scenari end


# Here there will be a loop for each scenario
### THIS IS NORMAL NATS 
Out_Community<- Community_Scenarios[[1]]

### FIRST LOOP - Scenarios 
Out_Community <- as.data.frame(t(Out_Community)) %>% mutate(Site_ID=1:nrow(.),.before=V1) %>% 
                 pivot_longer(cols = 2:ncol(.))
Out_Community<- Out_Community%>%filter(value>0)
TypeNAT="Rand"
fake_output <- data.frame()
fake_output_traits_str <- data.frame()
for (ind_year in 1:1) {
cat("We are at Scenario", Scenarios[ind_year],"__________________________________________________","\n")
  
### SECOND LOOP - Number of randomly selected lakes
for (iteration in 1:10) {
cat("We are at iteration", iteration,"__________________________________________________","\n")
  
  
  if(TypeNAT=="Rand"){
  # We randomly select lake names
  select_lakes <- sample(seq(1:ncol(b[[1]])), # The vector we want to select things from 
                         size =length(unique(Out_Community$Site_ID)), # The number of elements that we want to select
                         replace = F) # If we can repeat or not  
  }  

  if(TypeNAT=="Env"){
  Chosen_Beg <- sample(seq(1:ncol(b[[1]])), size =1,  replace = F)
  Chosen_OneS <- Chosen_Beg
  Chosen_One<- Rand_Filt[Chosen_Beg]
  Sel_Neig <- c()
  Sel_Neig[1] <- Chosen_Beg
  for (Chosing_Lakes in 2:ncol(b[[1]])) {
  Comm_to_Sel <- seq(1:ncol(b[[1]]))[-Chosen_OneS]
  Filt_Diff <- Rand_Filt[-Chosen_OneS]
  Sel_Neig[Chosing_Lakes] <- Comm_to_Sel[order(abs(Chosen_One-Rand_Filt[-Chosen_OneS]))[1]]
  Chosen_OneS <- c(Sel_Neig)
  Chosen_One<- mean(Rand_Filt[Chosen_OneS])
  }
  select_lakes <- Chosen_OneS
  }
  
  if(TypeNAT=="Dist"){
    lake_geo <- Lakes_coord[,1:3]
    lake_to_start <- sample(lake_geo$Lake,1)

    Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
    Llista_llacs[1] <- Lakes_coord$Lake[which(Lakes_coord$Lake==lake_to_start)]#posem el 1r llac al vector
    lake_geo[which(Lakes_coord$Lake==lake_to_start),1] <- "Chosen_One"
    #Chosen_One <- lake_to_start
    for (Chosing_Lakes in 2:nrow(Lakes_coord)) {
      xy <- lake_geo[,2:3]  
      Dist_Matr <- geosphere::distm(xy)
      Chosen_OneS <- which(lake_geo$Lake=="Chosen_One") # Posició del Orig lake
      Chosen_col <- Dist_Matr[,Chosen_OneS]#aillem columna de l'escollit
      Dist_Diff <-Chosen_col[-Chosen_OneS]#treiem la posicio que es 0
      Dist_Neigh <- which(Dist_Diff==min(Dist_Diff)) # Trobem la distànica minima entre Original i altres
      Closer <- Dist_Neigh[1]
      d <- subset(lake_geo, Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                      lake_geo$Lake[Closer]))#filtrem els dos llacs més propers
      new_lake <- data.frame(Lake="Chosen_One",
                             summarise(d,Lon=mean(Lon),Lat=mean(Lat)))
      Llista_llacs[Chosing_Lakes] <- lake_geo$Lake[Closer]
      
      lake_geo <- lake_geo %>% filter(!Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                                 lake_geo$Lake[Closer])) %>% 
        bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    }
select_lakes <- c()
for (Pos_Lak in 1:length(Lakes_coord$Lake)) {
select_lakes[Pos_Lak] <- which(Lakes_coord$Lake==Llista_llacs[Pos_Lak])
}} # If Spa ending 

  if(TypeNAT=="Dis_Env"){
    Chosen_Beg <- sample(seq(1:ncol(b[[1]])), size =1,  replace = F)
    Chosen_OneS <- Chosen_Beg
    Chosen_One<- Rand_Filt[Chosen_Beg]
    Sel_Neig <- c()
    Sel_Neig[1] <- Chosen_Beg
    for (Chosing_Lakes in 2:ncol(b[[1]])) {
      Comm_to_Sel <- seq(1:ncol(b[[1]]))[-Chosen_OneS]
      Filt_Diff <- Rand_Filt[-Chosen_OneS]
      Sel_Neig[Chosing_Lakes] <- Comm_to_Sel[order(abs(Chosen_One-Rand_Filt[-Chosen_OneS]))[length(Filt_Diff)]]
      Chosen_OneS <- c(Sel_Neig)
      Chosen_One<- mean(Rand_Filt[Chosen_OneS])
    }
    select_lakes <- Chosen_OneS
  }

  if(TypeNAT=="Dis_Dist"){
    lake_geo <- Lakes_coord[,1:3]
    lake_to_start <- sample(lake_geo$Lake,1)
    
    Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
    Llista_llacs[1] <- Lakes_coord$Lake[which(Lakes_coord$Lake==lake_to_start)]#posem el 1r llac al vector
    lake_geo[which(Lakes_coord$Lake==lake_to_start),1] <- "Chosen_One"
    #Chosen_One <- lake_to_start
    for (Chosing_Lakes in 2:nrow(Lakes_coord)) {
      xy <- lake_geo[,2:3]  
      Dist_Matr <- geosphere::distm(xy)
      Chosen_OneS <- which(lake_geo$Lake=="Chosen_One") # Posició del Orig lake
      Chosen_col <- Dist_Matr[,Chosen_OneS]#aillem columna de l'escollit
      Dist_Diff <-Chosen_col[-Chosen_OneS]#treiem la posicio que es 0
      Dist_Neigh <- which(Dist_Diff==max(Dist_Diff)) # Trobem la distànica minima entre Original i altres
      Closer <- Dist_Neigh[1]
      d <- subset(lake_geo, Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                      lake_geo$Lake[Closer]))#filtrem els dos llacs més propers
      new_lake <- data.frame(Lake="Chosen_One",
                             summarise(d,Lon=mean(Lon),Lat=mean(Lat)))
      Llista_llacs[Chosing_Lakes] <- lake_geo$Lake[Closer]
      
      lake_geo <- lake_geo %>% filter(!Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                                 lake_geo$Lake[Closer])) %>% 
        bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    }
    select_lakes <- c()
    for (Pos_Lak in 1:length(Lakes_coord$Lake)) {
      select_lakes[Pos_Lak] <- which(Lakes_coord$Lake==Llista_llacs[Pos_Lak])
    }}
for (selected_lakes in c(1,3,5,6,8,10,21,56)) {
cat("We have seleted", selected_lakes,"lakes","__________________________________________________","\n")
### THIRD LOOP - We will repeat the same thing  several times  
real_select_lakes <- select_lakes[1:selected_lakes]
macros_lakes_list_temp_temp <- Out_Community %>% filter(Site_ID%in%real_select_lakes)

# Check presence
traits_ind <- macros_lakes_list_temp_temp %>% ungroup() %>% # just in case 
        left_join(Trait_Matrix, by=c("name"="Sp_Name"),multiple ="all") %>% # join traits with the genus of the traits database
        na.omit() 
sp_rich <- nrow(traits_ind)
traits_ind_abund <- data.frame()
for (Row_ID in 1:nrow(traits_ind)) {
  Out_sel_row <- data.frame()
  for (Row_Abun in 1:ceiling(traits_ind$value[Row_ID])) {
    Selec_Row <- traits_ind[Row_ID,] 
    Out_sel_row <- bind_rows(Out_sel_row,Selec_Row)
  }
  traits_ind_abund <- bind_rows(traits_ind_abund,Out_sel_row)
}  

traits_ind_abund <-traits_ind_abund %>%dplyr::select(c(4:ncol(.))) %>% na.omit()

source("function_to_NATs.R")
# We have our database! We can calculate the network! 
NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind_abund)

edge_dens <- igraph::edge_density(NATs_Output$Graph)
mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
gr_Stre <- igraph::graph.strength(NATs_Output$Graph)
transit_W <- mean(igraph::transitivity(NATs_Output$Graph, type = "barrat"),na.rm=TRUE)
transit <- igraph::transitivity(NATs_Output$Graph, type = "global")


temp_output <- data.frame(#"Type_NATS"=Type_of_NATs,
                          "Scenario"=Scenarios[ind_year],
                          "n_sites"=selected_lakes,
                          "iter"=iteration,
                          "richness"=sp_rich,
                          edge_dens,
                          mean_grStre,
                          transit_W,
                          transit)
fake_output <- bind_rows(fake_output,temp_output)
temp_output_traits_str <- data.frame(#"Type_NATS"=Type_of_NATs,
                                "Scenario"=Scenarios[ind_year],
                                 "n_sites"=selected_lakes,
                                 "iter"=iteration,
                                  gr_Stre, Tra=colnames(Trait_Matrix)[2:ncol(Trait_Matrix)])
fake_output_traits_str <- bind_rows(fake_output_traits_str,temp_output_traits_str)
    }# End of selected_lakes
  }# End of iteration
}# End of ind_year
save(fake_output,file = "fake_NATs_dist_prova1_disp10km.RData")
save(fake_output_traits_str,file = "fake_NATs_dist_prova1_disp10km.RData")

temp_output %>% 
  group_by(Scenario,iter,n_sites) %>% 
  summarise(edge_dens=mean(edge_dens),mean_grStre=mean(mean_grStre)) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=edge_dens))


load("fake_NATs_rand_prova1_disp300km.RData")
rand_lakes <- fake_output
load("fake_NATs_env_prova1_disp300km.RData")
env_lakes <- fake_output
load("fake_NATs_dist_prova1_disp300km.RData")
dist_lakes <- fake_output
load("fake_NATs_Disenv_prova1_disp300km.RData")
DisEnv_lakes <- fake_output
load("fake_NATs_Disdist_prova1_disp300km.RData")
DisDist_lakes <- fake_output

load("fake_NATs_rand_prova1_disp10km.RData")
rand_lakes <- fake_output
load("fake_NATs_dist_prova1_disp10km.RData")
dist_lakes <- fake_output
load("fake_NATs_Disdist_prova1_disp1000km.RData")
DisDist_lakes <- fake_output

full_fake_output <- bind_rows(
  rand_lakes%>% mutate(TypeNAT="Rand"),
  #env_lakes%>% mutate(TypeNAT="Env"),
  dist_lakes%>% mutate(TypeNAT="Dist"))
  #DisEnv_lakes%>% mutate(TypeNAT="Dis_Env"),
  #DisDist_lakes%>% mutate(TypeNAT="Dis_Dist"))

rand_lakes$mean_grStre-dist_lakes$mean_grStre

unique(full_fake_output$TypeNAT)
colnames(full_fake_output)
full_fake_output %>%
  pivot_longer(cols = 4:8) %>% filter(name=="mean_grStre") %>%
  filter(TypeNAT%in%c("Rand", "Dist")) %>% 
  group_by(TypeNAT,n_sites) %>% 
  summarise(Mean_val=mean(value),sd_val=sd(value)) %>% 
  ggplot(aes(y=Mean_val, x=n_sites, group=TypeNAT))+ 
  geom_line(size=1, aes(linetype=TypeNAT, colour= TypeNAT))+ scale_linetype_manual(values=c("solid","solid","solid","solid", "dashed"))+
  scale_x_continuous(breaks = seq(10, 50, by = 10)) +scale_colour_manual(values=c("#481567FF","#33638DFF","#238A8DFF","#55C667FF","grey"))+
  #geom_errorbar(aes(ymin=Mean_val-sd_val,ymax=Mean_val+sd_val))+
  theme_classic()
  #+facet_wrap(.~Year) + labs(y="Mean strength", x= "Number of sites")+theme(legend.position = "none")


load("fake_NATs_rand_prova1_disp10km.RData")
rand_traits <- fake_output_traits_str
load("fake_NATs_dist_prova1_disp10km.RData")
dist_traits <- fake_output_traits_str

rand_traits %>% rename("Ran_gr_Stre"="gr_Stre") %>% 
  left_join(dist_traits %>% rename("Dist_gr_Stre"="gr_Stre"),
            by=c("Scenario","n_sites","iter","Tra")) %>% 
  mutate(Diff=Ran_gr_Stre-Dist_gr_Stre) %>% 
  group_by(Scenario,n_sites,iter) %>% 
  summarise(mean_Diff=mean(Diff)) %>% 
  ggplot()+
  geom_point(aes(x=n_sites,y=mean_Diff))

bind_cols(rand_traits,dist_traits)



