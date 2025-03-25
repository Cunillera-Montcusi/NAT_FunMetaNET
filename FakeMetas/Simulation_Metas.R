
# Beginning ####
library(tidyverse)# library(cooccur)
library(geosphere)
library(doParallel)
library(parallel)

#Lakes_coord <- readxl::read_excel("data/Lakes/lake_coord_surface.xlsx")
Lakes_coord <- readxl::read_excel("NAT/lake_coord_surface.xlsx")

Gamma_div <- 210 # Warning: A change in GAMMA diversity might impact some values linked to the filters (divided by 210)
J.freshwater<-rep(150,nrow(Lakes_coord)) # J is the size of each community. Constant in this case
id_NOmodule <- rep(1,nrow(Lakes_coord)) # Modules just mean if we want some sites to belong to the same module. 
pool_200 <- rep(1,Gamma_div) # Distribution of the species pool #rlnorm(n = 200,5,1) 
Meta_t0 <- matrix(nrow = length(pool_200), ncol =nrow(Lakes_coord), 1) #Previous Metacommunity (for considering time relevance)


# 1. Replicates ####
cl <- detectCores() #Number of cores in computer
registerDoParallel(cl)

Number_Of_Replicates <- 15
out <- foreach(Replicates=1:Number_Of_Replicates)%dopar%{
#for (Replicates in 1:6) {
library(tidyverse);library(geosphere)# Somehow we need to recharge the packages again.
source("NAT/H2020_Lattice_expKernel_Jenv_TempMeta_DispStr.R")
source("NAT/function_to_Cooccur.R")  

Community_Scenarios <- list()
Rand_Filt_list <- list()

# Scenari reffers to how many scenarios do we use. So far four: Env, Spa, Both, Null
# 2. Scenari ####
Scenarios <- c("Env","Spa", "Both", "Null")

for (Scenari in 1:length(Scenarios)) {
# 2.1. Filter assignation and ellaboration ####
# We generate a filter that benefits some packs of 15 species and punishes the others. 
# A filter that benefits is a "high" value -- 0.99 #Never write 1!
# A filter that punishes is a "low" value -- 0.1 
FilterS_to_Assign <- list()
for (filter_sp in seq(from=1,to=Gamma_div,by=15)) {
  filter_Sp_end <- filter_sp+14 # Set positions to sequence from begining to end
  Filter<- rep(0.2,Gamma_div) # Create a vector with "low" filters
  Filter[filter_sp:filter_Sp_end] <- 0.99 # Add to the "good" species the good filters
  
  Pos_To_Store <-which(seq(from=1,to=Gamma_div,by=15)==filter_sp) # Look for which position you are in the loop
  FilterS_to_Assign[[Pos_To_Store]] <- Filter # Store the filter values for each species
}

# 2.2. Traits matrix #### 
# We just create a trait matrix that will simulate changes in the 
## traits of our species and that will correspond to the filters that we will be using. 
Trait_Matrix <- matrix(nrow = Gamma_div,ncol=118,data = 0)
for (Sp_Tr in seq(1,Gamma_div,3)) {
  Sp_Tr_End <- Sp_Tr+2
  Col_Beg <- Sp_Tr_End/3
  Col_End <- Col_Beg+14
  Trait_Matrix[Sp_Tr:Sp_Tr_End,Col_Beg:Col_End] <- 1
}
Trait_Matrix <- as.data.frame(Trait_Matrix) %>% mutate(Sp_Name=colnames(as.data.frame(t(Trait_Matrix))),.before=V1)
colnames(Trait_Matrix) <- c("Sp_Name",paste("Tra",seq(1:(ncol(Trait_Matrix)-1)),sep = "_"))

# We assign and create different "trait profiles" based on the distributed traits. We sum different sections of traits 
# from species and create a table containing the sum of traits for those species. 
# Having a determined value will represent the belonging of that species to a "functional group" and will define 
# a link with a environmental filter (in the following lines)
Trait_Filter_Profile1 <- c();Trait_Filter_Profile2 <- c()
Trait_Filter_Profile3 <- c();Trait_Filter_Profile4 <- c()
Trait_Filter_Profile5 <- c();Trait_Filter_Profile6 <- c()
for (Spp_Trait in 1:nrow(Trait_Matrix)) {
Trait_Filter_Profile1[Spp_Trait] <- sum(Trait_Matrix[Spp_Trait,2:ceiling(118/6)])  
Trait_Filter_Profile2[Spp_Trait] <- sum(Trait_Matrix[Spp_Trait,ceiling(118/6):(ceiling(118/6)*2)])  
Trait_Filter_Profile3[Spp_Trait] <- sum(Trait_Matrix[Spp_Trait,(ceiling(118/6)*2):(ceiling(118/6)*3)])  
Trait_Filter_Profile4[Spp_Trait] <- sum(Trait_Matrix[Spp_Trait,(ceiling(118/6)*3):(ceiling(118/6)*4)])  
Trait_Filter_Profile5[Spp_Trait] <- sum(Trait_Matrix[Spp_Trait,(ceiling(118/6)*4):(ceiling(118/6)*5)])  
Trait_Filter_Profile6[Spp_Trait] <- sum(Trait_Matrix[Spp_Trait,(ceiling(118/6)*5):((118/6)*6)])  
}
Trait_Performance <- rbind(Trait_Filter_Profile1,Trait_Filter_Profile2,Trait_Filter_Profile3,Trait_Filter_Profile4,
          Trait_Filter_Profile5,Trait_Filter_Profile6)

# We create a "non-filter" scenario where everywhere and everyone is 0.99 
Filter_NO_Filter <- matrix(ncol=nrow(Lakes_coord),nrow =length(rep(1,Gamma_div)) ,data=0.99)

# Here we select the filters that we want that can be selected randomly or manually to set specific 
## environmental selection (e.g. only 1 type of filter is present everywhere for example... or whatever)
#Rand_Filt <- sample(rep(1:14,56),size=56,replace = F) # We randomly select different tipes of filters to assign 

# Based onthe previous selection we assign the corresponding filters in a table that will be later used 
## in the coalescent model to set environmental impact.
Filter_Scen <- matrix(ncol=nrow(Lakes_coord),nrow =length(rep(1,Gamma_div)) ,data=0.05)
Perc_Of_Filter_Assignment <- 0.95
Sites_With_Filter <- sample(1:ncol(Filter_Scen),size = ncol(Filter_Scen)*Perc_Of_Filter_Assignment,replace = F)
Order_Sites_With_Filter <- c(round(seq(from=1,to=length(Sites_With_Filter),by=length(Sites_With_Filter)/6)),length(Sites_With_Filter))
for (Assign_Filter in 1:nrow(Trait_Performance)) {
Spp_to_Assing <- which(Trait_Performance[Assign_Filter,]>5)
Sites_to_Filter <- Sites_With_Filter[Order_Sites_With_Filter[Assign_Filter]:Order_Sites_With_Filter[Assign_Filter+1]]
Filter_Scen[Spp_to_Assing,Sites_to_Filter] <- FilterS_to_Assign[[Assign_Filter]][Spp_to_Assing]
}
# Random filter is only summarising the filters in each lake. The name is arbitrary.
#for(Sites in 1:56){
#  print(length(which(Filter_Scen[,Sites]==0.99)))
#}
Rand_Filt <-apply(Filter_Scen,2,mean)

# We set different dispersal abilities for the species. We distributied them like this so 
# species with different trait configurations have different dispersal abilities
Disp_Str <- c(rep(1,Gamma_div/3),rep(2,Gamma_div/3),rep(3,Gamma_div/3))

# 2.3. Distance and dispersal ####
# Distance matrix corresponding to distances between lakes
Dist_True_Matr <- distm(Lakes_coord[,2:3])/1000
Dist_Matr_High <- Dist_True_Matr
Dist_Matr_Mid <- ifelse(Dist_True_Matr>250,100000,Dist_True_Matr)
Dist_Matr_Low <- ifelse(Dist_True_Matr>25,100000,Dist_True_Matr)

Dist_Matr_Null <- ifelse(Dist_True_Matr>0,1,0)

Dist_Matr_List_Disp <- list(Dist_Matr_High,Dist_Matr_Mid,Dist_Matr_Low)
Dist_Matr_List_Null <- list(Dist_Matr_Null,Dist_Matr_Null,Dist_Matr_Null)

# KEY POINT- Scenarios order ####
# D50 corresponding to the distance at which probability is 50%
# Scenarios order is Env, Spa, Both, Null 

dispersal_test <- c(1, # No dispersal differences - Null dispersal
                    100, # Dispersal of 100 - Dispersal activated
                    100, # Dispersal of 100 - Dispersal activated
                    1) # No dispersal differences - Null dispersal
Dist_Matr_List <- list(Dist_Matr_List_Null, # No distance differences - Null dispersal
                       Dist_Matr_List_Disp, # Distance with dispersal abilities - Dispersal  activated
                       Dist_Matr_List_Disp, # Distance with dispersal abilities - Dispersal  activated
                       Dist_Matr_List_Null) # No distance differences - Null dispersal
Applied_Filter <- list(Filter_Scen, # Differences in filters - Environment  activated
                       Filter_NO_Filter, # NO differences in filters - Null environment
                       Filter_Scen, # Differences in filters - Environment  activated
                       Filter_NO_Filter) # NO differences in filters - Null environment

# 2.4 Metacommunity model ####
# Metacommunity simulation under different Filter, Distance and Dispersal linked to Scenari (Env, Spa, Both, Null)

output <- H2020_Coalescent.and.lottery.exp.Kernel.J_TempMtcom_tempIT(
    Meta.pool = pool_200, # Species pool
    m.pool = 0.001, # Regional dispersal which is always constant 
    # Size of the communities (AKA: number of individuals/population contained in each community)
    Js = J.freshwater-ceiling((J.freshwater*((2/(apply(Applied_Filter[[Scenari]],2,mean))))/100)), 
    id.module = id_NOmodule, # id of modules if there are some - NOT used for us
    filter.env = Applied_Filter[[Scenari]], # Pollution scenarios (created at 2. Pollution assignation.R)
    Disp_Strat=c(rep(1,Gamma_div/3),rep(2,Gamma_div/3),rep(3,Gamma_div/3)),
    M.dist =list(Dist_Matr_List[[Scenari]][[1]],
                 Dist_Matr_List[[Scenari]][[2]],
                 Dist_Matr_List[[Scenari]][[3]]), # Distance matrix which corresponds to the STconmat (created at 1. OCnet - STconmat.R)
    D50 = dispersal_test[Scenari], # Dispersal distance scenario 
    m.max = 1, # Maximum migration
    tempo_imp = 0, # Relvance of "temporal" effect
    temp_Metacom = Meta_t0, # Metacommunity at time 0 (all species are equally favored)
    temp_it = 0, # Number of temporal iterations
    id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=pool_200, # If there are some communit. that should be fixed
    Lottery=T, 
    it=500, 
    prop.dead.by.it=0.05, # Lottery parameters, nº iterations and proportion of dead organisms  
    id.obs=1:nrow(Lakes_coord)) # Information if we would like to keep specific results only

Community_Scenarios[[Scenari]] <- output[[2]] #Out_Community
Rand_Filt_list[[Scenari]] <- Rand_Filt
}# Scenari end

# Once we have obtained the metacommunities for all the scenarios we can proceed to calculate NATs on each one of them. 
# For each metacommunity we are going to calculate NATs by adding communities based on:
# - Environmental similarity
# - Distance similarity 
# - Random 
# These values should illustrate the different responses of communities to the four different scenarios being
# different or similar to "randomly" add communities and identifying when one or two drivers are having 
# a stronger role in shaping functional assembly. 

# 3. NATS ####
# NATS calculation for each metacommunity
fake_output_Scenari<- data.frame()
fake_output_traits_str_Scenari <- data.frame()
for (Scenari in 1:length(Scenarios)) {
Out_Community<- Community_Scenarios[[Scenari]]

### FIRST LOOP - Scenarios 
Out_Community <- as.data.frame(t(Out_Community)) %>% mutate(Site_ID=1:nrow(.),.before=V1) %>% 
                 pivot_longer(cols = 2:ncol(.))
Out_Community<- Out_Community%>%filter(value>0)

TypeNAT=c("Rand","Dis_Env","Dis_Dist","Env","Dist")

fake_output_Total<- data.frame()
fake_output_traits_str_Total <- data.frame()
for (Type_NATs in 1:length(TypeNAT)) {
fake_output <- data.frame()
fake_output_traits_str <- data.frame()
for (ind_year in 1:1) {
cat("We are at", TypeNAT[Type_NATs],
    "from Scenari",Scenari,"/",length(Scenarios), 
    "and Replicate",Replicates,"/",Number_Of_Replicates,"___","\n")
  
### SECOND LOOP - Number of randomly selected lakes
for (iteration in 1:15) {
cat("We are at iteration", iteration,"__________________________________________________","\n")
  
  if(TypeNAT[[Type_NATs]]=="Rand"){
  # We randomly select lake names
  select_lakes <- sample(seq(1:ncol(output[[2]])), # The vector we want to select things from 
                         size =length(unique(Out_Community$Site_ID)), # The number of elements that we want to select
                         replace = F) # If we can repeat or not  
  }  

  if(TypeNAT[[Type_NATs]]=="Env"){
  Chosen_Beg <- sample(seq(1:ncol(output[[2]])), size =1,  replace = F)
  Chosen_OneS <- Chosen_Beg
  Chosen_One<- Rand_Filt_list[[Scenari]][Chosen_Beg]
  Sel_Neig <- c()
  Sel_Neig[1] <- Chosen_Beg
  for (Chosing_Lakes in 2:ncol(output[[2]])) {
  Comm_to_Sel <- seq(1:ncol(output[[2]]))[-Chosen_OneS]
  Filt_Diff <- Rand_Filt_list[[Scenari]][-Chosen_OneS]
  Sel_Neig[Chosing_Lakes] <- Comm_to_Sel[order(abs(Chosen_One-Rand_Filt_list[[Scenari]][-Chosen_OneS]))[1]]
  Chosen_OneS <- c(Sel_Neig)
  Chosen_One<- mean(Rand_Filt_list[[Scenari]][Chosen_OneS])
  }
  select_lakes <- Chosen_OneS
  }
  
  if(TypeNAT[[Type_NATs]]=="Dist"){
    lake_geo <- Lakes_coord[,1:3]
    lake_to_start <- sample(lake_geo$Lake,1)

    Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
    Llista_llacs[1] <- Lakes_coord$Lake[which(Lakes_coord$Lake==lake_to_start)]#posem el 1r llac al vector
    lake_geo[which(Lakes_coord$Lake==lake_to_start),1] <- "Chosen_One"
    #Chosen_One <- lake_to_start
    for (Chosing_Lakes in 2:nrow(Lakes_coord)) {
      xy <- lake_geo[,2:3]  
      Dist_Matr <- distm(xy)
      Chosen_OneS <- which(lake_geo$Lake=="Chosen_One") # Posició del Orig lake
      Chosen_col <- Dist_Matr[,Chosen_OneS]#aillem columna de l'escollit
      Dist_Diff <-Chosen_col[-Chosen_OneS]#treiem la posicio que es 0
      Dist_Neigh <- which(Dist_Diff==min(Dist_Diff)) # Trobem la distànica minima entre Original i altres
      Closer <- Dist_Neigh[1]
      d <- subset(lake_geo, Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                    c(lake_geo[-Chosen_OneS,1])$Lake[Closer]))#filtrem els dos llacs més propers
      new_lake <- data.frame(Lake="Chosen_One",
                             summarise(d,Lon=mean(Lon),Lat=mean(Lat)))
      Llista_llacs[Chosing_Lakes] <- c(lake_geo[-Chosen_OneS,1])$Lake[Closer]
      
      lake_geo <- lake_geo %>% filter(!Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                               c(lake_geo[-Chosen_OneS,1])$Lake[Closer])) %>% 
        bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    }
select_lakes <- c()
for (Pos_Lak in 1:length(Lakes_coord$Lake)) {
select_lakes[Pos_Lak] <- which(Lakes_coord$Lake==Llista_llacs[Pos_Lak])
}} # If Spa ending 

  if(TypeNAT[[Type_NATs]]=="Dis_Env"){
    Chosen_Beg <- sample(seq(1:ncol(output[[2]])), size =1,  replace = F)
    Chosen_OneS <- Chosen_Beg
    Chosen_One<- Rand_Filt_list[[Scenari]][Chosen_Beg]
    Sel_Neig <- c()
    Sel_Neig[1] <- Chosen_Beg
    for (Chosing_Lakes in 2:ncol(output[[2]])) {
      Comm_to_Sel <- seq(1:ncol(output[[2]]))[-Chosen_OneS]
      Filt_Diff <- Rand_Filt_list[[Scenari]][-Chosen_OneS]
      Sel_Neig[Chosing_Lakes] <- Comm_to_Sel[order(abs(Chosen_One-Rand_Filt_list[[Scenari]][-Chosen_OneS]))[length(Filt_Diff)]]
      Chosen_OneS <- c(Sel_Neig)
      Chosen_One<- mean(Rand_Filt_list[[Scenari]][Chosen_OneS])
    }
    select_lakes <- Chosen_OneS
  }

  if(TypeNAT[[Type_NATs]]=="Dis_Dist"){
    lake_geo <- Lakes_coord[,1:3]
    lake_to_start <- sample(lake_geo$Lake,1)
    
    Llista_llacs <- c() #obrim un vestor on guradar els llacs per ordre
    Llista_llacs[1] <- Lakes_coord$Lake[which(Lakes_coord$Lake==lake_to_start)]#posem el 1r llac al vector
    lake_geo[which(Lakes_coord$Lake==lake_to_start),1] <- "Chosen_One"
    #Chosen_One <- lake_to_start
    for (Chosing_Lakes in 2:nrow(Lakes_coord)) {
      xy <- lake_geo[,2:3]  
      Dist_Matr <- distm(xy)
      Chosen_OneS <- which(lake_geo$Lake=="Chosen_One") # Posició del Orig lake
      Chosen_col <- Dist_Matr[,Chosen_OneS]#aillem columna de l'escollit
      Dist_Diff <-Chosen_col[-Chosen_OneS]#treiem la posicio que es 0
      Dist_Neigh <- which(Dist_Diff==max(Dist_Diff)) # Trobem la distànica máxima entre Original i altres
      Closer <- Dist_Neigh[1]
      d <- subset(lake_geo, Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                    c(lake_geo[-Chosen_OneS,1])$Lake[Closer]))#filtrem els dos llacs més propers
      new_lake <- data.frame(Lake="Chosen_One",
                             summarise(d,Lon=mean(Lon),Lat=mean(Lat)))
      Llista_llacs[Chosing_Lakes] <- c(lake_geo[-Chosen_OneS,1])$Lake[Closer]
      
      lake_geo <- lake_geo %>% filter(!Lake%in%c(lake_geo$Lake[Chosen_OneS],
                                               c(lake_geo[-Chosen_OneS,1])$Lake[Closer])) %>% 
        bind_rows(new_lake) #treiem de la llista els dos llacs dels quals hem fer mean, i afegim "new lake" (mean dels dos llacs)
    }
    select_lakes <- c()
    for (Pos_Lak in 1:length(Lakes_coord$Lake)) {
      select_lakes[Pos_Lak] <- which(Lakes_coord$Lake==Llista_llacs[Pos_Lak])
    }
    }# If Dis_Spa ending
  
for (selected_lakes in c(1,5,9,11,13,15,20,35,45,56)) {#
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

source("NAT/function_to_NATs.R")

# We have our database! We can calculate the network! 
NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind_abund)

edge_dens <- igraph::edge_density(NATs_Output$Graph)
mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
gr_Stre <- igraph::graph.strength(NATs_Output$Graph)
transit_W <- mean(igraph::transitivity(NATs_Output$Graph, type = "barrat"),na.rm=TRUE)
transit <- igraph::transitivity(NATs_Output$Graph, type = "global")

connected <- igraph::induced_subgraph(NATs_Output$Graph, 
                                      igraph::V(NATs_Output$Graph)[igraph::degree(NATs_Output$Graph) > 0])
N_modules <- length(unique(igraph::membership(igraph::cluster_louvain(connected))))

temp_output <- data.frame("Replicates"=Replicates,
                          "Scenario"=Scenarios[Scenari],
                          "Type_NATS"=TypeNAT[[Type_NATs]],
                          "n_sites"=selected_lakes,
                          "iter"=iteration,
                          "richness"=sp_rich,
                          edge_dens,
                          mean_grStre,
                          N_modules,
                          transit_W,
                          transit)
fake_output <- bind_rows(fake_output,temp_output)
temp_output_traits_str <- data.frame("Replicates"=Replicates,
                                     "Scenario"=Scenari,
                                     "Type_NATS"=TypeNAT[[Type_NATs]],
                                     "n_sites"=selected_lakes,
                                     "iter"=iteration,
                                      "Tra_gr_Stre"=gr_Stre, 
                                      "Tra"=colnames(Trait_Matrix)[2:ncol(Trait_Matrix)])
fake_output_traits_str <- bind_rows(fake_output_traits_str,temp_output_traits_str)
    }# End of selected_lakes
  }# End of iteration
}# End of ind_year

fake_output_Total <- bind_rows(fake_output_Total,fake_output)
fake_output_traits_str_Total <- bind_rows(fake_output_traits_str_Total,fake_output_traits_str)
} # End NATs

fake_output_Scenari<- bind_rows(fake_output_Scenari,fake_output_Total)
fake_output_traits_str_Scenari <- bind_rows(fake_output_traits_str_Scenari,fake_output_traits_str_Total)
}# End of Scenari

out <- list("Out_NAT"=fake_output_Scenari,"Out_NAT_Trait"=fake_output_traits_str_Scenari)
out
}# End of Replicates

save(out,file = "NAT/OUT_fake_NATs.RData")
