# Beginning ####
library(tidyverse)
library(geosphere)
library(doParallel)
library(parallel)

Final_Output <- list()

All_ecosyst <- c("lake","river")
for (Type_of_ecosyst in 1:length(All_ecosyst)) {
ecosyst <- All_ecosyst[Type_of_ecosyst]
  
  if(ecosyst=="lake"){
  # Traits database (canviar a filtrat per generes)
  traits <-read.csv2("ClustNAT/data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
  traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=3,1,0)#all afiliations higher than 1 have a 1
  
  # We transform and caclulate the mean for each Genus and ensure that the values are equal 1
  traits <- traits %>% group_by(Genus..if.description.at.this.level.) %>% 
                       mutate_if(is.numeric, ~mean(.)) %>% 
                       mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
                       summarise_if(is.numeric, mean, na.rm = TRUE)
  
  # Load dataset from retromed
  macros_lakes_list <- readxl::read_excel("ClustNAT/data/macros_lakes_list.xlsx") %>%
                        filter(inv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.
  }

  if(ecosyst=="river"){
  # Traits database (canviar a filtrat per generes)
  traits <-read.csv("ClustNAT/data/tachet.traits.def_mod.csv", header=TRUE, sep=";", na.strings="")
  
  # Transform traits to 1 or 0 (losing affiliations)
  traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=3,1,0)#all afiliations higher than 1 have a 1
  
  summary(apply(traits[,11:ncol(traits)],1, sum))
  traits_gen <- traits %>% group_by(Genus..if.description.at.this.level.) %>% 
    mutate_if(is.numeric, ~mean(.)) %>% 
    mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)
  traits_fam <- traits %>% group_by(Family) %>% 
    mutate_if(is.numeric, ~mean(.)) %>% 
    mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
    summarise_if(is.numeric, mean, na.rm = TRUE)
  
  # Load dataset from retromed
  macros_lakes_list <- readxl::read_excel("ClustNAT/data/swed_river_list_s_ph.xlsx") %>%
    filter(indv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.
  }

# Sequence of additions
Addition_sequence<-c(1:15,seq(17,(length(unique(macros_lakes_list$site))-1),5),length(unique(macros_lakes_list$site)))
#Addition_sequence <- c(1,3,5,7,9,11,13,15,20,25,35,45,56)
#Number_of_iterations<-length(unique(macros_lakes_list$site))
Number_of_iterations <- 56

# All years 
years <- unique(macros_lakes_list$year)
years <- years

### FIRST LOOP - Scenarios 
TypeNAT=c("Random","Environment","Distance")#,"Dis_Distance","Dis_Environment")
  
cl <- detectCores() #Number of cores in computer
registerDoParallel(cl)

out <- foreach(Type_NATs=1:length(TypeNAT))%:%foreach(ind_year=1:length(years))%dopar% {
#for (Type_NATs in 1:length(TypeNAT)) {
  NAT_output <- data.frame()
  NAT_output_traits_str <- data.frame()
  
  # We create a matrix to store the names of the selected lakes for later carry the dbFD. 
  # We need to create a matrix in order to set the number of columns a priori (that will be the total lenght of possible habitat names)
  LakesMergedLakes <- matrix(ncol = length(unique(macros_lakes_list$site))+3, data = NA)
  
  library(tidyverse);library(geosphere)# Somehow we need to recharge the packages again.
  # We charge the functions to calculate all the different patterns
  source("ClustNAT/Functions/function_to_NATs.R"); source("ClustNAT/Functions/function_to_Cooccur.R")  
  source("ClustNAT/Functions/function_to_ENV_SIM.R");  source("ClustNAT/Functions/function_to_DIST_SIM.R")
  source("ClustNAT/Functions/function_to_ENV_DISIM.R");  source("ClustNAT/Functions/function_to_DIST_DISIM.R")
  
#for (ind_year in 1:length(years)) {
    #cat("We are at", TypeNAT[Type_NATs], "and year", ind_year, "of", length(years),"___","\n")
    
    macros_lakes_list_temp <- macros_lakes_list %>% filter(year==years[ind_year])
    # All lakes of that year
    year_lakes <- unique(macros_lakes_list_temp$site)
    
    ### SECOND LOOP - Number of randomly selected lakes
    for (iteration in 1:Number_of_iterations) {
      #cat("We are at iteration", iteration,"__________________________________________________","\n")

      # We randomly select lake names
      if(TypeNAT[Type_NATs]=="Random"){
        select_lakes <- sample(year_lakes, # The vector we want to select things from 
                               size =length(year_lakes), # The number of elements that we want to select
                               replace = F) # If we can repeat or not  
      }
      
      # We select lake names by similar environments
      if(TypeNAT[Type_NATs]=="Environment"){
        year_to_start <- sample(year_lakes,1)
        if(ecosyst=="lake"){select_lakes <- fun_to_ENV_SIM(ref_year = years[ind_year],orig_lake = year_to_start)}
        if(ecosyst=="river"){select_lakes <- fun_to_ENV_SIM_riv(ref_year = years[ind_year],orig_lake = year_to_start)}
      }# Envi
      # We select lake names by DISsimilar environments
      if(TypeNAT[Type_NATs]=="Dis_Environment"){
        year_to_start <- sample(year_lakes,1)
        if(ecosyst=="lake"){select_lakes <- fun_to_ENV_DISIM(ref_year = years[ind_year],orig_lake = year_to_start)}
        if(ecosyst=="river"){select_lakes <- fun_to_ENV_DISIM_riv(ref_year = years[ind_year],orig_lake = year_to_start)}
      }
      # We select lake names by similar distance
      if(TypeNAT[Type_NATs]=="Distance"){
        year_to_start <- sample(year_lakes,1)
        if(ecosyst=="lake"){select_lakes <- fun_to_DIST_SIM(orig_lake = year_to_start)}
        if(ecosyst=="river"){select_lakes <- fun_to_DIST_SIM_riv(orig_lake = year_to_start)}
      }
      # We select lake names by different distance
      if(TypeNAT[Type_NATs]=="Dis_Distance"){
        year_to_start <- sample(year_lakes,1)
        if(ecosyst=="lake"){select_lakes <- fun_to_DIST_DISIM(orig_lake = year_to_start)}
        if(ecosyst=="river"){select_lakes <- fun_to_DIST_DISIM_riv(orig_lake = year_to_start)}
      }
      
      for (selected_lakes in Addition_sequence) {#
        #cat("We have seleted", selected_lakes,"lakes","__________________________________________________","\n")
        ### THIRD LOOP - We will repeat the same thing  several times  
        real_select_lakes <- select_lakes[1:selected_lakes]
        macros_lakes_list_temp_temp <- macros_lakes_list_temp %>% filter(site%in%real_select_lakes)
        
        if(ecosyst=="lake"){
        # Check presence
        db_indv.l <- macros_lakes_list_temp_temp%>% filter(genus%in%traits$Genus..if.description.at.this.level.) # Filter Lakes genus from trait database 
        
        traits_ind <- db_indv.l %>% 
          ungroup() %>% # just in case 
          #group_by(genus) %>% # we group by genus 
          #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
          left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all")  # join traits with the genus of the traits database
        }

        if(ecosyst=="river"){
          # Check presence
          db_indv.l1 <- macros_lakes_list_temp_temp %>% filter(genus_def%in%traits_gen$Genus..if.description.at.this.level.) # Filter Lakes genus from trait database 
          db_indv.l2 <- macros_lakes_list_temp_temp%>% filter(family%in%traits_fam$Family) # Filter Lakes genus from trait database 
          
          traits_ind1 <- db_indv.l1 %>% rename("inv.l"="indv.l") %>% 
            ungroup() %>% # just in case 
            #group_by(genus) %>% # we group by genus 
            #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
            left_join(traits_gen, by=c("genus_def"="Genus..if.description.at.this.level."),multiple ="all") %>%  # join traits with the genus of the traits database
            dplyr::select(c(5,7,9:ncol(.)))
          traits_ind2 <- db_indv.l2 %>% rename("inv.l"="indv.l")%>% 
            ungroup() %>% # just in case 
            #group_by(genus) %>% # we group by genus 
            #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
            left_join(traits_fam, by=c("family"="Family"),multiple ="all") %>%
            dplyr::select(c(5,7,9:ncol(.)))
          
          traits_ind<- bind_rows(traits_ind1, traits_ind2)
          traits_ind <- unique(traits_ind)
        }
        
        
        # Normal NATs
        sp_rich <- nrow(traits_ind)
        traits_ind_abund <- data.frame()
        for (Row_ID in 1:nrow(traits_ind)) {
          Out_sel_row <- data.frame()
          for (Row_Abun in 1:ceiling(traits_ind$inv.l[Row_ID])) {
            Selec_Row <- traits_ind[Row_ID,] 
            Out_sel_row <- bind_rows(Out_sel_row,Selec_Row)
          }
          traits_ind_abund <- bind_rows(traits_ind_abund,Out_sel_row)
        }  
        
        traits_ind_abund <-traits_ind_abund %>%dplyr::select(c(6:ncol(.))) %>% na.omit()
        Names_Traits <- colnames(traits_ind_abund)
        # We have our database! We can calculate the network! 
        NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind_abund)
        
        edge_dens <- igraph::edge_density(NATs_Output$Graph)
        mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
        gr_Stre <- igraph::graph.strength(NATs_Output$Graph)
        g <- igraph::graph.adjacency(1/(NATs_Output$Adj_Table),weighted = TRUE)
        betw_Stre <- igraph::betweenness(g) 
        clos_Stre <- igraph::harmonic_centrality(g)
        
        transit_W <- mean(igraph::transitivity(NATs_Output$Graph, type = "barrat"),na.rm=TRUE)
        transit <- igraph::transitivity(NATs_Output$Graph, type = "global")
        
        connected <- igraph::induced_subgraph(NATs_Output$Graph, 
                                              igraph::V(NATs_Output$Graph)[igraph::degree(NATs_Output$Graph) > 0])
        N_modules <- length(unique(igraph::membership(igraph::cluster_louvain(connected))))
        
        temp_output <- data.frame("Year"=years[ind_year],
                                  "Type_NATS"=TypeNAT[[Type_NATs]],
                                  "n_sites"=selected_lakes,
                                  "iter"=iteration,
                                  "richness"=sp_rich,
                                  edge_dens,
                                  mean_grStre,
                                  N_modules,
                                  transit_W,
                                  transit)
        NAT_output <- bind_rows(NAT_output,temp_output)
        
        temp_output_traits_str <- data.frame("Year"=years[ind_year],
                                             "Type_NATS"=TypeNAT[[Type_NATs]],
                                             "n_sites"=selected_lakes,
                                             "iter"=iteration,
                                             gr_Stre, 
                                             betw_Stre,
                                             clos_Stre,
                                             "Tra"=Names_Traits)
        NAT_output_traits_str <- bind_rows(NAT_output_traits_str,temp_output_traits_str)
        
        # We obtain the IDs of the loop and hte names of the lakes used to built the NATs
        out_Names <- c(years[ind_year],selected_lakes,iteration,unique(macros_lakes_list_temp_temp$site))
        # We store this information in the matrix that we created
        LakesMergedLakes <- rbind(LakesMergedLakes,c(out_Names,rep(NA,(ncol(LakesMergedLakes)-length(out_Names)))))
        
      }# End of selected_lakes
    }# End of iteration
    out <- list("Out_NAT"=NAT_output,"Out_NAT_Trait"=NAT_output_traits_str, "LakMergLak"=LakesMergedLakes)
    out
  }# Parallel ending

NAT_output_TEMP <- data.frame()
NAT_output_traits_str_TEMP <- data.frame()
LakesMergedLakes <- matrix(ncol = length(unique(macros_lakes_list$site))+3, data = NA)
for (Type_NATs in 1:length(TypeNAT)) {
for (ind_year in 1:length(years)) {
NAT_output_TEMP <- bind_rows(NAT_output_TEMP,out[[Type_NATs]][[ind_year]]$Out_NAT)
NAT_output_traits_str_TEMP <- bind_rows(NAT_output_traits_str_TEMP,out[[Type_NATs]][[ind_year]]$Out_NAT_Trait)
LakesMergedLakes <- rbind(LakesMergedLakes,out[[Type_NATs]][[ind_year]]$LakMergLak)
  }# ind_years 
}# Type_NATs

Final_O_put <- list("Out_NAT"=NAT_output_TEMP,"Out_NAT_Trait"=NAT_output_traits_str_TEMP, "LakMergLak"=LakesMergedLakes)
Final_Output[[Type_of_ecosyst]] <- Final_O_put
}# ecosyst
names(Final_Output) <- All_ecosyst

save(Final_Output,file =  "ClustNAT/ClusterNATs.RData")
