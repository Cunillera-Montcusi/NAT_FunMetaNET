library(cooccur) # calculate coocurrance matrices 
library(tidyverse) # managing data.frames and organize and edit them
library(qgraph) # qgraph is a package for generating graphs and containing different types of graph settings
library(network) # a package to manage and create "network" like objects to plot
library(ggnetwork) # plotting networks in a ggplot environment
library(viridis) # colours for colour-sensitive persons 
library(igraph) # classic and mostly used package for network calculation 

source("function_to_NATs.R")
source("function_to_ENV_SIM.R")
source("function_to_ENV_DISIM.R")
source("function_to_DIST_SIM.R")
source("function_to_DIST_DISIM.R")

# Traits database (canviar a filtrat per generes)
traits <-read.csv("data/tachet.traits.def_mod.csv", header=TRUE, sep=";", na.strings="")
traits <- traits[,c(1:10,57:73)] #feeding gr
traits <- traits[,c(1:17,49:56)]#size i loc
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
macros_lakes_list <- readxl::read_excel("data/Rivers/swed_river_list_s_ph.xlsx") %>%
  filter(indv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.


# All years 
years <- unique(macros_lakes_list$year)
years <- years[c(11)]
# We create the data.frame where we will store everything during the loops
output <- data.frame()
output_traits_str <-  data.frame()
# We create a matrix to store the names of the selected lakes for later carry the dbFD. 
# We need to create a matrix in order to set the number of columns a priori (that will be the total lenght of possible habitat names)
LakesMergedLakes <- matrix(ncol = length(unique(macros_lakes_list$site))+4, data = NA)

#NATs_Netw <- list()
Type_of_NATs <- "Dis_Distance"
### FIRST LOOP - Years
for (ind_year in 1:length(years)) {
  cat("We are at year", years[ind_year],"__________________________________________________","\n")
  # We create a "temporary" file filtered according to the year selected. 
  macros_lakes_list_temp <- macros_lakes_list %>% filter(year==years[ind_year])
  # All lakes of that year
  year_lakes <- unique(macros_lakes_list_temp$site)
  
  #Iter_NATs_Netw <- list()
  ### THIRD LOOP - We will repeat the same thing  several times  
  for (iteration in 1:5) {
    cat("We are at iteration", iteration,"__________________________________________________","\n")
    
    # We randomly select lake names
    if(Type_of_NATs=="Random"){
      select_lakes <- sample(year_lakes, # The vector we want to select things from 
                             size =length(year_lakes), # The number of elements that we want to select
                             replace = F) # If we can repeat or not  
    }
    
    # We select lake names by similar environments
    if(Type_of_NATs=="Environment"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_ENV_SIM_riv(ref_year = years[ind_year],orig_lake = year_to_start)
    }
    # We select lake names by DISsimilar environments
    if(Type_of_NATs=="Dis_Environment"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_ENV_DISIM_riv(ref_year = years[ind_year],orig_lake = year_to_start)
    }
     #We select lake names by similar distance
    if(Type_of_NATs=="Distance"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_DIST_SIM_riv(orig_lake = year_to_start)
    }
    # We select lake names by different distance
    if(Type_of_NATs=="Dis_Distance"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_DIST_DISIM_riv(orig_lake = year_to_start)
    }
    #Sel_lak_NATs_Netw <- list()
    ### SECOND LOOP - Number of randomly selected lakes #c(1,3,5,7,9,11,13,15,20,25,39))
    for (selected_lakes in c(1,3,5,7,11,15,20,39)){#seq(1,length(year_lakes),10)) {
      cat("We have seleted", selected_lakes,"lakes","__________________________________________________","\n")  
      
      real_select_lakes <- select_lakes[1:selected_lakes]
      
      macros_lakes_list_temp_temp <- macros_lakes_list_temp %>% filter(site%in%real_select_lakes)
      
      
      # Check presence
      db_indv.l1 <- macros_lakes_list_temp_temp %>% filter(genus_def%in%traits_gen$Genus..if.description.at.this.level.) # Filter Lakes genus from trait database 
      db_indv.l2 <- macros_lakes_list_temp_temp%>% filter(family%in%traits_fam$Family) # Filter Lakes genus from trait database 
     
      traits_ind1 <- db_indv.l1 %>% 
        ungroup() %>% # just in case 
        #group_by(genus) %>% # we group by genus 
        #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
        left_join(traits_gen, by=c("genus_def"="Genus..if.description.at.this.level."),multiple ="all") %>%  # join traits with the genus of the traits database
        select(c(5,7,9:ncol(.)))
      traits_ind2 <- db_indv.l2 %>% 
        ungroup() %>% # just in case 
        #group_by(genus) %>% # we group by genus 
        #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
        left_join(traits_fam, by=c("family"="Family"),multiple ="all") %>%
        select(c(5,7,9:ncol(.)))
      
      traits_ind<- bind_rows(traits_ind1, traits_ind2)
      traits_ind <- unique(traits_ind)
      
      sp_rich <- nrow(traits_ind)
      traits_ind_abund <- data.frame()
      for (Row_ID in 1:nrow(traits_ind)) {
        Out_sel_row <- data.frame()
        for (Row_Abun in 1:ceiling(traits_ind$indv.l[Row_ID])) {
          Selec_Row <- traits_ind[Row_ID,] 
          Out_sel_row <- bind_rows(Out_sel_row,Selec_Row)
    }
        traits_ind_abund <- bind_rows(traits_ind_abund,Out_sel_row)
      }  
      
      traits_ind_abund <-traits_ind_abund %>%select(c(3:ncol(.))) %>% # select the columns with only traits
        na.omit() #%>% # NA elimination
      #tibble::column_to_rownames("genus") # we attach genus as rownames
      
      # We have our database! We can calculate the network! 
      NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind_abund)
      
      # We have two outputs of the function: 
      ## The table that is the "newtork" in matricial format NATs_Output$Adj_Table
      ## The graph object ready to be used with the "igraph" package NATs_Output$Graph 
      
      # We now calculate different metrics that can be used to characterize the compelxity of the network 
      ## Check the definition on the help for the igraph package. Overall we are looking at how much relevance 
      ## the nodes or the whole network are/is having. 
      edge_dens <- igraph::edge_density(NATs_Output$Graph)
      mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
      transit_W <- mean(igraph::transitivity(NATs_Output$Graph, type = "barrat"),na.rm=TRUE)
      transit <- igraph::transitivity(NATs_Output$Graph, type = "global")
      inter_even <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "interaction evenness"))
      module <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "modularity"))
      clust <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "cluster coefficient")[1])
      nest_w <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "weighted nestedness"))
      NODF_w <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "weighted NODF"))
      Connec_w <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "weighted connectance"))
      Mod.A <- length(bipartite::listModuleInformation(bipartite::computeModules(NATs_Output$Adj_Table)))
      gr_Stre <- igraph::graph.strength(NATs_Output$Graph)
      ## Here we have calculated everything that we needed already, so now it is time to "close" the whole process
      # we will first create a data.frame with all the information 
      
      temp_output <- data.frame("Year"=years[ind_year],
                                "n_sites"=selected_lakes,
                                "iter"=iteration,
                                edge_dens,
                                mean_grStre,
                                transit_W,
                                transit,inter_even,
                                module,clust, nest_w, NODF_w, Connec_w, Mod.A)
      
      output <- bind_rows(output,temp_output)
      temp_output_traits_str <- data.frame(#"Type_NATS"=Type_of_NATs,
        "Year"=years[ind_year],
        "n_sites"=selected_lakes,
        "iter"=iteration,
        gr_Stre, Tra=colnames(traits_ind)[3:ncol(traits_ind)])
      output_traits_str <- bind_rows(output_traits_str,temp_output_traits_str)
      #write.csv2(output, file="Result_NATs.csv")
      
      # We obtain the IDs of the loop and the names of the lakes used to built the NATs
      out_Names <- c(years[ind_year],selected_lakes,iteration,sp_rich,unique(macros_lakes_list_temp_temp$site))
      # We store this information in the matrix that we created
      LakesMergedLakes <- rbind(LakesMergedLakes,c(out_Names,rep(NA,(ncol(LakesMergedLakes)-length(out_Names)))))
      #Selcted_Lakes_order <- c(1,3,6,11,20,50)
      #Sel_lak_NATs_Netw[[which(Selcted_Lakes_order==selected_lakes)]] <- NATs_Output
    }# End of selected_lakes
    #Iter_NATs_Netw[[iteration]] <-Sel_lak_NATs_Netw
  }# End of iteration
  #NATs_Netw[[ind_year]] <-Iter_NATs_Netw
}# End of ind_year

# We polish and arrange the matrix in order to have a "nice" matrix.  
LakesMergedLakes <- LakesMergedLakes[-1,]
colnames(LakesMergedLakes) <-c("Year","n_sites","it","richness",rep("Lake_Name",(ncol(LakesMergedLakes)-4))) 
#NATs_Netw[[1]][[1]][[1]]$Adj_Table
save(output_traits_str,file = "NATs_riv_DisDist_CanvisFeedingTR17.RData")
save(output,file = "NATs_riv_Rand_CanvisLocomTR.RData")
save(LakesMergedLakes,file = "riv_for_FD_DisDist_def.RData")
out_env <- output

output %>%
  pivot_longer(cols = 4:ncol(.)) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>% filter(name=="edge_dens") %>%
  ggplot()+ 
  geom_jitter(aes(y=value, x=n_sites, colour=as.factor(Year)), width = 0.2)+
  geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()

output %>%
  group_by(Year, iter) %>% 
  mutate(Mean_val=mean(edge_dens)) %>% 
  ggplot(aes(y=edge_dens, x=n_sites, colour=as.factor(Year)))+ 
  geom_jitter(width = 0.2)+
  geom_smooth(aes(linetype=as.factor(Year)),method="loess",se=F)+ 
  theme_classic()




load("Rand_obs_NATs.RData")
Rand_Out <- output
load("Dis_obs_NATs.RData")
Dist_Out <- output
load("Env_obs_NATs.RData")
Env_Out <- output

full_output <- bind_rows(
  Rand_Out %>% mutate(TypeNAT="Rand"),
  Dist_Out %>% mutate(TypeNAT="Dist"),
  Env_Out%>% mutate(TypeNAT="Env"))

unique( full_output$Year)

full_output %>%
  #filter(Year%in%c(1995,2004)) %>% 
  group_by(TypeNAT,Year,n_sites) %>% 
  mutate(Mean_val=mean(module)) %>% 
  ggplot(aes(y=Mean_val, x=n_sites, shape=TypeNAT,colour=as.factor(iter)))+ 
  #geom_jitter(width = 0.2)+
  geom_line(aes(linetype=TypeNAT))+
  #geom_smooth(aes(linetype=TypeNAT),method="loess",se=F)+ 
  theme_classic()+facet_wrap(.~Year)


source("function_to_NATs.R")
source("function_to_ENV_SIM.R")
#source("function_to_DIST_SIM.R")
#source("function_to_ENV_DISIM_riv.R")

# Traits database (canviar a filtrat per generes)
traits <-read.csv("data/tachet.traits.def_mod.csv", header=TRUE, sep=";", na.strings="")

# Transform traits to 1 or 0 (losing affiliations)
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=3,1,0)#all afiliations higher than 1 have a 1

traits_gen <- traits %>% group_by(Genus..if.description.at.this.level.) %>% 
  mutate_if(is.numeric, ~mean(.)) %>% 
  mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)
traits_fam <- traits %>% group_by(Family) %>% 
  mutate_if(is.numeric, ~mean(.)) %>% 
  mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

# Load dataset from retromed
macros_lakes_list <- readxl::read_excel("data/Rivers/swed_river_list_s_ph.xlsx") %>%
  filter(indv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.


# All years 
years <- unique(macros_lakes_list$year)
#years <- years[c(2,5,10)]
# We create the data.frame where we will store everything during the loops
output <- data.frame()

# We create a matrix to store the names of the selected lakes for later carry the dbFD. 
# We need to create a matrix in order to set the number of columns a priori (that will be the total lenght of possible habitat names)
LakesMergedLakes <- matrix(ncol = length(unique(macros_lakes_list$site))+4, data = NA)

#NATs_Netw <- list()
Type_of_NATs <- "Distance"
### FIRST LOOP - Years
for (ind_year in 1:length(years)) {
  cat("We are at year", years[ind_year],"__________________________________________________","\n")
  # We create a "temporary" file filtered according to the year selected. 
  macros_lakes_list_temp <- macros_lakes_list %>% filter(year==years[ind_year])
  # All lakes of that year
  year_lakes <- unique(macros_lakes_list_temp$site)
  
  #Iter_NATs_Netw <- list()
  ### THIRD LOOP - We will repeat the same thing  several times  
  for (iteration in 1:20) {
    cat("We are at iteration", iteration,"__________________________________________________","\n")
    
    # We randomly select lake names
    if(Type_of_NATs=="Random"){
      select_lakes <- sample(year_lakes, # The vector we want to select things from 
                             size =length(year_lakes), # The number of elements that we want to select
                             replace = F) # If we can repeat or not  
    }
    
    # We select lake names by similar environments
    if(Type_of_NATs=="Environment"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_ENV_SIM_riv(ref_year = years[ind_year],orig_lake = year_to_start)
    }
    # We select lake names by DISsimilar environments
    if(Type_of_NATs=="Dis_Environment"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_ENV_DISIM_riv(ref_year = years[ind_year],orig_lake = year_to_start)
    }
    # We select lake names by similar distance
    if(Type_of_NATs=="Distance"){
    year_to_start <- sample(year_lakes,1)
    select_lakes <- fun_to_DIST_SIM(orig_lake = year_to_start)
    }
    #Sel_lak_NATs_Netw <- list()
    ### SECOND LOOP - Number of randomly selected lakes
    for (selected_lakes in c(1,3,5,7,9,11,13,15,20,25,39)){#seq(1,length(year_lakes),10)) {
      cat("We have seleted", selected_lakes,"lakes","__________________________________________________","\n")  
      
      real_select_lakes <- select_lakes[1:selected_lakes]
      
      macros_lakes_list_temp_temp <- macros_lakes_list_temp %>% filter(site%in%real_select_lakes)
      
      
      # Check presence
      db_indv.l1 <- macros_lakes_list_temp_temp %>% filter(genus_def%in%traits_gen$Genus..if.description.at.this.level.) # Filter Lakes genus from trait database 
      db_indv.l2 <- macros_lakes_list_temp_temp%>% filter(family%in%traits_fam$Family) # Filter Lakes genus from trait database 
      
      traits_ind1 <- db_indv.l1 %>% 
        ungroup() %>% # just in case 
        #group_by(genus) %>% # we group by genus 
        #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
        left_join(traits_gen, by=c("genus_def"="Genus..if.description.at.this.level."),multiple ="all") %>%  # join traits with the genus of the traits database
        select(c(5,7,9:ncol(.)))
      traits_ind2 <- db_indv.l2 %>% 
        ungroup() %>% # just in case 
        #group_by(genus) %>% # we group by genus 
        #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
        left_join(traits_fam, by=c("family"="Family"),multiple ="all") %>%
        select(c(5,7,9:ncol(.)))
      
      traits_ind<- bind_rows(traits_ind1, traits_ind2)
      traits_ind <- unique(traits_ind)
      
      sp_rich <- nrow(traits_ind)
      traits_ind_abund <- data.frame()
      for (Row_ID in 1:nrow(traits_ind)) {
        Out_sel_row <- data.frame()
        for (Row_Abun in 1:ceiling(traits_ind$indv.l[Row_ID])) {
          Selec_Row <- traits_ind[Row_ID,] 
          Out_sel_row <- bind_rows(Out_sel_row,Selec_Row)
        }
        traits_ind_abund <- bind_rows(traits_ind_abund,Out_sel_row)
      }  
      
      traits_ind_abund <-traits_ind_abund %>%select(c(3:ncol(.))) %>% # select the columns with only traits
        na.omit() #%>% # NA elimination
      #tibble::column_to_rownames("genus") # we attach genus as rownames
      
      # We have our database! We can calculate the network! 
      NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind_abund)
      
      # We have two outputs of the function: 
      ## The table that is the "newtork" in matricial format NATs_Output$Adj_Table
      ## The graph object ready to be used with the "igraph" package NATs_Output$Graph 
      
      # We now calculate different metrics that can be used to characterize the compelxity of the network 
      ## Check the definition on the help for the igraph package. Overall we are looking at how much relevance 
      ## the nodes or the whole network are/is having. 
      edge_dens <- igraph::edge_density(NATs_Output$Graph)
      mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
      transit_W <- mean(igraph::transitivity(NATs_Output$Graph, type = "barrat"),na.rm=TRUE)
      transit <- igraph::transitivity(NATs_Output$Graph, type = "global")
      inter_even <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "interaction evenness"))
      module <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "modularity"))
      clust <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "cluster coefficient")[1])
      nest_w <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "weighted nestedness"))
      NODF_w <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "weighted NODF"))
      Connec_w <- as.numeric(bipartite::networklevel(NATs_Output$Adj_Table,index = "weighted connectance"))
      Mod.A <- length(bipartite::listModuleInformation(bipartite::computeModules(NATs_Output$Adj_Table)))
      ## Here we have calculated everything that we needed already, so now it is time to "close" the whole process
      # we will first create a data.frame with all the information 
      
      temp_output <- data.frame("Year"=years[ind_year],
                                "n_sites"=selected_lakes,
                                "iter"=iteration,
                                edge_dens,
                                mean_grStre,
                                transit_W,
                                transit,inter_even,
                                module,clust, nest_w, NODF_w, Connec_w, Mod.A)
      
      output <- bind_rows(output,temp_output)
      #write.csv2(output, file="Result_NATs.csv")
      
      # We obtain the IDs of the loop and the names of the lakes used to built the NATs
      out_Names <- c(years[ind_year],selected_lakes,iteration,sp_rich,unique(macros_lakes_list_temp_temp$site))
      # We store this information in the matrix that we created
      LakesMergedLakes <- rbind(LakesMergedLakes,c(out_Names,rep(NA,(ncol(LakesMergedLakes)-length(out_Names)))))
      #Selcted_Lakes_order <- c(1,3,6,11,20,50)
      #Sel_lak_NATs_Netw[[which(Selcted_Lakes_order==selected_lakes)]] <- NATs_Output
    }# End of selected_lakes
    #Iter_NATs_Netw[[iteration]] <-Sel_lak_NATs_Netw
  }# End of iteration
  #NATs_Netw[[ind_year]] <-Iter_NATs_Netw
}# End of ind_year

# We polish and arrange the matrix in order to have a "nice" matrix.  
LakesMergedLakes <- LakesMergedLakes[-1,]
colnames(LakesMergedLakes) <-c("Year","n_sites","it","richness",rep("Lake_Name",(ncol(LakesMergedLakes)-4))) 
#NATs_Netw[[1]][[1]][[1]]$Adj_Table
save(output,file = "NATs_riv_dist_def.RData")
save(LakesMergedLakes,file = "riv_for_FD_dist_def.RData")



