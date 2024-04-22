library(cooccur) # calculate coocurrance matrices 
library(tidyverse) # managing data.frames and organize and edit them
library(qgraph) # qgraph is a package for generating graphs and containing different types of graph settings
library(network) # a package to manage and create "network" like objects to plot
library(ggnetwork) # plotting networks in a ggplot environment
library(viridis) # colours for colour-sensitive persons 
library(igraph) # classic and mostly used package for network calculation 

source("function_to_NATs.R")
source("function_to_ENV_SIM.R")
source("function_to_DIST_SIM.R")
source("function_to_ENV_DISIM.R")

# Traits database (canviar a filtrat per generes)
traits <-read.csv("data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits <- traits[,c(1:10,57:73)] #feeding gr
traits <- traits[,c(1:17,49:56)]#size i loc
traits <- traits[,c(1:10,110:112,123:128)]#pH i trphState
# Transform traits to 1 or 0 (losing affiliations)
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=3,1,0)#all afiliations higher than 1 have a 1

traits <- traits %>% group_by(Genus..if.description.at.this.level.) %>% 
  mutate_if(is.numeric, ~mean(.)) %>% 
  mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)


# Load dataset from retromed
macros_lakes_list <- readxl::read_excel("data/Lakes/macros_lakes_list.xlsx") %>%
  filter(inv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.

# All years 
sites <- unique(macros_lakes_list$site)
sites <- sites[c(5,10,20,30,40)]

# We create the data.frame where we will store everything during the loops
output <- data.frame()

# We create a matrix to store the names of the selected lakes for later carry the dbFD. 
# We need to create a matrix in order to set the number of columns a priori (that will be the total lenght of possible habitat names)
LakesMergedLakes <- matrix(ncol = length(unique(macros_lakes_list$year))+3, data = NA)
#NATs_Netw <- list()
Type_of_NATs <- "Temporal"
### FIRST LOOP - Years
for (ind_site in 1:length(sites)) {
  cat("We are at site", ind_site,"of",length(sites),"__________________________________________________","\n")
  # We create a "temporary" file filtered according to the site selected. 
  macros_lakes_list_temp <- macros_lakes_list %>% filter(site==sites[ind_site])
  # All years of that site
  year_lakes <- unique(macros_lakes_list_temp$year)
  
  #Iter_NATs_Netw <- list()
  ### THIRD LOOP - We will repeat the same thing  several times  
  #for (iteration in 1:2) {
    #cat("We are at iteration", iteration,"__________________________________________________","\n")
    
    # We randomly select lake names
    if(Type_of_NATs=="Random"){
      select_lakes <- sample(year_lakes, # The vector we want to select things from 
                             size =length(year_lakes), # The number of elements that we want to select
                             replace = F) # If we can repeat or not  
    }
    
    # We select lake names by similar environments
    if(Type_of_NATs=="Environment"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_ENV_SIM(ref_year = years[ind_year],orig_lake = year_to_start)
    }
    # We select lake names by DISsimilar environments
    if(Type_of_NATs=="Dis_Environment"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_ENV_DISIM(ref_year = years[ind_year],orig_lake = year_to_start)
    }
    # We select lake names by similar distance
    if(Type_of_NATs=="Distance"){
      year_to_start <- sample(year_lakes,1)
      select_lakes <- fun_to_DIST_SIM(orig_lake = year_to_start)
    }
   if(Type_of_NATs=="Temporal"){
     select_years <- year_lakes
     
  }
    #Sel_lak_NATs_Netw <- list()
    ### SECOND LOOP - Number of randomly selected lakes
    for (selected_years in seq(1,length(select_years),1)){
      cat("We have seleted", selected_years,"years","__________________________________________________","\n")  
      
      real_select_years <- select_years[1:selected_years]
      
      
      macros_lakes_list_temp_temp <- macros_lakes_list_temp %>% filter(year%in%real_select_years)
      
      
      # Check presence
      db_indv.l <- macros_lakes_list_temp_temp%>% filter(genus%in%traits$Genus..if.description.at.this.level.) # Filter Lakes genus from trait database 
      
      traits_ind <- db_indv.l %>% 
        ungroup() %>% # just in case 
        #group_by(genus) %>% # we group by genus 
        #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
        left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all")  # join traits with the genus of the traits database
      
      traits_ind_abund <- data.frame()
      rich <- nrow(traits_ind)
      for (Row_ID in 1:nrow(traits_ind)) {
        Out_sel_row <- data.frame()
        for (Row_Abun in 1:ceiling(traits_ind$inv.l[Row_ID])) {
          Selec_Row <- traits_ind[Row_ID,] 
          Out_sel_row <- bind_rows(Out_sel_row,Selec_Row)
        }
        traits_ind_abund <- bind_rows(traits_ind_abund,Out_sel_row)
      }  
      
      traits_ind_abund <-traits_ind_abund %>%select(c(6:ncol(.))) %>% # select the columns with only traits
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
      
      temp_output <- data.frame("Site"=sites[ind_site],
                                "n_years"=selected_years,
                                edge_dens,
                                mean_grStre,
                                transit_W,
                                transit,
                                inter_even,
                                module,
                                clust,
                                nest_w, 
                                NODF_w, 
                                Connec_w, Mod.A)
      
      output <- bind_rows(output,temp_output)
      write.csv2(output, file="Result_NATs_temporal.csv")
      
      # We obtain the IDs of the loop and hte names of the lakes used to built the NATs
      out_Names <- c(sites[ind_site],selected_years,rich,unique(macros_lakes_list_temp_temp$year))
      # We store this information in the matrix that we created
      LakesMergedLakes <- rbind(LakesMergedLakes,c(out_Names,rep(NA,(ncol(LakesMergedLakes)-length(out_Names)))))
      #Selcted_Lakes_order <- seq(1,length(select_years),1)
      #Sel_lak_NATs_Netw[[which(Selcted_Lakes_order==selected_years)]] <- NATs_Output
    }# End of selected_years
  #NATs_Netw[[ind_site]] <-Sel_lak_NATs_Netw
}# End of ind_site

# We polish and arrange the matrix in order to have a "nice" matrix.  
LakesMergedLakes <- LakesMergedLakes[-1,]
colnames(LakesMergedLakes) <-c("Site","n_year","richness",rep("Lake_Name",(ncol(LakesMergedLakes)-3))) 
#NATs_Netw[[1]][[1]][[1]]$Adj_Table
save(LakesMergedLakes,file = "lakes_for_FD_cumTemporal_def.RData")
save(output, file="Nats_lakes_cumTemp_def.RData")


output %>%
  pivot_longer(cols = 3:ncol(.)) %>% 
  #filter(Site%in%c("Allgjuttern","Alsjön")) %>%
  group_by(Site,n_years,name) %>% 
  mutate(Mean_val=mean(value)) %>% filter(name=="edge_dens")%>%
  ggplot()+ 
  geom_point(aes(y=value, x=n_years, colour=as.factor(Site)))+
  geom_line(aes(y=Mean_val, x=n_years, colour=as.factor(Site)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()

output %>%
  #group_by(Year,n_sites) %>% 
  #mutate(Mean_val=mean(edge_dens)) %>% 
  ggplot(aes(y=Connec_w, x=n_years))+ 
  geom_jitter(width = 0.2)+
  geom_smooth(method="loess",se=F)+ 
  theme_classic()


save(output,file = "NATs_sumTemporal_lakes.RData")


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


