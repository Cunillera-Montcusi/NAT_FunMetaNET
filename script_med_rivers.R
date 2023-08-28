library(cooccur) # calculate coocurrance matrices 
library(tidyverse) # managing data.frames and organize and edit them
library(qgraph) # qgraph is a package for generating graphs and containing different types of graph settings
library(network) # a package to manage and create "network" like objects to plot
library(ggnetwork) # plotting networks in a ggplot environment
library(viridis) # colours for colour-sensitive persons 
library(igraph) # classic and mostly used package for network calculation 

# Traits database (canviar a filtrat per generes)
traits <-read.csv("data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")

# Transform traits to 1 or 0 (losing affiliations)
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)#all afiliations higher than 1 have a 1

traits <- traits %>% group_by(Family) %>% 
                     mutate_if(is.numeric, ~mean(.)) %>% 
                     mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
                     summarise_if(is.numeric, mean, na.rm = TRUE)
               
# Load dataset from retromed
macros_river_list <- readxl::read_excel("data/Lakes/macros_lakes_list.xlsx") %>%
  filter(indv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.

# All years 
years <- unique(macros_river_list$year)
# We create the data.frame where we will store everything during the loops
output <- data.frame()

### FIRST LOOP - Years
for (ind_year in 1:length(years)) 
# We create a "temporary" file filtered according to the year selected. 
macros_river_list_temp <- macros_river_list %>% filter(year==years[ind_year])
# All lakes of that year
year_rivers <- unique(macros_river_list_temp$site)
  
### SECOND LOOP - Number of randomly selected lakes
for (rand_selection in seq(1,length(year_rivers),2)) {
    
### THIRD LOOP - We will repeat the same thing  several times  
for (iteration in 1:10) {
      
# We randomly select lake names
select_lakes <- sample(year_lakes, # The vector we want to select things from 
                             size =rand_selection, # The number of elements that we want to select
                             replace = F) # If we can repeat or not
      
macros_lakes_list_temp_temp <- macros_lakes_list_temp %>% filter(site%in%select_lakes)

# Check presence
db_indv.l <- macros_river_list %>% filter(taxa%in%traits$Family) # Filter river families from trait database 

head(db_indv.l) #Revisar: falten 82 generes per assignar

db_indv.l <- db_indv.l %>% filter(year=="1997") # We filter to obtain the ponds of the year 1997

traits_ind <- db_indv.l %>% 
  ungroup() %>% # just in case 
  group_by(taxa) %>% # we group by genus 
  summarise(tot_ab=sum(indv.l)) %>% # we sum the abundance of genus in the whole network
  left_join(traits, by=c("taxa"="Family"),multiple ="all") %>% # join traits with the genus of the traits database
  group_by(taxa) %>% # regroup by genus
  filter(n()<2) %>% # filter those species that are listed 2 times (by mistake), check which are these in case that there are some. It means mistakes.
  select(c(1,12:ncol(.))) %>% # select the columns with only traits
  na.omit() %>% # NA elimination
  tibble::column_to_rownames("taxa") # we attach genus as rownames

head(traits_ind)


#run analysis (takes a bit)
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
# this script comes from the supplementary material of the "original" NAT's paper (DOI: 10.1002/eap.2010)
# Note that we transpose the matrix (t())

head(cooccur.species_ccr)

# Generate an adjacency table from the coocurrence matrix 
## Create a blank table with only 0 
adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {# do the following for each spp
  # Summarise the strenght of the link for each genus by calulating the mean probability
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  # filter the results from the coocurrence table for a species. You isolate all the species with who the 
  ## targeted species is coocuring
  spp_to <- cooccur.species_ccr$results %>%filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  # You add the values into the "adj.table" that was empty and now is being filled
  adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
}
# We transpose the results to make the matrix symetrical (same from & to links)
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]

# We built a graph from the adjacency table
g <- igraph::graph.adjacency(adj_table,mode="undirected",weighted = TRUE)

#Variables
edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
edge_dens <- igraph::edge_density(g)
mean_grStre <- mean(igraph::graph.strength(g))
sd_grStre <- sd(igraph::graph.strength(g))
medi_grStre <- median(igraph::graph.strength(g))
min_grStre <- min(igraph::graph.strength(g))
max_grStre <- max(igraph::graph.strength(g))

