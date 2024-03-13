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

traits <- traits %>% group_by(Genus..if.description.at.this.level.) %>% 
                     mutate_if(is.numeric, ~mean(.)) %>% 
                     mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
                     summarise_if(is.numeric, mean, na.rm = TRUE)


# Load dataset from retromed
macros_lakes_list <- readxl::read_excel("data/Lakes/macros_lakes_list.xlsx") %>%
                     filter(inv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.

# All years 
years <- unique(macros_lakes_list$year)
# We create the data.frame where we will store everything during the loops
output <- data.frame()

# We create a matrix to store the names of the selected lakes for later carry the dbFD. 
# We need to create a matrix in order to set the number of columns a priori (that will be the total lenght of possible habitat names)
LakesMergedLakes <- matrix(ncol = length(unique(macros_lakes_list$site))+3, data = NA)

### FIRST LOOP - Years
for (ind_year in 1:length(years)) {
cat("We are at year", years[ind_year],"__________________________________________________","\n")
# We create a "temporary" file filtered according to the year selected. 
macros_lakes_list_temp <- macros_lakes_list %>% filter(year==years[ind_year])
# All lakes of that year
year_lakes <- unique(macros_lakes_list_temp$site)

### SECOND LOOP - Number of randomly selected lakes
for (rand_selection in seq(1,length(year_lakes),10)) {
cat("We have seleted", rand_selection,"lakes","__________________________________________________","\n")
### THIRD LOOP - We will repeat the same thing  several times  
for (iteration in 1:2) {
cat("We are at iteration", iteration,"__________________________________________________","\n")
# We randomly select lake names
select_lakes <- sample(year_lakes, # The vector we want to select things from 
                       size =rand_selection, # The number of elements that we want to select
                       replace = F) # If we can repeat or not

macros_lakes_list_temp_temp <- macros_lakes_list_temp %>% filter(site%in%select_lakes)


# Check presence
# db_indv.l <- macros_lakes_list_temp_temp %>% filter(genus%in%traits$Genus..if.description.at.this.level.) # Filter Lakes genus from trait database 
# 
# traits_ind <- db_indv.l %>% 
#   ungroup() %>% # just in case 
#   #group_by(genus) %>% # we group by genus 
#   #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
#   left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% # join traits with the genus of the traits database
#   select(c(6:ncol(.))) %>% # select the columns with only traits
#   na.omit() #%>% # NA elimination
#   #tibble::column_to_rownames("genus") # we attach genus as rownames
# 
# #run analysis (takes a bit)
# cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
# 
# # Generate an adjacency table from the coocurrence matrix 
# ## Create a blank table with only 0 
# adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
# for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {# do the following for each spp
#   # Summarise the strenght of the link for each genus by calulating the mean probability
#   link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur),.groups = "drop")
#   # filter the results from the coocurrence table for a species. You isolate all the species with who the 
#   ## targeted species is coocuring
#   spp_to <- cooccur.species_ccr$results %>%filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
#   # You add the values into the "adj.table" that was empty and now is being filled
#   adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
# }
# # We transpose the results to make the matrix symetrical (same from & to links)
# adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
# 
# # We built a graph from the adjacency table
# g <- igraph::graph.adjacency(adj_table,mode="undirected",weighted = TRUE)
# 
# #variables
# edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
# edge_dens <- igraph::edge_density(g)
# mean_grStre <- mean(igraph::graph.strength(g))
# sd_grStre <- sd(igraph::graph.strength(g))
# medi_grStre <- median(igraph::graph.strength(g))
# min_grStre <- min(igraph::graph.strength(g))
# max_grStre <- max(igraph::graph.strength(g))
# 
### Here we have calculated everything that we needed already, so now it is time to "close" the whole process
## we will first create a data.frame with all the information 

# temp_output <- data.frame("Year"=years[ind_year],
#                          "n_sites"=rand_selection,
#                          "iter"=iteration,
#                          edge_dens_man,
#                          edge_dens,
#                          mean_grStre,
#                          sd_grStre,
#                          medi_grStre,
#                          min_grStre, 
#                          max_grStre)
#
#output <- bind_rows(output,temp_output)
#write.csv2(output, file="Result_NATs_lake.csv")

# We obtain the IDs of the loop and hte names of the lakes used to built the NATs
out_Names <- c(years[ind_year],rand_selection,iteration,unique(macros_lakes_list_temp_temp$site))
# We store this information in the matrix that we created
LakesMergedLakes <- rbind(LakesMergedLakes,c(out_Names,rep(NA,(ncol(LakesMergedLakes)-length(out_Names)))))
}# End of iteration
}# End of rand_selection
}# End of ind_year

# We polish and arrange the matrix in order to have a "nice" matrix.  
LakesMergedLakes <- LakesMergedLakes[-1,]
colnames(LakesMergedLakes) <-c("Year","n_sites","it",rep("Lake_Name",(ncol(LakesMergedLakes)-3))) 

output %>%
  pivot_longer(cols = 4:ncol(.)) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=n_sites, colour=as.factor(Year)), width = 0.2)+
  geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()

# FUNCTIONAL BIODIVERSITY

# Each row of the LakesMergedLakes correspond to one of these combinations  
nrow(LakesMergedLakes)

# We charge the dataset and we eliminate the Lake and Year columns to make it more real to what we will have ;) 
sp_sites <- read.csv2("data/DataAnnaTest/fuzzy_traits.csv",dec = ".") %>% mutate(Year=as.character(Year)) %>% mutate_if(is.numeric,~ifelse(.>0,1,0)) %>% mutate(Year=as.numeric(Year))
load("data/DataAnnaTest/dis_traits_lakes.RData")
# What do we want now? We want to select the ponds that are listed in each row of our table no? So we need to make a "loop"
# where for each row we will select the listed lakes, filter them from the sp_sites and sum their abundance values

# We create the output dataframe that will contain all the combinations of communities
Combis_Communities <- data.frame()
for (combi in 1:nrow(LakesMergedLakes)) {
  Lakes_Combin <- LakesMergedLakes[combi,which(is.na(LakesMergedLakes[combi,])==F)] # We extract The row that corresponds to a determined combination 
  # Note that we select only the values that are NOT an NA
  
  sp_sites_Combin <- sp_sites %>% filter(Year==as.numeric(LakesMergedLakes[combi,1])) %>% # Filter the sp_sites by the year 
    filter(Lake%in%Lakes_Combin[4:length(Lakes_Combin)]) %>%  # Filter the lakes by the lakes present in the combination
    select(-c(Year,lake_year,Lake)) # We "deselect" the three cathegorical metrics
  
  # So, we have generated a dataframe that contains, all the samples (rows) from a year and a combination of ponds. Now we just need to 
  # sum by columns and transform the result into 1/0. Finally, bind it to the output dataframe and DONE! :D
  
  LakesMergedLakes_id<-LakesMergedLakes[combi,1:3]
  LakesMergedLakes_comb<-ifelse(apply(sp_sites_Combin,2,sum)>0,1,0) #%>% mutate_if(is.character, as.numeric)
  LakesMerged<-c(LakesMergedLakes_id,LakesMergedLakes_comb)
  LakesMerged[4:length(LakesMerged)]<-as.numeric(LakesMerged[4:length(LakesMerged)])
  Combis_Communities <- bind_rows(Combis_Communities,LakesMerged)
}

# Each one of the rows corresponds a each one of the used combinations and can be identified with the Year, n_sites and iteration.
# Therefore we can link it with the values of the NATs

# In case of need we could also "merge" the three identifiers: 

Combis_Communities<-Combis_Communities %>% mutate(Combi_ID=paste(Year,n_sites,it, sep="_"),.before = Year) %>% select(-c(Year,n_sites,it))

# PISTA!: Sempre que combinis "noms" com has fet amb el lake_year assegurat d'unir-los amb un "_" o un ".", d'aquesta 
# manera, en cas que volguessis separar-los seria fàcil perquè podries identificar per on "partir-los"
# La funció strsplit() va molt bé per fer això però necessita un identificador per on separar ;)

Combis_Communities <- data.frame(Combis_Communities, row.names = 1)%>% mutate_if(is.character, as.numeric)#per posar el nom de les localitats 

resFD <- dbFD( 
  dis_traits, 
  Combis_Communities, #fem el LOG si treballem amb abun, si fem P/A no cal
  corr = "cailliez",
  w.abun=TRUE,
  stand.FRic = TRUE, 
  m=7) #m=número de diM

important.indices <- cbind(resFD$nbsp, resFD$FRic, resFD$FEve, resFD$FDiv,
                           resFD$FDis, resFD$RaoQ)
colnames(important.indices) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")

bind_cols(LakesMergedLakes[,1:3], important.indices)%>%
  pivot_longer(cols = 4:ncol(.)) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=n_sites, colour=as.factor(Year)), width = 0.2)+
  geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()





