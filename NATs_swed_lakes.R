library(cooccur) # calculate coocurrance matrices 
library(tidyverse) # managing data.frames and organize and edit them
library(qgraph) # qgraph is a package for generating graphs and containing different types of graph settings
library(network) # a package to manage and create "network" like objects to plot
library(ggnetwork) # plotting networks in a ggplot environment
library(viridis) # colours for colour-sensitive persons 
library(igraph) # classic and mostly used package for network calculation 

source("function_to_NATs.R")

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
for (iteration in 1:3) {
cat("We are at iteration", iteration,"__________________________________________________","\n")
# We randomly select lake names
select_lakes <- sample(year_lakes, # The vector we want to select things from 
                       size =rand_selection, # The number of elements that we want to select
                       replace = F) # If we can repeat or not

macros_lakes_list_temp_temp <- macros_lakes_list_temp %>% filter(site%in%select_lakes)


# Check presence
db_indv.l <- macros_lakes_list_temp_temp %>% filter(genus%in%traits$Genus..if.description.at.this.level.) # Filter Lakes genus from trait database 
 
traits_ind <- db_indv.l %>% 
              ungroup() %>% # just in case 
              #group_by(genus) %>% # we group by genus 
              #summarise(tot_ab=sum(inv.l)) %>% # we sum the abundance of genus in the whole network
              left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% # join traits with the genus of the traits database
              select(c(6:ncol(.))) %>% # select the columns with only traits
              na.omit() #%>% # NA elimination
              #tibble::column_to_rownames("genus") # we attach genus as rownames

# We have our database! We can calculate the network! 
NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind )

# We have two outputs of the function: 
## The table that is the "newtork" in matricial format NATs_Output$Adj_Table
## The graph object ready to be used with the "igraph" package NATs_Output$Graph 

# We now calculate different metrics that can be used to characterize the compelxity of the network 
## Check the definition on the help for the igraph package. Overall we are looking at how much relevance 
## the nodes or the whole network are/is having. 
edge_dens <- igraph::edge_density(NATs_Output$Graph)
mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
sd_grStre <- sd(igraph::graph.strength(NATs_Output$Graph))
medi_grStre <- median(igraph::graph.strength(NATs_Output$Graph))
min_grStre <- min(igraph::graph.strength(NATs_Output$Graph))
max_grStre <- max(igraph::graph.strength(NATs_Output$Graph))

## Here we have calculated everything that we needed already, so now it is time to "close" the whole process
# we will first create a data.frame with all the information 

temp_output <- data.frame("Year"=years[ind_year],
                          "n_sites"=rand_selection,
                          "iter"=iteration,
                          edge_dens,
                          mean_grStre,
                          sd_grStre,
                          medi_grStre,
                          min_grStre, 
                          max_grStre)

output <- bind_rows(output,temp_output)
write.csv2(output, file="Result_NATs.csv")

# We obtain the IDs of the loop and hte names of the lakes used to built the NATs
out_Names <- c(years[ind_year],rand_selection,iteration,unique(macros_lakes_list_temp_temp$site))
# We store this information in the matrix that we created
LakesMergedLakes <- rbind(LakesMergedLakes,c(out_Names,rep(NA,(ncol(LakesMergedLakes)-length(out_Names)))))
}# End of iteration
}# End of rand_selection
}# End of ind_year

# We polish and arrange the matrix in order to have a "nice" matrix.  
LakesMergedLakes <- LakesMergedLakes
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






