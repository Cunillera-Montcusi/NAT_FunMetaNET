library(cooccur) # calculate coocurrance matrices 
library(tidyverse) # managing data.frames and organize and edit them
library(qgraph) # qgraph is a package for generating graphs and containing different types of graph settings
library(network) # a package to manage and create "network" like objects to plot
library(ggnetwork) # plotting networks in a ggplot environment
library(viridis) # colours for colour-sensitive persons 
library(igraph) # classic and mostly used package for network calculation 


# GENERALISTS #### 
# Traits database (canviar a filtrat per generes)
traits_ori <-read.csv("data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")

# Transform traits to 1 or 0 (losing affiliations)
traits_ori[,11:ncol(traits_ori)] <- ifelse(traits_ori[,11:ncol(traits_ori)]>=1,1,0)#all afiliations higher than 1 have a 1

traits_ori <- traits_ori %>% group_by(Genus..if.description.at.this.level.) %>% 
  mutate_if(is.numeric, ~mean(.)) %>% 
  mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

sum_traits <- as.numeric(apply(traits_ori[,2:ncol(traits_ori)],2,sum))
traits_ori[which(is.na(traits_ori[,(which(is.na(sum_traits)==T)+1)])==T),(which(is.na(sum_traits)==T)+1)] <- 0


par(mfrow=c(4,8))
Tra_to_plot_ori <- as.numeric(apply(traits_ori[,2:ncol(traits_ori)],2,sum))
plot(Tra_to_plot_ori,ylim = c(0,460),
xlab = "Traits",ylab = "Species",main = "How many species have the same trait")
abline(h = mean(Tra_to_plot_ori),col="red",lwd=3)
abline(h = quantile(Tra_to_plot_ori,probs = 0.75),col="red")
abline(h = quantile(Tra_to_plot_ori,probs = 0.25),col="red")

Spp_to_plot_ori <- as.numeric(apply(traits_ori[,2:ncol(traits_ori)],1,sum))
plot(Spp_to_plot_ori,ylim = c(0,118),
xlab = "Species",ylab = "Traits",main = "How many traits has each species")
abline(h = mean(Spp_to_plot_ori),col="red",lwd=3)
abline(h = quantile(Spp_to_plot_ori,probs = 0.75),col="red")
abline(h = quantile(Spp_to_plot_ori,probs = 0.25),col="red")


# We create the data.frame where we will store everything during the loops
output <- data.frame()

# Generalist vs Specialist traits generator
Special_Traits <- c(c(1,10,30,50,70,90,110),c(1,10,30,50,70,90,110))
Homogen_Spp <- c(rep(30,(length(Special_Traits)/2)),rep(450,(length(Special_Traits)/2)))
for (number in c(1:length(Special_Traits))) {
  
traits <- traits_ori
Trait_packs <- sample(2:ncol(traits),118,replace = F)
Spp_packs <- sample(1:nrow(traits),460,replace = F)

Trait_Paks_seq <- seq(1,118,Special_Traits[number])
Trait_pack_matrix <- matrix(nrow = (length(Trait_Paks_seq)-1),ncol = (Special_Traits[number]),data = NA)
for (Trait_rows in 1:(length(Trait_Paks_seq)-1)) {
Trait_pack_matrix[Trait_rows,] <- Trait_packs[(Trait_Paks_seq[Trait_rows]-1):(Trait_Paks_seq[Trait_rows+1]-1)][1:Special_Traits[number]]}


Spp_Paks_seq <- seq(1,460,Homogen_Spp[number])
Spp_pack_matrix <- matrix(nrow = (length(Spp_Paks_seq)-1),ncol = Homogen_Spp[number],data = NA)
for (Spp_rows in 1:(length(Spp_Paks_seq)-1)) {
Spp_pack_matrix[Spp_rows,] <- Spp_packs[(Spp_Paks_seq[Spp_rows]-1):(Spp_Paks_seq[Spp_rows+1]-1)][1:Homogen_Spp[number]]}

if ((nrow(Trait_pack_matrix)-nrow(Spp_pack_matrix))>0) {leng_to_change <-nrow(Spp_pack_matrix)}else{
  leng_to_change <-nrow(Trait_pack_matrix)}

for (final_rows in 1:leng_to_change) {
traits[Spp_pack_matrix[final_rows,],Trait_pack_matrix[final_rows,]] <- 1
}

Tra_to_plot <- as.numeric(apply(traits[,2:ncol(traits)],2,sum))
plot(Tra_to_plot,ylim = c(0,460),
     xlab = "Traits",ylab = "Species",main = "How many species have the same trait")
abline(h = mean(Tra_to_plot),col="red",lwd=3)
abline(h = quantile(Tra_to_plot,probs = 0.75),col="red")
abline(h = quantile(Tra_to_plot,probs = 0.25),col="red")

Spp_to_plot <- as.numeric(apply(traits[,2:ncol(traits)],1,sum))
plot(Spp_to_plot,ylim = c(0,118),
     xlab = "Species",ylab = "Traits",main = "How many traits has each species")
abline(h = mean(Spp_to_plot),col="red",lwd=3)
abline(h = quantile(Spp_to_plot,probs = 0.75),col="red")
abline(h = quantile(Spp_to_plot,probs = 0.25),col="red")





# Load dataset from retromed
macros_lakes_list <- readxl::read_excel("data/Lakes/macros_lakes_list.xlsx") %>%
  filter(inv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.

# All years 
years <- unique(macros_lakes_list$year)

### FIRST LOOP - Years
for (ind_year in 1:1) {#:length(years)) {
cat("We are at year", years[ind_year],"__________________________________________________","\n")
# We create a "temporary" file filtered according to the year selected. 
macros_lakes_list_temp <- macros_lakes_list %>% filter(year==years[ind_year])
# All lakes of that year
year_lakes <- unique(macros_lakes_list_temp$site)
  
### SECOND LOOP - Number of randomly selected lakes
for (rand_selection in c(1,3,4,6,8,10,13,seq(15,length(year_lakes),10))) {
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

#run analysis (takes a bit)
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)

# Generate an adjacency table from the coocurrence matrix 
## Create a blank table with only 0 
adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {# do the following for each spp
  # Summarise the strenght of the link for each genus by calulating the mean probability
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur),.groups = "drop")
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

# Layout plot (distance between nodes based on the links weight)
l <- igraph::layout_with_fr(g,start.temp = 1500)#start.temp is what allows more or less "movement" along the axis. High values make more "sparse" graphs but also they are "far" from reality. Smaller values make more "dense" graphs that provide less visual information. Balanxce between both levels must be achieved. 
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis) # Packages

n<- network(as.matrix(igraph::as_adjacency_matrix(g)), directed=F) # Newtork creation
# Setting the parameters that we want to use. You can add as many parameters as you wish. 
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector # Weighted degree
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum) # Trait frequency
n %e% "Edg_values" <-igraph::E(g)$weight # Edge weights (coocurrence probability)

# The plot per se
ggplot(n, 
       layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+ # We need to set the layout that we created. Remember that this is just to have an x and y coordinates. 
  geom_edges(aes(color =Edg_values, size=Edg_values),arrow=arrow(angle = 20),curvature = 0.15)+ # Information and aesthetics of the edges
  geom_nodes(aes(size=FreqTRAIT, fill=Mean_Prob),color="black" ,shape=21)+ # Information and aesthetics of the edges
  # From here onwards everythign is like a normal ggplot... 
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())

#variables
edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
edge_dens <- igraph::edge_density(g)
mean_grStre <- mean(igraph::graph.strength(g))
sd_grStre <- sd(igraph::graph.strength(g))
medi_grStre <- median(igraph::graph.strength(g))
min_grStre <- min(igraph::graph.strength(g))
max_grStre <- max(igraph::graph.strength(g))

NAT_Modul <- igraph::modularity(g,membership = 1:118)

### Here we have calculated everything that we needed already, so now it is time to "close" the whole process
## we will first create a data.frame with all the information 

temp_output <- data.frame("Year"=years[ind_year],
                          "n_sites"=rand_selection,
                          "iter"=iteration,
                          "Number"=number,
                          "Traits_up_0.75"=length(which(Tra_to_plot>quantile(Tra_to_plot_ori,0.75))),
                          "Species_up_0.75"=length(which(Spp_to_plot>quantile(Spp_to_plot_ori,0.75))),
                          "Avrg_Trait_SPP"=mean(Spp_to_plot),
                          edge_dens_man,
                          edge_dens,
                          mean_grStre,
                          sd_grStre,
                          medi_grStre,
                          min_grStre, 
                          max_grStre,
                          NAT_Modul)

output <- bind_rows(output,temp_output)
}# End of iteration
}# End of rand_selection
}# End of ind_year



}# End of number


output %>%
  pivot_longer(cols = 7:ncol(.)) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=n_sites, colour=as.factor(Year)), width = 0.2)+
  geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()

gridExtra::grid.arrange(
output %>% 
  pivot_longer(cols = 8:ncol(.)) %>% 
  filter(name=="NAT_Modul") %>% 
  group_by(Year,n_sites,Number) %>% 
  mutate(Mean_val=mean(value)) %>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Avrg_Trait_SPP,2))), width = 2)+
  geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Avrg_Trait_SPP,2))))+
  scale_color_viridis(discrete = T)+
  theme_classic()+labs(colour="Avrg_Trait_SPP",title="Modularity"),

output %>% 
  pivot_longer(cols = 8:ncol(.)) %>% 
  filter(name=="edge_dens") %>% 
  group_by(Year,n_sites,Number) %>% 
  mutate(Mean_val=mean(value)) %>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Avrg_Trait_SPP,2))), width = 2)+
  geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Avrg_Trait_SPP,2))))+
  scale_color_viridis(discrete = T)+
  theme_classic()+labs(colour="Avrg_Trait_SPP",title="Edge density"),
ncol=2)

gridExtra::grid.arrange(
  output %>% 
    pivot_longer(cols = 8:ncol(.)) %>% 
    filter(name=="NAT_Modul") %>% 
    group_by(Year,n_sites,Number) %>% 
    mutate(Mean_val=mean(value)) %>% 
    ggplot()+ 
    geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Species_up_0.75,2))), width = 2)+
    geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Species_up_0.75,2))))+
    scale_color_viridis(discrete = T)+
    theme_classic()+labs(colour="Species avobe 57 traits",title="Modularity"),
  
  output %>% 
    pivot_longer(cols = 8:ncol(.)) %>% 
    filter(name=="edge_dens") %>% 
    group_by(Year,n_sites,Number) %>% 
    mutate(Mean_val=mean(value)) %>% 
    ggplot()+ 
    geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Species_up_0.75,2))), width = 2)+
    geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Species_up_0.75,2))))+
    scale_color_viridis(discrete = T)+
    theme_classic()+labs(colour="Species avobe 57 traits",title="Edge density"),
  ncol=2)


# SPECIALIST ####
# Traits database (canviar a filtrat per generes)
traits_ori <-read.csv("data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")

# Transform traits to 1 or 0 (losing affiliations)
traits_ori[,11:ncol(traits_ori)] <- ifelse(traits_ori[,11:ncol(traits_ori)]>=1,1,0)#all afiliations higher than 1 have a 1

traits_ori <- traits_ori %>% group_by(Genus..if.description.at.this.level.) %>% 
  mutate_if(is.numeric, ~mean(.)) %>% 
  mutate_if(is.numeric, ~ifelse(.<0.5,0,1)) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

sum_traits <- as.numeric(apply(traits_ori[,2:ncol(traits_ori)],2,sum))
traits_ori[which(is.na(traits_ori[,(which(is.na(sum_traits)==T)+1)])==T),(which(is.na(sum_traits)==T)+1)] <- 0


par(mfrow=c(4,8))
Tra_to_plot_ori <- as.numeric(apply(traits_ori[,2:ncol(traits_ori)],2,sum))
plot(Tra_to_plot_ori,ylim = c(0,460),
     xlab = "Traits",ylab = "Species",main = "How many species have the same trait")
abline(h = mean(Tra_to_plot_ori),col="red",lwd=3)
abline(h = quantile(Tra_to_plot_ori,probs = 0.75),col="red")
abline(h = quantile(Tra_to_plot_ori,probs = 0.25),col="red")

Spp_to_plot_ori <- as.numeric(apply(traits_ori[,2:ncol(traits_ori)],1,sum))
plot(Spp_to_plot_ori,ylim = c(0,118),
     xlab = "Species",ylab = "Traits",main = "How many traits has each species")
abline(h = mean(Spp_to_plot_ori),col="red",lwd=3)
abline(h = quantile(Spp_to_plot_ori,probs = 0.75),col="red")
abline(h = quantile(Spp_to_plot_ori,probs = 0.25),col="red")

traits <- traits_ori
traits[,2:ncol(traits)] <- 0

as.numeric(apply(traits[,2:ncol(traits)],1,sum))

# We create the data.frame where we will store everything during the loops
output <- data.frame()

groups_of_tratis <- c(10,15,20,25,30,40,45,seq(50,118,10))
# Generalist vs Specialist traits generator
for (number in c(1:length(groups_of_tratis))) {
  
start <-seq(1,nrow(traits),groups_of_tratis[number])
if(start[length(start)]!=nrow(traits)){start <-c(seq(1,nrow(traits),groups_of_tratis[number]),nrow(traits))}
sequences_list <- list()
for (grouping in 1:(length(start)-1)) {
sequences_list[[grouping]] <- seq(start[grouping],(start[grouping+1]-1),1)
}

tr_start <- round(seq(2,  (ncol(traits)-1),(ncol(traits)-1)/length(sequences_list)))
tr_sequences_list <- list()
for (grouping in 1:(length(tr_start)-1)) {
tr_sequences_list[[grouping]] <- seq(tr_start[grouping],(tr_start[grouping+1]-1),1)
} 

if ((length(sequences_list)-length(tr_sequences_list))>0) {leng_to_change <-length(tr_sequences_list)}else{
  leng_to_change <-nrow(sequences_list)}

for (final_rows in 1:leng_to_change) {
  traits[sequences_list[[final_rows]],tr_sequences_list[[final_rows]]] <- 1
}


  Tra_to_plot <- as.numeric(apply(traits[,2:ncol(traits)],2,sum))
  plot(Tra_to_plot,ylim = c(0,460),
       xlab = "Traits",ylab = "Species",main = "How many species have the same trait")
  abline(h = mean(Tra_to_plot),col="red",lwd=3)
  abline(h = quantile(Tra_to_plot,probs = 0.75),col="red")
  abline(h = quantile(Tra_to_plot,probs = 0.25),col="red")
  
  Spp_to_plot <- as.numeric(apply(traits[,2:ncol(traits)],1,sum))
  plot(Spp_to_plot,ylim = c(0,118),
       xlab = "Species",ylab = "Traits",main = "How many traits has each species")
  abline(h = mean(Spp_to_plot),col="red",lwd=3)
  abline(h = quantile(Spp_to_plot,probs = 0.75),col="red")
  abline(h = quantile(Spp_to_plot,probs = 0.25),col="red")
  
  
  
  
  
  # Load dataset from retromed
  macros_lakes_list <- readxl::read_excel("data/Lakes/macros_lakes_list.xlsx") %>%
    filter(inv.l>0) # There are some lakes with a 0 abundance for the taxons. Must be removed.
  
  # All years 
  years <- unique(macros_lakes_list$year)
  
  ### FIRST LOOP - Years
  for (ind_year in 1:1) {#:length(years)) {
    cat("We are at year", years[ind_year],"__________________________________________________","\n")
    # We create a "temporary" file filtered according to the year selected. 
    macros_lakes_list_temp <- macros_lakes_list %>% filter(year==years[ind_year])
    # All lakes of that year
    year_lakes <- unique(macros_lakes_list_temp$site)
    
    ### SECOND LOOP - Number of randomly selected lakes
    for (rand_selection in c(3,4,6,8,10,13,seq(15,length(year_lakes),10))) {
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
        
        #run analysis (takes a bit)
        cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
        
        # Generate an adjacency table from the coocurrence matrix 
        ## Create a blank table with only 0 
        adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
        for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {# do the following for each spp
          # Summarise the strenght of the link for each genus by calulating the mean probability
          link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur),.groups = "drop")
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
        
        # Layout plot (distance between nodes based on the links weight)
        l <- igraph::layout_with_fr(g,start.temp = 1500)#start.temp is what allows more or less "movement" along the axis. High values make more "sparse" graphs but also they are "far" from reality. Smaller values make more "dense" graphs that provide less visual information. Balanxce between both levels must be achieved. 
        colnames(l) <- c("x","y")
        library(network);library(ggnetwork);library(viridis) # Packages
        
        n<- network(as.matrix(igraph::as_adjacency_matrix(g)), directed=F) # Newtork creation
        # Setting the parameters that we want to use. You can add as many parameters as you wish. 
        n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector # Weighted degree
        n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum) # Trait frequency
        n %e% "Edg_values" <-igraph::E(g)$weight # Edge weights (coocurrence probability)
        
        # The plot per se
        ggplot(n, 
               layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+ # We need to set the layout that we created. Remember that this is just to have an x and y coordinates. 
          geom_edges(aes(color =Edg_values, size=Edg_values),arrow=arrow(angle = 20),curvature = 0.15)+ # Information and aesthetics of the edges
          geom_nodes(aes(size=FreqTRAIT, fill=Mean_Prob),color="black" ,shape=21)+ # Information and aesthetics of the edges
          # From here onwards everythign is like a normal ggplot... 
          scale_color_viridis(direction = -1,option = "D")+
          scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
          theme_classic()+
          theme(axis.text = element_blank(),
                axis.ticks = element_blank(),
                panel.background=element_blank())
        
        #variables
        edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
        edge_dens <- igraph::edge_density(g)
        mean_grStre <- mean(igraph::graph.strength(g))
        sd_grStre <- sd(igraph::graph.strength(g))
        medi_grStre <- median(igraph::graph.strength(g))
        min_grStre <- min(igraph::graph.strength(g))
        max_grStre <- max(igraph::graph.strength(g))
        
        NAT_Modul <- igraph::modularity(g,membership = 1:118)
        
        ### Here we have calculated everything that we needed already, so now it is time to "close" the whole process
        ## we will first create a data.frame with all the information 
        
        temp_output <- data.frame("Year"=years[ind_year],
                                  "n_sites"=rand_selection,
                                  "iter"=iteration,
                                  "Trait_grouping"=groups_of_tratis[number],
                                  "Traits_up_0.75"=length(which(Tra_to_plot>quantile(Tra_to_plot_ori,0.75))),
                                  "Species_up_0.75"=length(which(Spp_to_plot>quantile(Spp_to_plot_ori,0.75))),
                                  "Avrg_Trait_SPP"=mean(Spp_to_plot),
                                  edge_dens_man,
                                  edge_dens,
                                  mean_grStre,
                                  sd_grStre,
                                  medi_grStre,
                                  min_grStre, 
                                  max_grStre,
                                  NAT_Modul)
        
        output <- bind_rows(output,temp_output)
      }# End of iteration
    }# End of rand_selection
  }# End of ind_year
  
  
  
}# End of number


output %>%
  pivot_longer(cols = 7:ncol(.)) %>% 
  group_by(Year,n_sites,name) %>% 
  mutate(Mean_val=mean(value)) %>% 
  ggplot()+ 
  geom_jitter(aes(y=value, x=n_sites, colour=as.factor(Year)), width = 0.2)+
  geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(Year)))+
  facet_grid(name~.,scales = "free") + 
  theme_classic()

gridExtra::grid.arrange(
  output %>% 
    pivot_longer(cols = 8:ncol(.)) %>% 
    filter(name=="NAT_Modul") %>% 
    group_by(Year,n_sites,Trait_grouping) %>% 
    mutate(Mean_val=mean(value)) %>% 
    ggplot()+ 
    geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Trait_grouping,2))), width = 2)+
    geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Trait_grouping,2))))+
    scale_color_viridis(discrete = T)+
    theme_classic()+labs(colour="N_Trait_PerGroup",title="Modularity"),
  
  output %>% 
    pivot_longer(cols = 8:ncol(.)) %>% 
    filter(name=="edge_dens") %>% 
    group_by(Year,n_sites,Trait_grouping) %>% 
    mutate(Mean_val=mean(value)) %>% 
    ggplot()+ 
    geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Trait_grouping,2))), width = 2)+
    geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Trait_grouping,2))))+
    scale_color_viridis(discrete = T)+
    theme_classic()+labs(colour="N_Trait_PerGroup",title="Edge density"),
  ncol=2)

gridExtra::grid.arrange(
  output %>% 
    pivot_longer(cols = 8:ncol(.)) %>% 
    filter(name=="NAT_Modul") %>% 
    group_by(Year,n_sites,Trait_grouping) %>% 
    mutate(Mean_val=mean(value)) %>% 
    ggplot()+ 
    geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Species_up_0.75,2))), width = 2)+
    geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Species_up_0.75,2))))+
    scale_color_viridis(discrete = T)+
    theme_classic()+labs(colour="Species avobe 57 traits",title="Modularity"),
  
  output %>% 
    pivot_longer(cols = 8:ncol(.)) %>% 
    filter(name=="edge_dens") %>% 
    group_by(Year,n_sites,Trait_grouping) %>% 
    mutate(Mean_val=mean(value)) %>% 
    ggplot()+ 
    geom_jitter(aes(y=value, x=n_sites, colour=as.factor(round(Species_up_0.75,2))), width = 2)+
    geom_line(aes(y=Mean_val, x=n_sites, colour=as.factor(round(Species_up_0.75,2))))+
    scale_color_viridis(discrete = T)+
    theme_classic()+labs(colour="Species avobe 57 traits",title="Edge density"),
  ncol=2)



