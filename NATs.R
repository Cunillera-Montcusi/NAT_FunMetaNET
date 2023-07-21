
source("function_to_NATs.R")

library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("data/Retromed/abun_macro_sp.RData")
traits <-read.csv("data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l.Total <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
  mutate(Network = case_when(
    str_detect(codi_lloc,"CLOT")~ "GUILS",
    str_detect(codi_lloc,"CANB")~ "ALBERA",
    str_detect(codi_lloc,"CANG")~ "ALBERA",
    str_detect(codi_lloc,"CANP")~ "ALBERA",
    str_detect(codi_lloc,"DELF")~ "ALBERA",
    str_detect(codi_lloc,"FAIG")~ "ALBERA",
    str_detect(codi_lloc,"GUTC")~ "ALBERA",
    str_detect(codi_lloc,"GUTR")~ "ALBERA",
    str_detect(codi_lloc,"SERR")~ "ALBERA",
    str_detect(codi_lloc,"TORC")~ "ALBERA",
    str_detect(codi_lloc,"TORG")~ "ALBERA",
    str_detect(codi_lloc,"TORP")~ "ALBERA",
    str_detect(codi_lloc,"B1")~ "GIARA",
    str_detect(codi_lloc,"B2")~ "GIARA",
    str_detect(codi_lloc,"B3")~ "GIARA",
    str_detect(codi_lloc,"E1")~ "GIARA",
    str_detect(codi_lloc,"E3")~ "GIARA",
    str_detect(codi_lloc,"M4")~ "GIARA",
    str_detect(codi_lloc,"M5")~ "GIARA",
    str_detect(codi_lloc,"M7")~ "GIARA",
    str_detect(codi_lloc,"S4")~ "GIARA",
    str_detect(codi_lloc,"S5")~ "GIARA",
    str_detect(codi_lloc,"S6")~ "GIARA",
    TRUE ~ "altres" )) 

# We isolate the values of each network
Network_names <- unique(indv.l.Total$Network)[1:3]
# We create the main table where we will store all the results
out_table_Networks <- data.frame()
# We create the list where we will store all the plots
plot_Networks_Total <- list()
for (Net_Net in 1:length(Network_names)) {
# We define the range of samples that we will add. From 1 to 10 in this case. This will define our "loops"
number_additions <- 1:10

# We filter for the corresponding network
indv.l.filt <- indv.l.Total%>%filter(Network==Network_names[[Net_Net]])

# We select the pond names from the network
Pond_names <- unique(as.vector(indv.l.filt %>% ungroup()%>% select(codi_lloc))[[1]])
# We creat the table to store the results
out_table <- data.frame()
plot_Networks <- list()
for (Add_samp in 1:length(number_additions)) {# We set a loop to do this action as many times as the "number of additions" 
for (iterati in 1:15) {# We set another loop to repeat the same action "x" number of times. 
# With these two loops we will do the following actions for each element of the "number_additions" several times.
## In other words: We will do NATs for only 1 sample 10 times, for 2 samples 10 times, etc ... 

# The function sample "samples" elements of a determined vector randomly. We use it to randomly select samples
Selecte_ponds <- sample(Pond_names,size = number_additions[Add_samp],replace = F) 

# We filter the dataframe to obtain only the species present in the ponds selected randomly.
indv.l.filt.filt <- indv.l.filt %>% filter(codi_lloc%in%Selecte_ponds)

# We now can create the matrix of spp x traits that we need for the calculation of the NAT
## WARNING: All the following steps should not be necessary in a "good" dataset. The current dataset is not 
## fully depurated and contains some mistakes (names repeated or other uncoherent things) that have not been 
## corrected as it is only intended to exemplify. 
traits_ind <- indv.l.filt.filt %>% ungroup() %>% # We remove all previous groupings
                                   group_by(genus) %>% # We group by genus
                                   summarise(tot_ab=sum(indv.l)) %>% # We sum all the abundances to avoid 
                                   # repeated species. 
                                   left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),
                                             multiple ="all") %>% # We merge the dataset of the traits. 
                                   # We merge the two databases based on their Genus names, any mismatch is due 
                                   # differences between the two databases names. 
                                   group_by(genus) %>% # We group again
                                   filter(n()<2) %>% # We filter some repeated columns with the same information
                                   dplyr::select(c(1,12:ncol(.))) %>% # We select only the columns of interest 
                                   na.omit() %>% # Eliminate (just in case) the NAs
                                   tibble::column_to_rownames("genus") # We attach the genus into the rownames


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


# We now save all these values into a dataframe that we will use to store all the information generated
## along the loop
out_temp <- data.frame("edge_dens"=edge_dens, # Result
                       "mean_grStre"=mean_grStre, # Result
                       "sd_grStre"=sd_grStre, # Result
                       "medi_grStre"=medi_grStre, # Result
                       "min_grStre"=min_grStre, # Result
                       "max_grStre"=max_grStre, # Result
                       "Network"=Network_names[[Net_Net]], # Information about which network we are using
                       "AddedPonds"=length(Selecte_ponds), # Information about how many ponds we used
                       "Iter"=iterati) 

# We attach this results to the previous loop. 
out_table <- bind_rows(out_table,out_temp)

#______________________
## ALL THIS CHUNK is just to create a plot and store it in a list. If we do not want plots we can just remove it. 
#______________________
# Layout plot (distance between nodes based on the links weight)
l <- igraph::layout_with_fr(NATs_Output$Graph,start.temp = 1500)#start.temp is what allows more or less "movement" along the axis. High values make more "sparse" graphs but also they are "far" from reality. Smaller values make more "dense" graphs that provide less visual information. Balanxce between both levels must be achieved. 
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis) # Packages
n<- network(as.matrix(igraph::as_adjacency_matrix(NATs_Output$Graph)), directed=F) # Newtork creation
# Setting the parameters that we want to use. You can add as many parameters as you wish. 
n %v% "Mean_Prob" <- igraph::eigen_centrality(NATs_Output$Graph)$vector # Weighted degree
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum) # Trait frequency
n %e% "Edg_values" <-igraph::E(NATs_Output$Graph)$weight # Edge weights (coocurrence probability)

# The plot per se
plot_Networks[[Add_samp]] <- ggplot(n,
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
#______________________



}# Iterations loop end
}# Add_samp loop end
# We attach the results to the main table where we store all results for each network
out_table_Networks <- bind_rows(out_table_Networks,out_table)
plot_Networks_Total[[Net_Net]] <- plot_Networks
}# Net_Net loop end 

#_________________________
# AND BIM! The results are all stored here:
#_________________________
out_table_Networks
# And all plots in plot_Networks_Total


#_________________________
#_________________________
## Let's continue with the GAMMA NAT
# We now carry the "GAMMA" NAT considering the 3 newtorks together
indv.l.filt.filt <- indv.l.Total%>%filter(Network%in%Network_names)

# We now can create the matrix of spp x traits that we need for the calculation of the NAT
## WARNING: All the following steps should not be necessary in a "good" dataset. The current dataset is not 
## fully depurated and contains some mistakes (names repeated or other uncoherent things) that have not been 
## corrected as it is only intended to exemplify. 
traits_ind <- indv.l.filt.filt %>% ungroup() %>% # We remove all previous groupings
  group_by(genus) %>% # We group by genus
  summarise(tot_ab=sum(indv.l)) %>% # We sum all the abundances to avoid 
  # repeated species. 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),
            multiple ="all") %>% # We merge the dataset of the traits. 
  # We merge the two databases based on their Genus names, any mismatch is due 
  # differences between the two databases names. 
  group_by(genus) %>% # We group again
  filter(n()<2) %>% # We filter some repeated columns with the same information
  dplyr::select(c(1,12:ncol(.))) %>% # We select only the columns of interest 
  na.omit() %>% # Eliminate (just in case) the NAs
  tibble::column_to_rownames("genus") # We attach the genus into the rownames


NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind )

edge_dens <- igraph::edge_density(NATs_Output$Graph)
mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
sd_grStre <- sd(igraph::graph.strength(NATs_Output$Graph))
medi_grStre <- median(igraph::graph.strength(NATs_Output$Graph))
min_grStre <- min(igraph::graph.strength(NATs_Output$Graph))
max_grStre <- max(igraph::graph.strength(NATs_Output$Graph))


out_GAMMA_NAT <- data.frame("edge_dens"=edge_dens, # Result
                       "mean_grStre"=mean_grStre, # Result
                       "sd_grStre"=sd_grStre, # Result
                       "medi_grStre"=medi_grStre, # Result
                       "min_grStre"=min_grStre, # Result
                       "max_grStre"=max_grStre # Result
                       )

#_________________________
#_________________________
## Let's continue with the MEGA_GAMMA NAT
# We now carry the "MEGA_GAMMA" NAT considering the the database! 

# The matrix is actually the tachet database... so no need of merging! But again, need of checking for
# duplicates and stuff, so far we eliminate them for pragmatism. 
traits_ind <-traits %>% 
  dplyr::select(c(4,11:ncol(.))) %>%  # Select the traits and the "genus" column
  group_by(Genus..if.description.at.this.level.) %>% # Group by genus
  summarise_if(.predicate = is.numeric,.funs = sum) %>% # We sum repeated traits... probably larvae
  na.omit() %>%
  tibble::column_to_rownames("Genus..if.description.at.this.level.") # We attach the genus into the rownames

traits_ind <- ifelse(traits_ind>=1,1,0)

NATs_Output <- fun_to_NATs(Spp_x_Traits_Matrix =traits_ind )

edge_dens <- igraph::edge_density(NATs_Output$Graph)
mean_grStre <- mean(igraph::graph.strength(NATs_Output$Graph))
sd_grStre <- sd(igraph::graph.strength(NATs_Output$Graph))
medi_grStre <- median(igraph::graph.strength(NATs_Output$Graph))
min_grStre <- min(igraph::graph.strength(NATs_Output$Graph))
max_grStre <- max(igraph::graph.strength(NATs_Output$Graph))


out_MEGA_GAMMA_NAT <- data.frame("edge_dens"=edge_dens, # Result
                            "mean_grStre"=mean_grStre, # Result
                            "sd_grStre"=sd_grStre, # Result
                            "medi_grStre"=medi_grStre, # Result
                            "min_grStre"=min_grStre, # Result
                            "max_grStre"=max_grStre # Result
                            )

# TIME TO PLOT! 
out_table_Networks
out_GAMMA_NAT
out_MEGA_GAMMA_NAT

gridExtra::grid.arrange(
  
  out_table_Networks %>% ggplot(aes(y=edge_dens,x=AddedPonds, colour=Network))+
    geom_point()+ geom_smooth(method = "loess", se=F)+
    geom_hline(yintercept = out_GAMMA_NAT$edge_dens, colour="red")+
    geom_hline(yintercept = out_MEGA_GAMMA_NAT$edge_dens, colour="blue")+
    theme_classic()+
    labs(y="Igraph EdgeDensity",title= "Igr EdgDens"),
  
  out_table_Networks %>% ggplot(aes(y=mean_grStre,x=AddedPonds, colour=Network))+
    geom_point()+ geom_smooth(method = "loess", se=F)+
    geom_hline(yintercept = out_GAMMA_NAT$mean_grStre, colour="red")+
    geom_hline(yintercept = out_MEGA_GAMMA_NAT$mean_grStre, colour="blue")+
    theme_classic()+
    labs(y="Mean graph srtength",title= "Mean GrStre"),
  
  out_table_Networks %>% ggplot(aes(y=sd_grStre,x=AddedPonds, colour=Network))+
    geom_point()+ geom_smooth(method = "loess", se=F)+
    geom_hline(yintercept = out_GAMMA_NAT$sd_grStre, colour="red")+
    geom_hline(yintercept = out_MEGA_GAMMA_NAT$sd_grStre, colour="blue")+
    theme_classic()+
    labs(y="sd graph srtength",title= "SD GrStre"),
  
  out_table_Networks %>% ggplot(aes(y=medi_grStre,x=AddedPonds, colour=Network))+
    geom_point()+ geom_smooth(method = "loess", se=F)+
    geom_hline(yintercept = out_GAMMA_NAT$medi_grStre, colour="red")+
    geom_hline(yintercept = out_MEGA_GAMMA_NAT$medi_grStre, colour="blue")+
    theme_classic()+
    labs(y="Median graph srtength",title= "Med GrStre"),
  
  out_table_Networks %>% ggplot(aes(y=min_grStre,x=AddedPonds, colour=Network))+
    geom_point()+ geom_smooth(method = "loess", se=F)+
    geom_hline(yintercept = out_GAMMA_NAT$min_grStre, colour="red")+
    geom_hline(yintercept = out_MEGA_GAMMA_NAT$min_grStre, colour="blue")+
    theme_classic()+
    labs(y="Min graph srtength",title= "Min GrStre"),
  
  out_table_Networks %>% ggplot(aes(y=max_grStre,x=AddedPonds, colour=Network))+
    geom_point()+ geom_smooth(method = "loess", se=F)+
    geom_hline(yintercept = out_GAMMA_NAT$max_grStre, colour="red")+
    geom_hline(yintercept = out_MEGA_GAMMA_NAT$max_grStre, colour="blue")+
    theme_classic()+
    labs(y="Max graph srtength",title= "Max GrStre"),
  
  nrow=3)






