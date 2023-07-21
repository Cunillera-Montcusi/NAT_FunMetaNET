#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
# WARNING !!!!
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________

# This is a very dirty and drafted script that should only be used by David Cunillera-Montcusí
# or someone able to understand his messy existential moves... 

# Even for him might be complex.... there are other scripts more polished and nice. 
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________
#_____________________________________________________

### TRAITS DATABASE ####

# Install biomonitoR
library(devtools)
install_github("alexology/biomonitoR", ref = "main", build_vignettes = TRUE)
library(biomonitoR)

# We need to check the matching between traits and our databases
traits <-read.csv("data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
# Transform traits to 1 or 0 (losing affiliations)
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

library(cooccur);library(tidyverse);library(qgraph)
library(network);library(ggnetwork);library(viridis);library(igraph)
#LOAD 2Coocurrence - species.csv
#load dataset
load("data/abun_macro_sp.RData")


# Check presence and filter the dataset according the existing genus
# Changing manually names and assigning the names to networks
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network=="ALBERA")#, codi_lloc== "SERR")

# Generating the trait presence data (spp x trait) to later run the coocurrence
traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

# Run coocurrence analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
class(cooccur.species_ccr)
summary(cooccur.species_ccr)
#all pairs probabilities analysis
prob.table(cooccur.species_ccr)

# Generate an adjacency table from the coocurrence matrix 
## Create a blank table with only 0 
adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {# do the following for each spp
# Summarise the strenght of the link for each genus by calulating the mean probability
link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
# filter the results from the coocurrence table for a species. You isolate all the species with who the 
## targeted species is coocuring
spp_to <- cooccur.species_ccr$results %>% 
          filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
# You add the values into the "adj.table" that was empty and now is being filled
adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
}
# We transpose the results to make the matrix symetrical (same from & to links)
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]

# We built a graph from the adjacency table
g <- igraph::graph_from_adjacency_matrix(adj_table,weighted = T)


# Plots

#Layout type 1 
# Layout plot (distància entre nodes=weight dels nodes)
minC <- rep(-Inf, length(igraph::V(g)))
maxC <- rep(Inf, length(igraph::V(g)))
minC[1] <- maxC[1] <- 0

l <- igraph::layout_with_fr(g,grid = "auto",  start.temp = 20,
                            minx = minC, maxx = maxC,
                            miny = minC, maxy = maxC,niter = 12)
colnames(l) <- c("x","y")

n<- network(as.matrix(igraph::as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-igraph::E(g)$weight

a <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15) +
  geom_nodes(aes(size=Mean_Prob, fill=FreqTRAIT),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())


# Layout plot (distància entre nodes=weight dels nodes)
l <- qgraph::qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=8*(vcount(g)^2),
                                               repulse.rad=(vcount(g)^3.1))

colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-E(g)$weight

b <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
  geom_nodes(aes(size=Mean_Prob, fill=FreqTRAIT),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())


png(filename = "NAts.png", width = 3000, height = 3000, res=300)
gridExtra::grid.arrange(a,b)
dev.off()

# Closennesss
#value_CANB
#values_TORP <- igraph::closeness(g)
#values_CLOTS15 <- igraph::closeness(g)
#values_CLOTS17 <- igraph::closeness(g)
#values_SERR <- igraph::closeness(g)
#values_DELF <- igraph::closeness(g)

values_CLOTS6 <- igraph::eigen_centrality(g)$vector

gridExtra::grid.arrange(a,b)





summary(lm(
values_CLOTS6~
apply(ifelse(traits_ind>0,1,0),2,sum)
))

data.frame(
"Clo"=values_CLOTS6,
"TraitFre"=apply(ifelse(traits_ind>0,1,0),2,sum)) %>% 
ggplot()+geom_point(aes(x=Clo,y=TraitFre))+
  geom_smooth(aes(x=Clo,y=TraitFre),method = "lm")+theme_classic()


plot(g)
V(g)
colnames(traits_ind)


gridExtra::grid.arrange(
bind_rows(
data.frame("Clo"=value_CANB,"ID"="CANB","Index"=seq(1:length(value_CANB))), 
data.frame("Clo"=values_TORP,"ID"="TORP","Index"=seq(1:length(values_TORP))),
data.frame("Clo"=values_SERR,"ID"="SERR","Index"=seq(1:length(values_SERR))),
data.frame("Clo"=values_DELF,"ID"="DELF","Index"=seq(1:length(values_DELF))),
data.frame("Clo"=values_CLOTS15,"ID"="CLOTS15","Index"=seq(1:length(values_CLOTS15))),
data.frame("Clo"=values_CLOTS6,"ID"="CLOTS6","Index"=seq(1:length(values_CLOTS6))),
data.frame("Clo"=values_CLOTS17,"ID"="CLOTS17","Index"=seq(1:length(values_CLOTS17)))
) %>% 
  filter(Clo<0.5)%>%ggplot()+geom_point(aes(x=Index, y=Clo, colour=ID),alpha=0.3)+theme_classic()+
  geom_smooth(aes(x=Index, y=Clo, colour=ID),method = "lm",se=F),

bind_rows(
  data.frame("Clo"=value_CANB,"ID"="CANB","Index"=seq(1:length(value_CANB))), 
  data.frame("Clo"=values_TORP,"ID"="TORP","Index"=seq(1:length(values_TORP))),
  data.frame("Clo"=values_SERR,"ID"="SERR","Index"=seq(1:length(values_SERR))),
  data.frame("Clo"=values_DELF,"ID"="DELF","Index"=seq(1:length(values_DELF))),
  data.frame("Clo"=values_CLOTS15,"ID"="CLOTS15","Index"=seq(1:length(values_CLOTS15))),
  data.frame("Clo"=values_CLOTS6,"ID"="CLOTS6","Index"=seq(1:length(values_CLOTS6))),
  data.frame("Clo"=values_CLOTS17,"ID"="CLOTS17","Index"=seq(1:length(values_CLOTS17)))
) %>%
  filter(Clo<0.5) %>% ggplot()+geom_violin(aes(x=ID, y=Clo, fill=ID))+theme_classic()
)


bind_rows(value_CANB,values_TORP)



igraph::components(g)

which(igraph::components(g)$membership==1)

igraph::modularity(igraph::decompose.graph(g)[[1]],membership=1:113)



g <- graph_ALBERA
igraph::modularity(igraph::decompose.graph(g)[[1]],membership=1:113)

g <- graph_GUILS
igraph::modularity(igraph::decompose.graph(g)[[2]],membership=1:107)




V(g)


igraph::closeness(g)
igraph::modularity(g)









setwd("C:/Users/David CM/Desktop")

library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("abun_macro_sp.RData")
traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network=="ALBERA", codi_lloc== c("SERR"))

traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

#run analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
class(cooccur.species_ccr)
summary(cooccur.species_ccr)
#all pairs probabilities analysis
#sink("pair_table.txt")
prob.table(cooccur.species_ccr)
#sink()

adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
  
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  for (col_ID in 1:unique(spp_to$sp2)) {
    
    weight <- as.numeric(link_weight %>% 
                           filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID],sp2==unique(spp_to$sp2)[col_ID]) %>% ungroup() %>%  dplyr::select(mean_prob))
    adj_table[unique(cooccur.species_ccr$results$sp1)[row_ID],unique(spp_to$sp2)[col_ID]] <- weight
    
  }
  
}
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
adj_table <- ifelse(adj_table<0.3,0,adj_table)

g <- graph_from_adjacency_matrix(adj_table,weighted = T)


graph_SERR<- g

# Layout plot (distància entre nodes=weight dels nodes)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-E(g)$weight

plot_SERR <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
  geom_nodes(aes(size=Mean_Prob, fill=FreqTRAIT),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())



setwd("C:/Users/David CM/Desktop")

library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("abun_macro_sp.RData")
traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network=="ALBERA", codi_lloc== "CANB")

traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

#run analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
class(cooccur.species_ccr)
summary(cooccur.species_ccr)
#all pairs probabilities analysis
#sink("pair_table.txt")
prob.table(cooccur.species_ccr)
#sink()

adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
  
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  for (col_ID in 1:unique(spp_to$sp2)) {
    
    weight <- as.numeric(link_weight %>% 
                           filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID],sp2==unique(spp_to$sp2)[col_ID]) %>% ungroup() %>%  dplyr::select(mean_prob))
    adj_table[unique(cooccur.species_ccr$results$sp1)[row_ID],unique(spp_to$sp2)[col_ID]] <- weight
    
  }
  
}
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
adj_table <- ifelse(adj_table<0.3,0,adj_table)

g <- graph_from_adjacency_matrix(adj_table,weighted = T)

graph_CANB <- g

# Layout plot (distància entre nodes=weight dels nodes)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-E(g)$weight

plot_CANB <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
  geom_nodes(aes(size=Mean_Prob, fill=FreqTRAIT),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())



setwd("C:/Users/David CM/Desktop")

library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("abun_macro_sp.RData")
traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network=="ALBERA", codi_lloc== "DELF")

traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

#run analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
class(cooccur.species_ccr)
summary(cooccur.species_ccr)
#all pairs probabilities analysis
#sink("pair_table.txt")
prob.table(cooccur.species_ccr)
#sink()

adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
  
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  for (col_ID in 1:unique(spp_to$sp2)) {
    
    weight <- as.numeric(link_weight %>% 
                           filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID],sp2==unique(spp_to$sp2)[col_ID]) %>% ungroup() %>%  dplyr::select(mean_prob))
    adj_table[unique(cooccur.species_ccr$results$sp1)[row_ID],unique(spp_to$sp2)[col_ID]] <- weight
    
  }
  
}
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
adj_table <- ifelse(adj_table<0.3,0,adj_table)

g <- graph_from_adjacency_matrix(adj_table,weighted = T)

graph_DELF <- g

# Layout plot (distància entre nodes=weight dels nodes)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-E(g)$weight

plot_DELF <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
  geom_nodes(aes(size=Mean_Prob, fill=FreqTRAIT),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())



setwd("C:/Users/David CM/Desktop")

library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("abun_macro_sp.RData")
traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network=="ALBERA", codi_lloc== "GUTC")

traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

#run analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
class(cooccur.species_ccr)
summary(cooccur.species_ccr)
#all pairs probabilities analysis
#sink("pair_table.txt")
prob.table(cooccur.species_ccr)
#sink()

adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
  
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  for (col_ID in 1:unique(spp_to$sp2)) {
    
    weight <- as.numeric(link_weight %>% 
                           filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID],sp2==unique(spp_to$sp2)[col_ID]) %>% ungroup() %>%  dplyr::select(mean_prob))
    adj_table[unique(cooccur.species_ccr$results$sp1)[row_ID],unique(spp_to$sp2)[col_ID]] <- weight
    
  }
  
}
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
adj_table <- ifelse(adj_table<0.3,0,adj_table)

g <- graph_from_adjacency_matrix(adj_table,weighted = T)

graph_GUTC <- g

# Layout plot (distància entre nodes=weight dels nodes)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-E(g)$weight

plot_GUTC <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
  geom_nodes(aes(size=Mean_Prob, fill=FreqTRAIT),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())


setwd("C:/Users/David CM/Desktop")

library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("abun_macro_sp.RData")
traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network=="ALBERA")#, codi_lloc== "SERR")

traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

#run analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
class(cooccur.species_ccr)
summary(cooccur.species_ccr)
#all pairs probabilities analysis
#sink("pair_table.txt")
prob.table(cooccur.species_ccr)
#sink()

adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
  
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  for (col_ID in 1:unique(spp_to$sp2)) {
    
    weight <- as.numeric(link_weight %>% 
                           filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID],sp2==unique(spp_to$sp2)[col_ID]) %>% ungroup() %>%  dplyr::select(mean_prob))
    adj_table[unique(cooccur.species_ccr$results$sp1)[row_ID],unique(spp_to$sp2)[col_ID]] <- weight
    
  }
  
}
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
adj_table <- ifelse(adj_table<0.3,0,adj_table)

g <- graph_from_adjacency_matrix(adj_table,weighted = T)

graph_ALBERA <- g


# Layout plot (distància entre nodes=weight dels nodes)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=10*(vcount(g)^2.7),repulse.rad=(vcount(g)^3.5))
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-E(g)$weight

plot_ALBERA <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1,) +
  geom_nodes(aes(size=Mean_Prob, fill=FreqTRAIT),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())


gridExtra::grid.arrange(
  plot_CANB+labs(title="CANB",subtitle=paste("nº links",length(E(graph_CANB)))),
  plot_DELF+labs(title="DELF",subtitle=paste("nº links",length(E(graph_DELF)))),
  plot_GUTC+labs(title="GUTC",subtitle=paste("nº links",length(E(graph_GUTC)))),
  plot_SERR+labs(title="SERR",subtitle=paste("nº links",length(E(graph_SERR)))),
  plot_ALBERA+labs(title="ALBERA",subtitle=paste("nº links",length(E(graph_ALBERA))))
  )


igraph::edge_density(graph_CANB)
length(E(graph_CANB))/(length(V(graph_CANB))*length(V(graph_CANB)))





igraph::edge_density(graph_DELF)
igraph::edge_density(graph_GUTC)
igraph::edge_density(graph_SERR)
igraph::edge_density(graph_ALBERA)

igraph::degree(graph_ALBERA)*graph.strength(graph_ALBERA)

E(graph_CANB)$weight







pond_name <- c("CANB","CANG","CANP","DELF","FAIG","GUTC","GUTR","SERR","TORC","TORG","TORP")
number_additions <- 1:10
plot_Networks <- list()
out_table <- data.frame()

for (Net_Net in 1:length(number_additions)) {
# Falta el loop per fer diverses vegades
for (iterati in 1:8) {
Selecte_ponds <- sample(pond_name,number_additions[Net_Net],replace = F) 
  
setwd("C:/Users/David CM/Desktop")
library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("abun_macro_sp.RData")
traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network=="ALBERA", codi_lloc%in%Selecte_ponds)

traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

#run analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)

adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
  
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  for (col_ID in 1:unique(spp_to$sp2)) {
    
    weight <- as.numeric(link_weight %>% 
                           filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID],sp2==unique(spp_to$sp2)[col_ID]) %>% ungroup() %>%  dplyr::select(mean_prob))
    adj_table[unique(cooccur.species_ccr$results$sp1)[row_ID],unique(spp_to$sp2)[col_ID]] <- weight
  }
}
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
#adj_table <- ifelse(adj_table<0.3,0,adj_table)

g <- graph_from_adjacency_matrix(adj_table,weighted = T)

# Layout plot (distància entre nodes=weight dels nodes)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=10*(vcount(g)^2.7),repulse.rad=(vcount(g)^3.5))
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-E(g)$weight

plot_Networks[[Net_Net]] <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
  geom_nodes(aes(size=FreqTRAIT, fill=Mean_Prob),color="black" ,shape=21)+
  #geom_label(aes(label=Names),label.size = 0.1)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())


edge_dens <- length(E(g))/(ncol(traits_ind)*ncol(traits_ind))
out_temp <- data.frame("EdG"=edge_dens,"AddedPonds"=length(Selecte_ponds),"Iter"=iterati)
out_table <- bind_rows(out_table,out_temp)

}
}


plot_Networks[[1]]
plot_Networks[[10]]


out_table %>% mutate(link=EdG*(ncol(traits_ind)*ncol(traits_ind)))



out_table %>% ggplot()+geom_point(aes(y=EdG,x=AddedPonds))+
  geom_smooth(aes(y=EdG,x=AddedPonds), method = "loess", se=F, colour="grey40")+
  geom_hline(yintercept = length(E(g_1))/(ncol(traits_ind)*ncol(traits_ind)), colour="red")+
  scale_y_continuous(limits = c(0,0.5))+
  theme_classic()

E(graph_ALBERA)

igraph::edge_density(g_1)
(length(E(g_1))/((118*117)/2))/2



# ACCROSS SYSTEMS! 

pond_name <- c("CANB","CANG","CANP","DELF","FAIG","GUTC","GUTR","SERR","TORC","TORG","TORP")
number_additions <- 1:10
plot_Networks <- list()
out_table <- data.frame()

for (Net_Net in 1:length(number_additions)) {
  # Falta el loop per fer diverses vegades
  for (iterati in 1:10) {
    Selecte_ponds <- sample(pond_name,number_additions[Net_Net],replace = F) 
    
    setwd("C:/Users/David CM/Desktop")
    library(cooccur);library(tidyverse)
    #LOAD 2Coocurrence - species.csv
    #load dataset
    load("abun_macro_sp.RData")
    traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
    traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)
    
    # Check presence
    indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
        TRUE ~ "altres" )) %>% 
      filter(Network=="ALBERA", codi_lloc%in%Selecte_ponds)
    
    traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
      left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
      group_by(genus) %>% filter(n()<2) %>% 
      select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")
    
    #run analysis
    cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
    
    adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
    for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
      link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
      spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
      adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
    }
    adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
    
    g <- igraph::graph_from_adjacency_matrix(adj_table,weighted = T)
    
    # Layout plot (distància entre nodes=weight dels nodes)
    minC <- rep(-Inf, length(igraph::V(g)))
    maxC <- rep(Inf, length(igraph::V(g)))
    minC[1] <- maxC[1] <- 0
    
    l <- igraph::layout_with_fr(g,grid = "auto",  start.temp = 30,
                                minx = minC, maxx = maxC,
                                miny = minC, maxy = maxC,niter = 12)
    colnames(l) <- c("x","y")
    library(network);library(ggnetwork);library(viridis)
    n<- network(as.matrix(igraph::as_adjacency_matrix(g)), directed=F, diag=T)
    n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
    n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
    n %v% "Names" <-colnames(traits_ind)
    n %e% "Edg_values" <-igraph::E(g)$weight
    
    plot_Networks[[Net_Net]] <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
      geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
      geom_nodes(aes(size=FreqTRAIT, fill=Mean_Prob),color="black" ,shape=21)+
      #geom_label(aes(label=Names),label.size = 0.1)+
      scale_color_viridis(direction = -1,option = "D")+
      scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background=element_blank())
    
    
    edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
    edge_dens <- igraph::edge_density(g)
    mean_grStre <- mean(igraph::graph.strength(g))
    sd_grStre <- sd(igraph::graph.strength(g))
    medi_grStre <- median(igraph::graph.strength(g))
    min_grStre <- min(igraph::graph.strength(g))
    max_grStre <- max(igraph::graph.strength(g))
    
    out_temp <- data.frame("EdG"=edge_dens_man,
                           "edge_dens_man"=edge_dens,
                           "mean_grStre"=mean_grStre,
                           "sd_grStre"=sd_grStre,
                           "medi_grStre"=medi_grStre,
                           "min_grStre"=min_grStre,
                           "max_grStre"=max_grStre,
                           "Network"="ALBERA",
                           "AddedPonds"=length(Selecte_ponds),"Iter"=iterati)
    
    out_table <- bind_rows(out_table,out_temp)
    
  }
}

out_table_ALBERA <- out_table

pond_name <- c("CLOTS15","CLOTS17","CLOTS18","CLOTS21","CLOTS22","CLOTS23","CLOTS4","CLOTS5","CLOTS6","CLOTSX3")
number_additions <- 1:10
plot_Networks <- list()
out_table <- data.frame()

for (Net_Net in 1:length(number_additions)) {
  # Falta el loop per fer diverses vegades
  for (iterati in 1:10) {
    Selecte_ponds <- sample(pond_name,number_additions[Net_Net],replace = F) 
    
    setwd("C:/Users/David CM/Desktop")
    library(cooccur);library(tidyverse)
    #LOAD 2Coocurrence - species.csv
    #load dataset
    load("abun_macro_sp.RData")
    traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
    traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)
    
    # Check presence
    indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
        TRUE ~ "altres" )) %>% 
      filter(Network=="GUILS", codi_lloc%in%Selecte_ponds)
    
    traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
      left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
      group_by(genus) %>% filter(n()<2) %>% 
      select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")
    
    #run analysis
    cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
    
    adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
    for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
      link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
      spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
      adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
    }
    adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
    
    g <- igraph::graph_from_adjacency_matrix(adj_table,weighted = T)
    
    # Layout plot (distància entre nodes=weight dels nodes)
    minC <- rep(-Inf, length(igraph::V(g)))
    maxC <- rep(Inf, length(igraph::V(g)))
    minC[1] <- maxC[1] <- 0
    
    l <- igraph::layout_with_fr(g,grid = "auto",  start.temp = 30,
                                minx = minC, maxx = maxC,
                                miny = minC, maxy = maxC,niter = 12)
    colnames(l) <- c("x","y")
    library(network);library(ggnetwork);library(viridis)
    n<- network(as.matrix(igraph::as_adjacency_matrix(g)), directed=F, diag=T)
    n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
    n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
    n %v% "Names" <-colnames(traits_ind)
    n %e% "Edg_values" <-igraph::E(g)$weight
    
    plot_Networks[[Net_Net]] <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
      geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
      geom_nodes(aes(size=FreqTRAIT, fill=Mean_Prob),color="black" ,shape=21)+
      #geom_label(aes(label=Names),label.size = 0.1)+
      scale_color_viridis(direction = -1,option = "D")+
      scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background=element_blank())
    
    
    edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
    edge_dens <- igraph::edge_density(g)
    mean_grStre <- mean(igraph::graph.strength(g))
    sd_grStre <- sd(igraph::graph.strength(g))
    medi_grStre <- median(igraph::graph.strength(g))
    min_grStre <- min(igraph::graph.strength(g))
    max_grStre <- max(igraph::graph.strength(g))
    
    out_temp <- data.frame("EdG"=edge_dens_man,
                           "edge_dens_man"=edge_dens,
                           "mean_grStre"=mean_grStre,
                           "sd_grStre"=sd_grStre,
                           "medi_grStre"=medi_grStre,
                           "min_grStre"=min_grStre,
                           "max_grStre"=max_grStre,
                           "Network"="GUILS",
                           "AddedPonds"=length(Selecte_ponds),"Iter"=iterati)
    
    out_table <- bind_rows(out_table,out_temp)
    
  }
}  
out_table_GUILS <- out_table

pond_name <- c("B1","B2","B3","E1","E3","M4","M5","M7","S4","S5","S6")  
number_additions <- 1:10
plot_Networks <- list()
out_table <- data.frame()

for (Net_Net in 1:length(number_additions)) {
  # Falta el loop per fer diverses vegades
  for (iterati in 1:10) {
    Selecte_ponds <- sample(pond_name,number_additions[Net_Net],replace = F) 
    
    setwd("C:/Users/David CM/Desktop")
    library(cooccur);library(tidyverse)
    #LOAD 2Coocurrence - species.csv
    #load dataset
    load("abun_macro_sp.RData")
    traits <-read.csv("tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
    traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)
    
    # Check presence
    indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
        TRUE ~ "altres" )) %>% 
      filter(Network=="GIARA", codi_lloc%in%Selecte_ponds)
    
    traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
      left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
      group_by(genus) %>% filter(n()<2) %>% 
      select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")
    
    #run analysis
    cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)
    
    adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
    for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
      link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
      spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
      adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
    }
    adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
    
    g <- igraph::graph_from_adjacency_matrix(adj_table,weighted = T)
    
    # Layout plot (distància entre nodes=weight dels nodes)
    minC <- rep(-Inf, length(igraph::V(g)))
    maxC <- rep(Inf, length(igraph::V(g)))
    minC[1] <- maxC[1] <- 0
    
    l <- igraph::layout_with_fr(g,grid = "auto",  start.temp = 30,
                                minx = minC, maxx = maxC,
                                miny = minC, maxy = maxC,niter = 12)
    colnames(l) <- c("x","y")
    library(network);library(ggnetwork);library(viridis)
    n<- network(as.matrix(igraph::as_adjacency_matrix(g)), directed=F, diag=T)
    n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
    n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
    n %v% "Names" <-colnames(traits_ind)
    n %e% "Edg_values" <-igraph::E(g)$weight
    
    plot_Networks[[Net_Net]] <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
      geom_edges(aes(color =Edg_values),arrow=arrow(angle = 20),curvature = 0.15, size=0.1) +
      geom_nodes(aes(size=FreqTRAIT, fill=Mean_Prob),color="black" ,shape=21)+
      #geom_label(aes(label=Names),label.size = 0.1)+
      scale_color_viridis(direction = -1,option = "D")+
      scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
      theme_classic()+
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background=element_blank())
    
    
    edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
    edge_dens <- igraph::edge_density(g)
    mean_grStre <- mean(igraph::graph.strength(g))
    sd_grStre <- sd(igraph::graph.strength(g))
    medi_grStre <- median(igraph::graph.strength(g))
    min_grStre <- min(igraph::graph.strength(g))
    max_grStre <- max(igraph::graph.strength(g))
    
    out_temp <- data.frame("EdG"=edge_dens_man,
                           "edge_dens_man"=edge_dens,
                           "mean_grStre"=mean_grStre,
                           "sd_grStre"=sd_grStre,
                           "medi_grStre"=medi_grStre,
                           "min_grStre"=min_grStre,
                           "max_grStre"=max_grStre,
                           "Network"="GIARA",
                           "AddedPonds"=length(Selecte_ponds),"Iter"=iterati)
    
    out_table <- bind_rows(out_table,out_temp)
    
  }
} 
out_table_GIARA <- out_table




setwd("C:/Users/David CM/Desktop")
library(cooccur);library(tidyverse)
#LOAD 2Coocurrence - species.csv
#load dataset
load("data/Retromed/abun_macro_sp.RData")
traits <-read.csv("data/tachet.traits.def.csv", header=TRUE, sep=";", na.strings="")
traits[,11:ncol(traits)] <- ifelse(traits[,11:ncol(traits)]>=1,1,0)

# Check presence
indv.l <- indv.l %>% filter(genus%in%traits$Genus..if.description.at.this.level.) %>% 
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
    TRUE ~ "altres" )) %>% 
  filter(Network%in%c("GIARA","ALBERA","GUILS"))

traits_ind <- indv.l %>% ungroup() %>% group_by(genus) %>% summarise(tot_ab=sum(indv.l)) %>% 
  left_join(traits, by=c("genus"="Genus..if.description.at.this.level."),multiple ="all") %>% 
  group_by(genus) %>% filter(n()<2) %>% 
  select(c(1,12:ncol(.))) %>% na.omit() %>% tibble::column_to_rownames("genus")

#run analysis
cooccur.species_ccr <- cooccur(mat = t(traits_ind), type = "spp_site", thresh = T, spp_names = TRUE)

adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)
for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {
  link_weight <- cooccur.species_ccr$results %>% group_by(sp1,sp2) %>% summarise(mean_prob=mean(prob_cooccur))
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
}
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]

g <- igraph::graph.adjacency(adj_table,mode="undirected",weighted = TRUE)
igraph::diameter(g)

# Layout plot (distància entre nodes=weight dels nodes)
minC <- rep(-Inf, length(igraph::V(g)))
maxC <- rep(Inf, length(igraph::V(g)))
minC[1] <- maxC[1] <- 0

l <- igraph::layout_with_fr(g,grid = "auto",  start.temp = 2500)
colnames(l) <- c("x","y")
library(network);library(ggnetwork);library(viridis)
n<- network(as.matrix(igraph::as_adjacency_matrix(g)), directed=F, diag=T)
n %v% "Mean_Prob" <- igraph::eigen_centrality(g)$vector
n %v% "FreqTRAIT" <-apply(ifelse(traits_ind>0,1,0),2,sum)
n %v% "Names" <-colnames(traits_ind)
n %e% "Edg_values" <-igraph::E(g)$weight

plot_Networks[[Net_Net]] <- ggplot(n, layout=as.matrix(l),aes(x = x, y = y, xend = xend, yend = yend))+
  geom_edges(aes(color =Edg_values, size=Edg_values),arrow=arrow(angle = 20),curvature = 0.15) +
  geom_nodes(aes(size=FreqTRAIT, fill=Mean_Prob),color="black" ,shape=21)+
  scale_color_viridis(direction = -1,option = "D")+
  scale_fill_viridis(direction = -1,option = "A")+labs(title="")+
  theme_classic()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background=element_blank())

edge_dens_man <- length(igraph::E(g))/(((ncol(traits_ind)*(ncol(traits_ind)-1)))/2)
edge_dens <- igraph::edge_density(g)
mean_grStre <- mean(igraph::graph.strength(g))
sd_grStre <- sd(igraph::graph.strength(g))
medi_grStre <- median(igraph::graph.strength(g))
min_grStre <- min(igraph::graph.strength(g))
max_grStre <- max(igraph::graph.strength(g))

out_temp <- data.frame("EdG"=edge_dens_man,
                       "edge_dens_man"=edge_dens,
                       "mean_grStre"=mean_grStre,
                       "sd_grStre"=sd_grStre,
                       "medi_grStre"=medi_grStre,
                       "min_grStre"=min_grStre,
                       "max_grStre"=max_grStre,
                       "Network"="MetaNAT",
                       "AddedPonds"=length(Selecte_ponds),"Iter"=iterati)

out_table <- bind_rows(out_table,out_temp)




out_table_MERGE <- bind_rows(out_table_GIARA,out_table_GUILS,out_table_ALBERA)


gridExtra::grid.arrange(

out_table_MERGE %>% ggplot(aes(y=EdG,x=AddedPonds, colour=Network))+
  geom_point()+ geom_smooth(method = "loess", se=F)+
  geom_hline(yintercept = out_temp$EdG, colour="red")+
  theme_classic()+
  labs(y="Manual EdgeDensity",title= "Man EdgDens"),

out_table_MERGE %>% ggplot(aes(y=edge_dens_man,x=AddedPonds, colour=Network))+
  geom_point()+ geom_smooth(method = "loess", se=F)+
  geom_hline(yintercept = out_temp$edge_dens_man, colour="red")+
  theme_classic()+
  labs(y="Igraph EdgeDensity",title= "Igr EdgDens"),

out_table_MERGE %>% ggplot(aes(y=mean_grStre,x=AddedPonds, colour=Network))+
  geom_point()+ geom_smooth(method = "loess", se=F)+
  geom_hline(yintercept = out_temp$mean_grStre, colour="red")+
  theme_classic()+
  labs(y="Mean graph srtength",title= "Mean GrStre"),

out_table_MERGE %>% ggplot(aes(y=sd_grStre,x=AddedPonds, colour=Network))+
  geom_point()+ geom_smooth(method = "loess", se=F)+
  geom_hline(yintercept = out_temp$sd_grStre, colour="red")+
  theme_classic()+
  labs(y="sd graph srtength",title= "SD GrStre"),

out_table_MERGE %>% ggplot(aes(y=medi_grStre,x=AddedPonds, colour=Network))+
  geom_point()+ geom_smooth(method = "loess", se=F)+
  geom_hline(yintercept = out_temp$medi_grStre, colour="red")+
  theme_classic()+
  labs(y="Median graph srtength",title= "Med GrStre"),

out_table_MERGE %>% ggplot(aes(y=min_grStre,x=AddedPonds, colour=Network))+
  geom_point()+ geom_smooth(method = "loess", se=F)+
  geom_hline(yintercept = out_temp$min_grStre, colour="red")+
  theme_classic()+
  labs(y="Min graph srtength",title= "Min GrStre"),

out_table_MERGE %>% ggplot(aes(y=max_grStre,x=AddedPonds, colour=Network))+
  geom_point()+ geom_smooth(method = "loess", se=F)+
  geom_hline(yintercept = out_temp$max_grStre, colour="red")+
  theme_classic()+
  labs(y="Max graph srtength",title= "Max GrStre"),

nrow=4)




plot_Networks[[2]]













