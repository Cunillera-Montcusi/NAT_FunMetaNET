
edge.weights <- function(community, network, weight.within = 100, weight.between = 1) {
  bridges <- igraph::crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights) 
}

G<- NATs_Output$Graph
# Detect isolated nodes
isolated <- igraph::V(G)[igraph::degree(G) == 0]
# Detect communities for connected nodes
connected <- igraph::induced_subgraph(G, igraph::V(G)[igraph::degree(G) > 0])
community <- igraph::cluster_louvain(connected)

V(g)$strength <- strength(g, weights = E(g)$weight)
V(g)$betweenness <- betweenness(g, weights = E(g)$weight)

igraph::V(connected)$module <- igraph::membership(community)
module_stats <- lapply(unique(igraph::V(connected)$module), function(mod) {
  # Get nodes in the module
  nodes_in_module <- igraph::V(connected)[module == mod]
  
  # Extract metrics
  strengths <- igraph::V(connected)[nodes_in_module]$strength
  betweennesses <- igraph::V(connected)[nodes_in_module]$betweenness
  
  # Summary statistics for the module
  list(
    module = mod,
    nodes = names(nodes_in_module),
    strength_mean = mean(strengths),
    strength_sd = sd(strengths),
    betweenness_mean = mean(betweennesses),
    betweenness_sd = sd(betweennesses)
  )
})


#  Plot the network 
# Modify weights to better visualise modules
igraph::E(connected)$weight <- igraph::E(connected)$weight*edge.weights(community, connected,10)

# Assign colors to communities
Module_Groups <- LETTERS[1:length(unique(community$membership))]
igraph::V(connected)$Module_Belonging <- Module_Groups[igraph::membership(community)]
#
# Layout for connected nodes
layout_connected <- igraph::layout_with_fr(connected)
#layout_combined <- rbind(
#  layout_connected, 
#  matrix(runif(2 * length(isolated), -1, 1), ncol = 2)  # Random layout for isolated
#)
#
#igraph::V(G)$Module_Belonging <- NA
#igraph::V(G)$Module_Belonging[igraph::degree(G) > 0] <- igraph::V(connected)$Module_Belonging
#igraph::V(G)$Module_Belonging[igraph::degree(G) == 0] <- NA  # Assign gray color to isolated nodes

ggnetwork::ggnetwork(igraph::simplify(connected),layout=layout_connected) %>%
  #filter(weight>900) %>% 
  #mutate(weight=ifelse(weight>900,weight,0.01)) %>% 
  ggplot()+
  ggnetwork::geom_edges(aes(x=x, xend=xend,y=y,yend=yend,alpha=weight),
                        curvature = 0.1)+
  ggnetwork::geom_nodes(aes(x=x,y=y,fill=Module_Belonging),shape=21,size=3)+
  viridis::scale_colour_viridis(direction = -1)+
  theme_void()+
  theme(legend.position = "none")

  