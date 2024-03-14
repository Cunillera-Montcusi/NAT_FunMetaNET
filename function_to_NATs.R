
fun_to_NATs <- function(Spp_x_Traits_Matrix){
require(igraph);require(cooccur)
#run analysis of coocurrence with the transposed matrix of Spp_x_Traits_Matrix
cooccur.species_ccr <- cooccur(mat = t(Spp_x_Traits_Matrix), 
                               type = "spp_site", thresh = T, spp_names = TRUE)

# Generate an adjacency table from the coocurrence matrix 
## Create a blank table with only 0 
adj_table <- matrix(nrow =cooccur.species_ccr$species,ncol = cooccur.species_ccr$species,data = 0)

for (row_ID in 1:length(unique(cooccur.species_ccr$results$sp1))) {# do the following for each spp
#cat("We are at trait interaction", row_ID, "of", length(unique(cooccur.species_ccr$results$sp1)), "\n")
  # Summarise the strenght of the link for each genus by calulating the mean probability
  link_weight <- cooccur.species_ccr$results %>% 
                 group_by(sp1) %>% 
                 summarise(mean_prob=mean(prob_cooccur))
  # filter the results from the coocurrence table for a species. You isolate all the species with who the 
  ## targeted species is coocuring
  spp_to <- cooccur.species_ccr$results %>% filter(sp1==unique(cooccur.species_ccr$results$sp1)[row_ID])
  # You add the values into the "adj.table" that was empty and now is being filled
  adj_table[unique(spp_to$sp1),spp_to$sp2] <- spp_to$prob_cooccur
}
# We transpose the results to make the matrix symetrical (same from & to links)
adj_table[lower.tri(adj_table)] <- t(adj_table)[lower.tri(adj_table)]
# We built a graph from the adjacency table
g <- igraph::graph.adjacency(adj_table,mode="undirected",weighted = TRUE)

output <- list("Adj_Table"=adj_table,"Graph"=g)
output
}
