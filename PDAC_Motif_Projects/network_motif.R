#Function to find network motifs in a directed graph
library(igraph)
findmotifs <- function(graph, motifsize){
  my_graph <- graph.motifs(graph, motifsize)
  number_motifs <- length(my_graph)
  for (i in 1:number_motifs){
    motif <- my_graph[i]
    if (motif>0) {
      motif_graph <- graph.isocreate(size=motifsize, number=i-1, directed=TRUE)
      edges <- E(motif_graph)
      print(paste0("This motif occurs", motif, "times:"))
      print(edges)
    }
  }
  
  
  
  
  
  
  
  
}
