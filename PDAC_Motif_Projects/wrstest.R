#Test the significance of the communities using Wilcoxon Rank-Sum test
community.significance.test <- function(g, vs) {
  if (is.directed(g)) stop("This method requires an undirected graph")
  subgraph <- induced.subgraph(g, vs)
  in.degrees <- degree(subgraph)
  out.degrees <- degree(g, vs) - in.degrees
  result <- wilcox.test(in.degrees, out.degrees) 
  return(result)
}