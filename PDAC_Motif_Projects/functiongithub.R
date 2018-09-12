library(igraph)
genes <- read.table("~/Desktop/hgnc_mapped.txt", header=TRUE, sep="\t")
g <- as.data.frame(genes)
y <- graph.data.frame(g, directed=FALSE, vertices=NULL)
y <- simplify(y)
z <- graph.data.frame(e, directed=FALSE, vertices=NULL)

################################################################
## Combine undirected interaction networks into a             ##
## single graph. Attributes will be merged when possible.     ##
################################################################
mergeGraphs = function(graphs, setSourceAttr = T, firstAsBase = F) {
  interactions = graph.empty()
  if(firstAsBase) {
    interactions = graphs[[1]]
    graphs = graphs[2:length(graphs)]
    vas = list.vertex.attributes(interactions)
    if(!("name" %in% vas) && "identifier" %in% vas) V(interactions)$name = V(interactions)$identifier
  }
  
  for(graphName in names(graphs)) {
    message("processing ", graphName)
    g = graphs[[graphName]]
    vas = list.vertex.attributes(g)
    if(!("name" %in% vas) && "identifier" %in% vas) {
      print("Setting identifier attribute as name")
      V(g)$name = V(g)$identifier
    }
    
    #Add nodes
    newv = V(g)[!(name %in% V(interactions)$name)]
    newnames = V(g)[newv]$name
    attr = list()
    attr$name = newnames
    
    interactions = add.vertices(
      interactions, length(newv), attr = attr
    )
    
    ## Merge node attributes
    newAttr = setdiff(list.vertex.attributes(g), list.vertex.attributes(interactions))
    for(an in newAttr) {
      interactions = set.vertex.attribute(interactions, an, V(interactions), "")
      interactions = set.vertex.attribute(interactions, an, V(interactions)[V(g)$name], get.vertex.attribute(g, an))
    }
    mergeAttr = intersect(list.vertex.attributes(g), list.vertex.attributes(interactions))
    for(an in mergeAttr) {
      if(length(newnames) > 0)
        interactions = set.vertex.attribute(interactions, an, V(interactions)[newnames], get.vertex.attribute(g, an, V(g)[newnames]))
      
      ovl = setdiff(V(g)$name, newnames)
      if(length(ovl) > 0) {
        ovl.attr = cbind(
          get.vertex.attribute(interactions, an, V(interactions)[ovl]),
          get.vertex.attribute(g, an, V(g)[ovl])
        )
        ovl.attr = apply(ovl.attr, 1, function(x) {
          xs = gsub("([\\.\\(\\)\\+]{1})", "\\\\\\1", x[2])
          
          if(length(x[[1]]) == 0 || is.na(x[[1]])) x[[2]]
          else if(length(x[[2]]) == 0 || is.na(x[[2]])) x[[1]]
          else if(
            x[[1]] != x[[2]] && 
              length(grep(paste(xs, ", ", sep=""), x[[1]])) == 0 &&
              length(grep(paste(xs, "$", sep=""), x[[1]])) == 0) {
            paste(x, collapse = ", ")
          } else {
            x[[1]]
          }
        })
        interactions = set.vertex.attribute(interactions, an, V(interactions)[ovl], ovl.attr)
      }
    }
    
    #Add edges
    newNodeIndex = as.numeric(V(interactions)[V(g)$name])
    edges = numeric(ecount(g)*2)
    el = get.edgelist(g, names = F)
    for(n in 1:nrow(el)) {
      if(n %% 1000 == 0) message(n, " out of ", ecount(g))
      edges[c(n*2-1,n*2)] = c(newNodeIndex[el[n,1]], newNodeIndex[el[n,2]])
    }
    attr = list()
    for(n in list.edge.attributes(g)) {
      attr[[n]] = as.character(get.edge.attribute(g, n))
    }
    if(setSourceAttr) attr$SourceFile = rep(graphName, ecount(g))
    attr$sourceEdge = as.numeric(E(g))
    interactions = add.edges(interactions, edges, attr = attr)
  }
  interactions
}