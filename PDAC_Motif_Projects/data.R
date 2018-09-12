karate <- graph.famous("Zachary")
fc <- infomap.community(karate)
E(karate)[crossing(fc, karate)]  # this would give the edges between communities 
plot.communities(fc,karate)


#Get actual vertex
library(igraph)
karate <- graph.famous("Zachary")
wckarate <- walktrap.community(karate)
V(karate)$name<-V(karate)
V(karate)$community<-membership(wckarate)

V(karate)$size<-sizes(wckarate)[membership(wckarate)]
graph_neighbours<-graph.neighborhood(karate,1,nodes=V(karate))
for(i in 1:vcount(karate)) {
  write.graph(graph_neighbours[[i]],file=c(as.character(V(karate)$name[i])),format="dot")
}

#vertex proportional to degree probability
g <- graph.tree(10, 3, mode="undirected")
sample(V(g), 1, prob=degree(g))
# Vertex sequence:
# [1] 9

#Or, if you want vertex names:
  
V(g)$name <- letters[1:10]
sample(V(g)$name, 1, prob=degree(g))
# [1] "b"


#Break up the graph into communities with 2-6 nodes
x <- which(sizes(wc) <= 6)
subg <- induced.subgraph(y, which(sizes(wc) <= 6))
#scc <- clusters(subg, "strong")
length(subg)
length(y)
#x2 <- which.max(sizes(subwc))
#subg2 <- induced.subgraph(subg, which(membership(subwc) ==x2))
#vcount(subg2)

#colors <- rainbow(max(membership(wc)))
#plot(subg,vertex.color=colors[membership(wc)], 
#     layout=layout.fruchterman.reingold)

#plot(subg2, layout=layout.kamada.kawai,vertex.size=7, vertex.label=V(y)$name, vertex.color=membership(wc), vertex.frame.color="red", edge.color="grey", edge.arrow.size=0.01, rescale=TRUE,vertex.label=TRUE, vertex.label.dist=0.0, vertex.label.cex=0.5, add=FALSE,   vertex.label.font=1)
