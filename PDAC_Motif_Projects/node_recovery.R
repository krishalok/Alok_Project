library(GraphAlignment)
library(igraph)
setwd("/Users/BC/Documents/spinglass/results/matching_communities/")


temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.delim)

run1 <- read.csv("community1/run1-comm4.csv", header=T)
run2 <- read.csv("community1/run2-comm24.csv", header=T)
run3 <- read.csv("community1/run3-comm14.csv", header=T)
run4 <- read.csv("community1/gamma-2-start1-stop-.8-comm-12.csv", header=T)
run5 <- read.csv("community1/run5-comm-3.csv", header=T)
run6 <- read.csv("community1/run6-comm-8.csv", header=T)
run7 <- read.csv("community1/run7-comm-19.csv", header=T)
run8 <- read.csv("community1/run8-comm-14.csv", header=T)
run9 <- read.csv("community1/run9-comm-15.csv", header=T)
run10 <- read.csv("community1/run10-comm-13.csv", header=T)


y1=as.matrix(run1[,2:3])
y2=as.matrix(run2[,2:3])
y3=as.matrix(run3[,2:3])
y4=as.matrix(run4[,2:3])
y5=as.matrix(run5[,2:3])
y6=as.matrix(run6[,2:3])
y7=as.matrix(run7[,2:3])
y8=as.matrix(run8[,2:3])
y9=as.matrix(run9[,2:3])
y10=as.matrix(run10[,2:3])

g1 <- graph.edgelist(y1, directed=FALSE)
g2 <- graph.edgelist(y2, directed=FALSE)
g3 <- graph.edgelist(y3, directed=FALSE)
g4 <- graph.edgelist(y4, directed=FALSE)
g5 <- graph.edgelist(y5, directed=FALSE)
g6 <- graph.edgelist(y6, directed=FALSE)
g7 <- graph.edgelist(y7, directed=FALSE)
g8 <- graph.edgelist(y8, directed=FALSE)
g9 <- graph.edgelist(y9, directed=FALSE)
g10 <- graph.edgelist(y10, directed=FALSE)


g_1 <- get.adjacency(g1)
g_2 <- get.adjacency(g2)
g_3 <- get.adjacency(g3)
g_4 <- get.adjacency(g4)
g_5 <- get.adjacency(g5)
g_6 <- get.adjacency(g6)
g_7 <- get.adjacency(g7)
g_8 <- get.adjacency(g8)
g_9 <- get.adjacency(g9)
g_10 <- get.adjacency(g10)

g_1.1 <- graph.adjacency(g_1)
g_2.1 <- graph.adjacency(g_2)
g_3.1 <- graph.adjacency(g_3)
g_4.1 <- graph.adjacency(g_4)
g_5.1 <- graph.adjacency(g_5)
g_6.1 <- graph.adjacency(g_6)
g_7.1 <- graph.adjacency(g_7)
g_8.1 <- graph.adjacency(g_8)
g_9.1 <- graph.adjacency(g_9)
g_10.1 <- graph.adjacency(g_10)

#Calculate the node similarity between two graphs
g_sim <- graph.intersection(g_1.1, g_2.1,g_3.1,g_4.1,g_5.1,g_6.1,g_7.1,g_8.1,g_9.1,g_10.1, byname = "auto", keep.all.vertices = FALSE)

adj_sim <- get.adjacency(g_sim, type="both")

g<- graph.adjacency(adj_sim,mode = "undirected",weighted = T)

edge_ara <- get.data.frame(g,what = "edges")

edge_ara<- edge_ara[order(abs(edge_ara$weight),decreasing = T),]

write.csv(edge_ara, "node_sim.csv")

