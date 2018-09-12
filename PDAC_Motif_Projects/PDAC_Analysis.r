
#Install the packages
install.packages("devtools")
devtools::install_github("hadley/multidplyr")
require(multidplyr)
library(parallel)
library(tidyverse)
library(igraph)
library(NetIndices)
cl <- detectCores()
cluster <- create_cluster(cores = cl)
cluster %>%
  # Assign libraries
  cluster_library("igraph") %>%
  cluster_library("tidyverse") %>%
  cluster_library("magrittr") %>%
  cluster_library("dplyr") %>%
  cluster_library("RColorBrewer") %>%
  cluster_library("NetIndices") %>%
  # Assign values (use this to load functions or data to each core)
  cluster_assign_value("GenInd", 6)





setwd("C:\\Users\\alokum\\Desktop\\PNAS\\Menuscript\\data")
library(data.table)
library(igraph)
library(NetIndices)
library(reshape2)
url <- "./PAAD_GE_Co_Network.csv"
proteinList <- fread(url)
t <- as.data.frame(proteinList)
t <- t[,1:3]
y <- graph.data.frame(t, directed=FALSE, vertices=(NULL))
pattern <- graph.full(3)
#my.graph <- grg.game(100, 0.2)        # just an example graph, use yours
iso <- subgraph_isomorphisms(pattern, y)   # takes a while
motifs <- lapply(iso, function (x) { induced_subgraph(y, x) })
capture.output(motifs, file = "PDAC_final_motifs.txt")
h <- do.call(cbind, motifs)
gc()
h.melt <- melt(h)
ggplot(h.melt,aes(x=value,fill=N))+geom_density()+facet_grid(N~.)+geom_rug()
## transform dataNA     NA 353174(total)   3180 (Trinagle)
graph.motifs(y, size = 3)
triad.census(y)
triangle <- graph.full(3)
graph.count.subisomorphisms.vf2(y, triangle) / 6
graph.motifs.no(y,size=3)
plot(y, layout = layout.reingold.tilford(y, root=3))
motifs(y, size = 3, cut.prob = rep(0, 3))
count_motifs(y, 3)
sample_motifs(y, 3)
edges <- E(y)
nodes <- V(y)
plot.igraph(y,
            layout=layout.lgl,
            vertex.size=2,
            vertex.label.cex=.5,
            edge.arrow.size=.5)
## transform data from https://assemblingnetwork.wordpress.com/2013/06/10/network-basics-with-r-and-igraph-part-ii-of-iii/
test.graph.adj<-get.adjacency(y,sparse=F)
result <- clusterMap(cluster, function1)

test.graph.properties<-GenInd(test.graph.adj)
test.graph.properties$N            #number of nodes

test.graph.properties$Ltot        #number of links

test.graph.properties$LD        #link density (average # of links per node)

test.graph.properties$C            #the connectance of the graph

in.deg.testgraph<-degree(y,v=V(y),mode="in")
out.deg.testgraph<-degree(y,v=V(y),mode="out")
all.deg.testgraph<-degree(y,v=V(y),mode="all")

# Degree distribution is the cumulative frequency of nodes with a given degree
# this, like degree() can be specified as "in", "out", or "all"
deg.distr<-degree.distribution(y,cumulative=T,mode="all")

# Using the power.law.fit() function I can fit a power law to the degree distribution
power<-power.law.fit(all.deg.testgraph)

# Then I can plot the degree distribution
plot(deg.distr,log="xy",
     ylim=c(0.01,10),
     bg="black",pch=20,
     xlab="Degree",
     ylab="Cumulative Frequency")

# And the expected power law distribution
lines(1:20,10*(1:20)^((-power$alpha)+1))


# Diameter is essentially the longest path between two vertices
diameter(y)
# Gives me the length of the diameter while

nodes.diameter<-get.diameter(y)
# Gives me the labels for each node that participates in the diameter

# I can look at the diameter graphically also
# First I will define the node and edge attributes
V(y)$color<-"skyblue"
# I want all the nodes to be skyblue
V(y)$size<-7
# I want all the nodes to be size=7
V(y)[nodes.diameter]$color<-"darkgreen"
V(y)[nodes.diameter]$size<-10
V(y)[nodes.diameter]$label.color<-"white"
# but the nodes in the diameter should be darkgreen and larger than the rest
# with a white label instead of black
# this will make the diameter pop out of the larger network
E(y)$color<-"grey"
# all non-diameter edges will be grey
E(y,path=nodes.diameter)$color<-"darkgreen"
E(y,path=nodes.diameter)$width<-2
# Edges in the diameter will be darkgreen and a little extra wide

# If you do not set the attributes of all of the nodes and edges then it will
# default such that you only see what you have defined

# Now when I plot the diameter will be larger than everything else, and darkgreen instead
# of grey/blue
par(mar=c(.1,.1,.1,.1))
plot.igraph(y,
            layout=layout.fruchterman.reingold,
            vertex.label.cex=.5,
            edge.arrow.size=.5)

# Clustering coefficient is the proportion of
# a nodes neighbors that can be reached by other neighbors
# in igraph this property is apparently called "transitivity"

transitivity(y)
# gives the clustering coefficient of the whole network

transitivity(y,type="local")
# gives the clustering coefficient of each node

# Betweenness is the number of shortest paths between two nodes that go through each node of interest

graph.betweenness<-betweenness(y,v=V(y))
graph.edge.betweenness<-edge.betweenness(y,e=E(y))

# Closeness refers to how connected a node is to its neighbors

graph.closeness<-closeness(y,vids=V(y))

# Clustering coefficient, betweenness, and closeness
# all describe the small world properties of the network.
# A network with small world properties is one in which
# it takes a relatively short path to get from one node to the next
# (e.g., six degrees of separation)

# Every graph can be decomposed into its component n-node subgraphs.
# In particular there are 13 unique ways to arrange 3 nodes in directed graphs.
# Here are the adjacency matrices for each of the 13 subgraphs
s1<-matrix(c(0,1,0,0,0,1,0,0,0),nrow=3,ncol=3)
s2<-matrix(c(0,1,1,0,0,1,0,0,0),nrow=3,ncol=3)
s3<-matrix(c(0,1,0,0,0,1,1,0,0),nrow=3,ncol=3)
s4<-matrix(c(0,0,1,0,0,1,0,0,0),nrow=3,ncol=3)
s5<-matrix(c(0,1,1,0,0,0,0,0,0),nrow=3,ncol=3)
d2<-matrix(c(0,1,1,1,0,1,0,0,0),nrow=3,ncol=3)
d1<-matrix(c(0,1,1,0,0,1,0,1,0),nrow=3,ncol=3)
d3<-matrix(c(0,0,1,1,0,0,1,0,0),nrow=3,ncol=3)
d4<-matrix(c(0,0,0,1,0,1,0,1,0),nrow=3,ncol=3)
d5<-matrix(c(0,1,1,0,0,1,1,0,0),nrow=3,ncol=3)
d6<-matrix(c(0,1,1,1,0,1,1,1,0),nrow=3,ncol=3)
d7<-matrix(c(0,1,1,1,0,1,1,0,0),nrow=3,ncol=3)
d8<-matrix(c(0,1,1,1,0,0,1,0,0),nrow=3,ncol=3)

# I then make the 13 matrices into a list
subgraph3.mat<-list(s1,s2,s3,s4,s5,d1,d2,d3,d4,d5,d6,d7,d8)
# And convert the matrices into graph objects
subgraph3.graph<-lapply(subgraph3.mat,graph.adjacency)

# Here I have created a simple for loop to go through the list of subgraphs
# and count how many times that subgraph appears in the larger test.graph
subgraph.count<-c()
for(i in 1:13){
  subgraph.count[i]<-
    graph.count.subisomorphisms.vf2(y,subgraph3.graph[[i]])
}

plot(y,type="o",lty=3, xlab="Subgraph",ylab="Frequency")
points(y,type="o",lty=2)

plot(subgraph.freq.otago~subgraph.freq.otago.p)
abline(a=0,b=1)