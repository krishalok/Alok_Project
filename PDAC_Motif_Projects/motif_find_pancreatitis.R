setwd("C:\\Users\\alokum\\Desktop\\Menuscript\\PNAS\\Papers\\Main_Motif_nature\\GA_model")
setwd("C:\\Users\\alokum\\Desktop\\Menuscript\\PNAS\\Menuscript\\data")
setwd("C:\\Users\\alokum\\Desktop\\PNAS\\Menuscript\\data")
rm(list=ls())
library(igraph)
library(iplots)
library(linkcomm)
library(ggplot2)
library(zoo)
library(lattice)
library(NetIndices)
#proteinList <- read.table("PAAD_GE_Co_Network.csv", header=TRUE, sep=",")
proteinList <- read.csv("pancreatitis.csv", header=TRUE, sep=",")
#proteinList <- read.table("ADEX_out_0.05.csv", header=TRUE, sep=",")
#proteinList <- read.table("Immunogenic_out_0.05.csv", header=TRUE, sep=",")
#proteinList <- read.table("Pancreatic_Progenitor_out_0.05.csv", header=TRUE, sep=",")
#proteinList <- read.table("Squamous_output_0.05.csv", header=TRUE, sep=",")
t <- as.data.frame(proteinList)
#############Calculate Z-scare#############
#################https://stackoverflow.com/questions/6148050/creating-z-scores
#t$Correlation<-scale(t$Correlation)
t <- t[,1:3]
y <- graph.data.frame(t, directed=FALSE, vertices=(NULL))

plot.igraph(graph.adjlist(t$InteractorA))
########################print motifs#########################
pattern <- graph.full(3)
#my.graph <- grg.game(100, 0.2)        # just an example graph, use yours
iso <- subgraph_isomorphisms(pattern, y)   # takes a while
motifs <- lapply(iso, function (x) { induced_subgraph(y, x) })
plot.igraph(motifs)
plot(motifs)
head(motifs)
capture.output(motifs, file = "normal_final_motifs.txt")
class(motifs)
## transform data
require(reshape2)
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

#####################just try
library(data.tree)
str(motifs,max.level=1)

# Function to find network motifs in a directed graph
findnetworkmotifs <- function(graph, motifsize)
{
  library(igraph)
  # Find network motifs in the graph "graph":
  mygraphmotifs <- graph.motifs(graph, motifsize)
  # Find which motifs occur:
  nummotifs <- length(mygraphmotifs)
  for (i in 1:nummotifs)
  {
    motif <- mygraphmotifs[i]
    if (motif > 0) # There are some occurrences of this motif
    {
      # Find out what the motif looks like:
      motifgraph <- graph.isocreate(size=motifsize, number=i-1, directed=TRUE)
      edges <- E(motifgraph)
      print(paste("This motif occurs",motif,"times:"))
      print(edges)
    }
  }
}




subGraph = graph.neighborhood(y, order = 1, V(y), mode = 'all')[[1]]
allMotifs = triad.census(subGraph)
removeNode = delete.vertices(subGraph, 'one')
node1Motifs = allMotifs - triad.census(removeNode)

#####count of motif

triads18 <- triad.census(y)
motifs18 <- graph.motifs(y,3)

degree(y)
plot(dendrogram, y)
dd <- degree(y)
max(dd)
dd %>% table() %>% barplot()
degree_distribution(y)

y %>% summary()

ggplot(motifs, aes(x = Seeds, y = Correct)) +
geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
geom_line() +
geom_point() +
labs(x = "Seeds", y = "Fraction correctly matched",
title = "Effect of seeding vs correlation strength")

motifs(y, 3)
count_motifs(y, 3)
sample_motifs(y, 4)
dyad_census(y)
triad_census(y)
########################################################################################################################
########################################################################################################################


#Extract the linking communities
lc <- getLinkCommunities(t, hcmethod = "average", use.all.edges=TRUE, edglim = 10^4, directed=TRUE, plot = TRUE)
#Find overlapping communities

oc <-get.community.overlaps(lc)

#Get nodes from a community
getNodesIn(lc, clusterids = 4)

integer.edgelist(t)

linkcomm2cytoscape(lc, interaction = "pp", ea = "temp.ea")

######################### couting motifs#############################
motifs(y, 3, cut.prob = rep(0, 3))
#g <- barabasi.game(100)
motifs(y,3)
count_motifs(y, 3)
sample_motifs(y, 3)
degree(y)
degree_distribution(y)

plot(y)
plot(y, edge.arrow.size=.4,vertex.label=NA)


# Plot the egree distribution for our network:
deg.dist <- degree_distribution(y, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")


############print vertices in motifs#########################





# Set layout options
#l <- layout.fruchterman.reingold(y)

# Plot graph and subgraph
#plot.igraph(x=motifs,layout=l)


library(igraph)

net <- graph_from_data_frame(d=proteinList, vertices=nodes, directed=T) 


class(y)
E(y)
V(y)
plot(y, edge.arrow.size=.4,vertex.label=NA)

edge.start <- ends(y, es=E(y), names=F)[,1]

edge.col <- V(y)$color[edge.start]

plot(y, edge.color=edge.col, edge.curved=.1)
ceb <- cluster_edge_betweenness(y) 


dendPlot(ceb, mode="hclust")


deg <- degree(y, mode="all")
plot(ceb, y) 
deg.dist <- degree_distribution(y, cumulative=T, mode="all")

plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      
      xlab="Degree", ylab="Cumulative Frequency")

plot(y, vertex.size=deg*3)

net.sym <- as.undirected(y, mode= "collapse",
                         
                         edge.attr.comb=list(weight="sum", "ignore"))
cliques(net.sym) # list of cliques       

sapply(cliques(net.sym), length) # clique sizes

largest_cliques(net.sym) # cliques with max number of nodes


transitivity(y, type="global")  # net is treated as an undirected network

transitivity(as.undirected(y, mode="collapse")) # same as above

transitivity(y, type="local")

triad_census(y) # for directed networks
cocitation(y)
####################################
net.sym <- as.undirected(y, mode= "collapse",
                         
                         edge.attr.comb=list(weight="sum", "ignore"))
cliques(net.sym) # list of cliques       

sapply(cliques(net.sym), length) # clique sizes

largest_cliques(net.sym) # cliques with max number of nodes

vcol <- rep("grey80", vcount(net.sym))

vcol[unlist(largest_cliques(net.sym))] <- "gold"

plot(as.undirected(net.sym), vertex.label=V(net.sym)$name, vertex.color=vcol)    


##############################################

source("https://bioconductor.org/biocLite.R")
biocLite("pandaR")
library(pandaR)
data(pandaToyData)
pandaResult <- panda(pandaToyData$motif, pandaToyData$expression, pandaToyData$ppi)
pandaResult <- panda(proteinList$InteractorA,proteinList$InteractorB,proteinList$Correlation)
pandaResult
topNet <- topedges(pandaResult, 1000)
topSubnet <- subnetwork(topNet, c("AR","ARID3A","ELK1"))
plotGraph(topSubnet)
panda.res1 <- with(pandaToyData, panda(motif, expression, ppi, hamming=1))
panda.res2 <- with(pandaToyData, panda(motif, expression + rnorm(prod(dim(expression)),sd=5), ppi, hamming=1))
plotZ(panda.res1, panda.res2,addLine=FALSE)



##############################################
bridge(y, communities = NULL, useCommunities = "all",directed = NULL, nodes = NULL, average = FALSE)
b2 <- bridge(y, communities="all")

triad_census(y)
scan_stat(y, k = 1, tau = 4, ell = 2)
co <- components(y, mode = "weak")
groups(co)[[2]]
dendrogram <- cluster_edge_betweenness(y)
dendrogram
compare_all <- function(cl1, cl2) {
  methods <- eval(as.list(args(compare))$method)
  vapply(methods, compare, 1.0, comm1 = cl1, comm2 = cl2)
}


par(mar=c(0,0,0,0)); plot_dendrogram(dendrogram, direction = "downwards")






