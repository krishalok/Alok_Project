
library(igraph)
library(NetIndices)
# Kim N. Mouritsen, Robert Poulin, John P. McLaughlin and David W. Thieltges. 2011.
# Food web including metazoan parasites for an intertidal ecosystem in New Zealand.
# Ecology 92:2006.

# Website: http://esapubs.org/archive/ecol/E092/173/

# Otago Harbour: intertidal mudflat
#otago.links.data<-read.csv("~/Desktop/Projects/FoodwebAlpha/Data/Otago_Data_Links.csv")
#otago.nodes.data<-read.csv("~/Desktop/Projects/FoodwebAlpha/Data/Otago_Data_Nodes.csv")
otago.links.data<-as_edgelist(y, names = TRUE)
otago.nodes.data<-V(y)

# Column names for data
colnames(otago.links.data)
colnames(otago.nodes.data)

# Convert the data into a graph object using the first 2 columns of the dataset as an edgelist
otago.graph<-graph.edgelist(as.matrix(otago.links.data[,1:2]))
# Create graph object of just predator prey links
#otago.graph.p<-graph.edgelist(as.matrix(otago.links.data[,1:2]))
otago.graph.p<-graph.edgelist(as.matrix(otago.links.data[1:18961,1:2]))


# Get the web into matrix form
otago.adjmatrix<-get.adjacency(otago.graph,sparse=F)
otago.adjmatrix.p<-get.adjacency(otago.graph.p,sparse=F)

# Get the basic network indices from the matrices with GenInd()
# with parellaization
ind.otago<-lapply(length(y),GenInd(otago.adjmatrix))
# without parellaization
ind.otago<-GenInd(otago.adjmatrix)
ind.otago.p<-GenInd(otago.adjmatrix.p)

ind.otago<-mclapply(length(y),GenInd(otago.adjmatrix))

# Now to plot these two webs to get a feel for what we are dealing with

par(mar=c(.1,.1,.1,.1))
plot.igraph(otago.graph,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,layout=layout.circle)
plot.igraph(otago.graph.p,vertex.label=NA,vertex.size=3,edge.arrow.size=.25,layout=layout.circle)


# The NetIndices package also has a function to get some of the trophic properties of the food web
# TrophInd() takes in an adjacency matrix and gives an output of the trophic level of each node,
# as well as an index of the degree of omnivory for each node

troph.otago<-TrophInd(otago.adjmatrix)
troph.otago.p<-TrophInd(otago.adjmatrix.p)

# An interesting aside, by adding parasites to the web it increases the trophic level of all species in
# this web.

plot(troph.otago[1:123,1]~troph.otago.p[,1],xlab="Level Without Parasites",ylab="Level With Parasites")
abline(a=0,b=1)




# An interesting use for this trophic level function is to then use trophic level as a plotting parameter.
# This way, I can plot the food web nodes according to trophic height. I think that this adds greatly to a plot
# of a food web, since you can gain more information about the trophic structure of the web by simply
# glancing at the plot.

# First we need to create a two-column matrix identifying the x and y values for each node.
layout.matrix.1<-matrix(
  nrow=length(V(otago.graph)),  # Rows equal to the number of vertices
  ncol=2
)
layout.matrix.1[,1]<-runif(length(V(otago.graph))) # randomly assign along x-axis
layout.matrix.1[,2]<-troph.otago$TL # y-axis value based on trophic level

layout.matrix.1p<-matrix(
  nrow=length(V(otago.graph.p)),  # Rows equal to the number of vertices
  ncol=2
)
layout.matrix.1p[,1]<-runif(length(V(otago.graph.p)))
layout.matrix.1p[,2]<-troph.otago.p$TL

# Now we can use these matrices to define the layout instead of using the circle layout

par(mar=c(.1,.1,.1,.1),mfrow=c(1,2))

plot.igraph(otago.graph,
            vertex.label.cex=.35,
            vertex.size=3,
            edge.arrow.size=.25,
            layout=layout.matrix.1)

plot.igraph(otago.graph.p,
            vertex.label.cex=.35,
            vertex.size=3,
            edge.arrow.size=.25,
            layout=layout.matrix.1p)

# I am still working on the best way to plot the nodes along the x-axis. You may notice that using
# runif() means that there is some chance that two nodes with the same trophic level
# will be right on top of one another


# It is also a bit interesting to see how the inclusion of parasites impacts community detection
wtc.otago<-walktrap.community(otago.graph)
wtc.otago.p<-walktrap.community(otago.graph.p)

par(mar=c(.1,.1,.1,.1),mfrow=c(1,2))

plot.igraph(otago.graph,
            vertex.label.cex=.35,
            vertex.size=3,
            edge.arrow.size=.25,
            layout=layout.matrix.1,
            mark.groups=wtc.otago$membership,
            mark.col="green")

plot.igraph(otago.graph.p,
            vertex.label.cex=.35,
            vertex.size=3,
            edge.arrow.size=.25,
            layout=layout.matrix.1p,
            mark.groups=wtc.otago.p$membership,
            mark.col="green")

# It is clear that the increase in the connectivity of the web with parasites has led to
# a larger densely connected community

# It is clear that the increase in the connectivity of the web with parasites has led to
# a larger densely connected community

# The degree distribution of a food web can tell us a lot about the amount of specialization and
# generalization in the web (in degree), as well as vulnerability (out degree)

deg.otago<-degree(otago.graph)
deg.otago.p<-degree(otago.graph.p)

# Using the degree distribution gives a better way to visualize any differences
# Looking at the in degree tells us about how general the diets of consumers are
dd.otago.in<-degree.distribution(otago.graph,mode="in",cumulative=T)
dd.otago.in.p<-degree.distribution(otago.graph.p,mode="in",cumulative=T)

# Out degree is a measure of the vulnerability of organisms, telling us how many consumers
# eat each species.
dd.otago.out<-degree.distribution(otago.graph,mode="out",cumulative=T)
dd.otago.out.p<-degree.distribution(otago.graph.p,mode="out",cumulative=T)

# And finally the degree ("all") simply tells us about how well connected that species is
# within the network
dd.otago<-degree.distribution(otago.graph,mode="all",cumulative=T)
dd.otago.p<-degree.distribution(otago.graph.p,mode="all",cumulative=T)

par(mfrow=c(2,2))
plot(dd.otago.in,xlim=c(0,80))
plot(dd.otago.out,xlim=c(0,80))
plot(dd.otago.in.p,xlim=c(0,80))
plot(dd.otago.out.p,xlim=c(0,80))

power.fit<-power.law.fit(deg.otago)
power.fit.p<-power.law.fit(deg.otago.p)

par(mfrow=c(1,2))
plot(dd.otago,log="xy")
lines(1:180,10*(1:180)^((-power.fit$alpha)+1))

plot(dd.otago.p,log="xy")
lines(1:100,10*(1:100)^((-power.fit.p$alpha)+1))

# I can look at the diameter of the two versions of the web
# For food webs the diameter is going to be the longest food chain
# since energy only flows in one direction, the diameter will read from
# basal species to top predator.

get.diameter(otago.graph)
get.diameter(otago.graph.p)

# I think that here it is interesting to note that the diameter of the predator-prey only
# food web (which we expect to be smaller) is not a subset of the diameter for the
# larger parasites included network

# The next few properties are all related to the small world-ness of the network:

transitivity(otago.graph)
transitivity(otago.graph.p)

# Betweenness is the number of shortest paths going through a specified node or edge

otago.between<-betweenness(otago.graph)
otago.between.p<-betweenness(otago.graph.p)

plot(otago.between[1:123]~otago.between.p)
abline(a=0,b=1)

otago.edge.between<-edge.betweenness(otago.graph)
otago.edge.between.p<-edge.betweenness(otago.graph.p)

closeness(otago.graph)

# Here are the adjacency matrices for each of the 13 subgraphs again
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

# Turn them into a convenient list
subgraph3.mat<-list(s1,s2,s3,s4,s5,d1,d2,d3,d4,d5,d6,d7,d8)
# And then into a list of graph objects
subgraph3.graph<-lapply(subgraph3.mat,graph.adjacency)

# Count the number of the 13 different 3-node subgraphs in the two webs
subgraph.freq.otago<-c()
subgraph.freq.otago.p<-c()
for(i in 1:13){
  subgraph.freq.otago[i]<-
    graph.count.subisomorphisms.vf2(otago.graph,subgraph3.graph[[i]])
  subgraph.freq.otago.p[i]<-
    graph.count.subisomorphisms.vf2(otago.graph.p,subgraph3.graph[[i]])
}

plot(subgraph.freq.otago,type="o",lty=3, xlab="Subgraph",ylab="Frequency")
points(subgraph.freq.otago.p,type="o",lty=2)

plot(subgraph.freq.otago~subgraph.freq.otago.p)
abline(a=0,b=1)

