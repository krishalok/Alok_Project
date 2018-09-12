setwd("C:\\Users\\alokum\\Desktop\\Menuscript\\PNAS\\Papers\\Main_Motif_nature\\GA_model")
setwd("C:\\Users\\alokum\\Desktop\\Menuscript\\PNAS\\Menuscript\\data")
setwd("C:\\Users\\alokum\\Desktop\\Menuscript\\PNAS\\Menuscript\\data")



rm(list=ls())
library(igraph)
library(iplots)
library(linkcomm)
library(ggplot2)
#proteinList <- read.table("PAAD_GE_Co_Network.csv", header=TRUE, sep=",")
proteinList <- read.csv("mets_final.csv", header=TRUE, sep=",")
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






