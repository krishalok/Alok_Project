rm(list=ls())
library(igraph)
library(iplots)
library(linkcomm)


proteinList <- read.table("~/Desktop/InWeb29.txt", header=TRUE, sep="\t")
t <- as.data.frame(proteinList)
y <- graph.data.frame(t, directed=FALSE, vertices=(NULL))

#Extract the linking communities
lc <- getLinkCommunities(t, hcmethod = "average", use.all.edges=TRUE, edglim = 10^4, directed=TRUE, plot = TRUE)
#Find overlapping communities

oc <-get.community.overlaps(lc)

#Get nodes from a community
getNodesIn(lc, clusterids = 4)

integer.edgelist(t)

linkcomm2cytoscape(lc, interaction = "pp", ea = "temp.ea")

