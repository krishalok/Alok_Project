library(igraph)
setwd("/Users/BC/Documents/spinglass/")
genes <- read.table("data/cc.txt", header=FALSE, sep="\t")
g=as.matrix(genes)
y <- graph.edgelist(g, directed=FALSE)#, weights=E(y)$weight)
snps <- read.table("~/Documents/NeXTProject/combined_MTLasso_mEQTL_eqtls.txt", header=TRUE, sep="\t")
snps <- as.data.frame(snps)
snps <- snps$SNP
snps
y <- simplify(y)

bad.vs<-V(y)[degree(y) == 1] 

# remove isolated nodes
y <-delete.vertices(y, bad.vs)
is.connected(y)
diameter(y, directed = FALSE, unconnected = FALSE, weights = NULL)
get.diameter(y, directed = FALSE, unconnected = FALSE, weights = NULL)

sg1 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg2 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg3 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg4 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg5 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg6 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg7 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg8 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg9 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)
sg10 <- spinglass.community(y, weights = NULL, update.rule="config",start.temp=1, stop.temp=.80, gamma =2)




#Compare community structures
comp1.10 <- compare(membership(sg1), membership(sg10), method="nmi")
comp1.9 <- compare(membership(sg1), membership(sg9), method="nmi")
comp1.8 <- compare(membership(sg1), membership(sg8), method="nmi")
comp1.7 <- compare(membership(sg1), membership(sg7), method="nmi")
comp1.6 <- compare(membership(sg1), membership(sg6), method="nmi")
comp1.5 <- compare(membership(sg1), membership(sg5), method="nmi")
comp1.4 <- compare(membership(sg1), membership(sg4), method="nmi")
comp1.3 <- compare(membership(sg1), membership(sg3), method="nmi")
comp1.2 <- compare(membership(sg1), membership(sg2), method="nmi")

# x <- which(farthest.nodes(y)<=7)
# x

sizes(sg)
membership(sg)
communities(sg)

#layout <- layout.kamada.kawai(y)
g <- induced.subgraph(graph=y,vids=unlist(neighborhood(graph=y,order=7,nodes=snps)))
g <- get.edgelist(g)
head(g)
#plot(g, vertex.color=membership(sg), layout=layout,vertex.size=9, 
#     vertex.label=V(y)$name, vertex.frame.color="red", edge.color="grey",
#     edge.arrow.size=0.01, rescale=TRUE,vertex.label=TRUE, 
#     vertex.label.dist=0.0,vertex.label.cex=0.1, add=FALSE,vertex.label.font=1)
edgelist <- tapply(seq_along(membership(sg3)), membership(sg3), function(xx) xx)
comList <- tapply(membership(sg3), membership(sg3), names)

length(comList)                         ## number of communities
comsize <- sapply(comList, length) #find out the size of each community
comsize
#vx <- induced.subgraph(neighborhood.size(y,order=8, nodes=snps))
#V(y)$name <- V(y)$name 
#vx<- graph.neighborhood(vx, 8, snps)
#vx<-get.edgelist(vx)
#vx

#Find out which communities are greater than 7 nodes  
bigComIndex <- which(comsize>7)
bigComIndex
h <- comList[bigComIndex] # Extract the membership withe the node names in each community

coms3 <- edgelist[bigComIndex[[1]]]    ## pull out first community
coms3
#coms4 <- edgelist[bigComIndex[[3]]]    ## pull out third community
comGraphBig1 <- induced.subgraph(y, coms3[[1]])
#ly1 <- layout[unlist(coms3[[2]]),] 
cc <-get.edgelist(comGraphBig1) #get the edgelist for the selected community
cc
write.csv(cc, "results/sg1_comm1.csv")

#par(mfrow=c(3,3))
#plot(dg[[11]], layout=layout.fruchterman.reingold,vertex.color="green", vertex.size=22, vertex.label.cex=1.5,main = "Community 11")
#plot(dg[[44]], layout=layout.fruchterman.reingold,vertex.color="red", vertex.size=22, vertex.label.cex=1.5, main = "Community 44")
#plot(dg[[43]], layout=layout.fruchterman.reingold,vertex.color="yellow", vertex.size=22, vertex.label.cex=1.5,main = "Community 43")
#plot(dg[[7]], layout=layout.fruchterman.reingold,vertex.color="cyan", vertex.size=22, vertex.label.cex=1.5, main = "Community 7")
#plot(dg[[20]], layout=layout.fruchterman.reingold,vertex.color="purple", vertex.size=22, vertex.label.cex=1.5, main = "Community 20")
#plot(dg[[17]], layout=layout.fruchterman.reingold,vertex.color="blue", vertex.size=22, vertex.label.cex=1.5,main = "Community 17")


  
  
  #Test walktrap community detection
  
wt1 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt2 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt3 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt4 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt5 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt6 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt7 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt8 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt9 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
wt10 <- walktrap.community(y, weights = NULL, membership=T, modularity=T)
sizes(wt1)
membership(wt1)
communities(wt1)
edgelist <- tapply(seq_along(membership(wt1)), membership(wt1), function(xx) xx)
comList <- tapply(membership(wt1), membership(wt1), names)
length(comList)                         ## number of communities
comsize <- sapply(comList, length)
comsize
bigComIndex <- which(comsize>7)
bigComIndex
#bigComIndex <- which(neighborhood(graph=y,order=7,nodes=snps))
h <- comList[bigComIndex]
h
coms3 <- edgelist[bigComIndex[[1]]]    ## pull out first community
coms3
#coms4 <- edgelist[bigComIndex[[3]]]    ## pull out third community
comGraphBig1 <- induced.subgraph(y, coms3[[1]])
#ly1 <- layout[unlist(coms3[[2]]),] 
cc <-get.edgelist(comGraphBig1)
cc
write.csv(cc, "results/walktrap/run 1/wt1_comm1.csv")


