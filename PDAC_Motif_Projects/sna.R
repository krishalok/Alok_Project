##----------------------------- Social Network Analysis ------------------------------

if (!require(network)) {install.packages("network"); require(network)}  ## basic SNA stuff, relational data
if (!require(igraph)) {install.packages("igraph"); require(igraph)}  ## fancy SNA stuff, main package
if (!require(tm)) {install.packages("tm"); require(tm)}  ## text mining package


##------------------------------ SNA Data Structure -------------------------------
dat3490 <- read.csv("http://people.fas.harvard.edu/~mair/psych3490/3490list.csv")
fmName <- as.character(dat3490[,2])
fmName                        ## first and middle names           
fmsplit <- strsplit(fmName, split = " ")     ## separate first and middle name in structure
fmsplit
fName <- sapply(fmsplit, function(xx) xx[1]) ## extract first name only and store in vector
fName
N <- length(fName)      ## number of participants

## --- binary (no weights), undirected network
set.seed(123)
vec01 <- sample(0:1, N*(N-1)/2, replace = TRUE)
X <- diag(N)
X[upper.tri(X)] <- vec01
X <- X + t(X) 
diag(X) <- 0
rownames(X) <- colnames(X) <- fName
X                   ## binary, symmetric adjacency matrix
isSymmetric(X)

net1 <- graph.adjacency(X, mode = "undirected")      ## getting the data in shape in order to use it in igraph
x11()
plot(net1)

el1 <- get.edgelist(net1, names = TRUE)    ## this would be the corresponding edge list
el1

##--- weighted, directed network
set.seed(123)
veck <- rpois(N*N, lambda = 0.5)   ## create weights (drawing from a Poisson distribution)
veck
hist(veck)
veck[c(4, 20)] <- 8                ## let's create some outliers


X <- matrix(veck, ncol = N)
diag(X) <- 0
rownames(X) <- colnames(X) <- fName
X        ## adjacency matrix (not symmetric, weights)
isSymmetric(X)

## Note: John -> Arthur and Aleksandr -> Nicolas strong tie

net2 <- graph.adjacency(X, weighted = TRUE)      ## getting the data in shape in order to use it in igraph
E(net2)                         ## edges
E(net2)$weight                  ## edge weights
V(net2)                         ## vertices
x11()
plot(net2,  edge.arrow.size = 0.5, edge.width = E(net2)$weight,
     layout = layout.fruchterman.reingold)

## plot.igraph provides numerous options for customizing the plot; see
?igraph.plotting

## interactive plot
tkplot(net2,  edge.arrow.size = 0.5, edge.width = E(net2)$weight)


##------------------------------ Centrality and Prestige ------------------------
## ------- undirected networks ------
## Florentine families: During the renaissance (14th-17th century) Florence
## was the banking center of Europe (mostly due to the Medicis). During this period
## of time many aristocratic families lived in Florence. Aristocratic family
## members tend to marry members from other aristocratic families. 
## In this example we explore the Florentine marriage network and compute 
## various network measures

data(flo)
flo                ## adjacency matrix, undirected --> symmetric
isSymmetric(flo)

floG <- graph.adjacency(flo, mode = "undirected")   ## let's just convert the adjacency matrix into an igraph object

x11()
plot(floG, vertex.size = 0, edge.arrow.size = 0.25,
     vertex.label.dist = 0.5)     

## let's remove the poor Pucci's
floG1 <- delete.vertices(floG, v = "Pucci")
x11()
plot(floG1, vertex.size = 0, edge.arrow.size = 0.25,
     vertex.label.dist = 0.5)  

E(floG1)               ## edges
V(floG1)               ## vertices
N <- length(V(floG1))  ## number of nodes

## -- Nodal degree:
degvec <- sort(degree(floG1), decreasing = TRUE)
degvec
mean(degvec)     ## useful if we want to compare different marriage networks 

## -- Density:
graph.density(floG1)  ## ~19% of all possible edges are present

## -- Centrality measures:
## Degree centrality
degCvec <- sort(degree(floG1), decreasing = TRUE)  ## degree centrality
degCvec
## since the degree depends on the number of actors we can standardize 
## it by the number of nodes - 1 (e.g. for comparing across networks)
degCvec/(N-1)

centralization.degree(floG1)  ## group degree centralization 

## Closeness centrality
sort(closeness(floG1), decreasing = TRUE)   ## closeness centrality
shortest.paths(floG1)     ## matrix with geodesics

centralization.closeness(floG1)

## Betweenness centrality
sort(betweenness(floG1, directed = FALSE, normalized = TRUE), decreasing = TRUE)  ## node betweenness
centralization.betweenness(floG1, directed = FALSE)  ## aggregate measure

ebet <- edge.betweenness(floG1, directed = FALSE) ## edge betweenness  
ebet
el <- get.edgelist(floG1)
el
data.frame(from = el[,1], to = el[,2], betweenness = ebet)

## ----- directed networks -----
## Let's now look at some directed trade data between countries
if (!require(SNAData)) {install.packages("SNAData"); require(SNAData)}  ## datasets
data(basicGoods)     ## it's a graphNEL object (general object structure for graphs) 

tradenet <- igraph.from.graphNEL(basicGoods)  ## convert it to an igraph object

x11()
plot(tradenet, vertex.size = 0, edge.arrow.size = 0.4,
     vertex.label.dist = 0.5)

get.adjacency(tradenet)   ## adjacency matrix
get.edgelist(tradenet)

## --- Indegree and Outdegree
## Degree Prestige (same as indegree)
sort(igraph:::degree(tradenet, mode = c("in")), decreasing = TRUE)  ## indegree (imports)

## Outdegree 
sort(igraph:::degree(tradenet, mode = c("out")), decreasing = TRUE) ## outdegree (exports)

## --- Centrality/Prestige
## Centrality as above (takes outgoing edges) 
sort(closeness(tradenet), decreasing = TRUE)   ## exports
## Prestige (takes ingoing edges)
sort(closeness(tradenet, mode = "in"), decreasing = TRUE)    ## imports


##-------------------------- Cohesive Subgroups -------------------------------
load(url("http://people.fas.harvard.edu/~mair/psych3490/statementsGOP.rda"))  

inspect(gopCorp.unique)   ## 254 statements

inspect(tm_filter(gopCorp.unique, FUN = function(x) any(grep("sanctity", x))))   ## as an example

##--- re-organize the data as edge list
gopCorp.unique1 <- tm_map(gopCorp.unique, tolower)
mystopwords <- c("beleive", "shld", "1", "-", "wenot", "conservatismthe", "etc", "im",
                 "fatherthe", "conservativebelieve", "governmentprolife2nd",
                 "amendmentand", "valuessmall", "ive", "4", "familyrepublican",
                 "-government", "1st", "believe", "belive", "still", "dont", 
                 "want","seen", "b", "w","can")
statementsGOP <- Corpus(VectorSource(lapply(gopCorp.unique1, removeWords, 
                                            c(mystopwords, stopwords("english")))))
slen <- unlist(tm_map(statementsGOP, function(tf) sum(termFreq(tf))))  ## statement length
statementsGOP <- statementsGOP[which(slen > 1)]                        ## select the statements with more than 1 word
statementsGOP

## function to get pairwise word structure suited for directed network graph
textsna <- function(str1) {
  splitState <- scan_tokenizer(str1)
  indmat <- matrix(c(1, rep(2:length(splitState), each = 2)), ncol = 2, 
                   nrow = (length(splitState) - 1))
  matrix(splitState[indmat], ncol = 2, byrow = TRUE)
}  
statementPairs <- do.call(rbind, lapply(statementsGOP, textsna))   ## warnings can be ignored

statementsGOP[[1]]
statementPairs[1:11,]   ## in row 11 the 2nd statement starts
## ----- end data preparation

## ------ create network and produce network plot
statementGraph <- graph.edgelist(statementPairs, directed = TRUE) ## make a graph with each row entry as a from-to vertex and an edge
E(statementGraph)$weight <- count.multiple(statementGraph)        ## edge frequencies -- count how often an edge has exactly the same tail and edge vertices  
statementGraph <- simplify(statementGraph, edge.attr.comb = 
                             list(weight = max, name = "concat", "ignore") ) ## keep multiplicity as edges
E(statementGraph)$weight

set.seed(123)
ly <- layout.fruchterman.reingold(statementGraph)    ## define layout (more concentric)
## ly <- layout.kamada.kawai(statementGraph)         ## more branches going out of the center


x11()
plot(statementGraph, layout = ly, vertex.size = 2, edge.color = "lightgray", 
     edge.arrow.size = 0.05, edge.curved = FALSE, 
     vertex.label = NA, 
     rescale = FALSE, 
     xlim = c(-530, 630), ylim = c(-440, 600),
     main= "Republican Statements (Large Communities)", asp = 0,
     margin = -1)


##--------- compute cliques ------------
cliqList <- cliques(statementGraph, min = 4)   ## find and all cliques (at least 4 nodes)
cliqList                   ## returns the node index only
V(statementGraph)$name           ## gives the names
cliqlistN <- lapply(cliqList, function(xx) {
  V(statementGraph)$name[xx]
})
cliqlistN

## let's extract the largest clique and produce a plot
cliqLarge <- largest.cliques(statementGraph)
cliqLarge
statementGraphSub <- induced.subgraph(statementGraph, cliqLarge[[1]])
x11()
plot(statementGraphSub, vertex.size = 2, edge.color = "gray",
     vertex.color= "black", vertex.label.color = "blue",
     edge.arrow.size = 1, edge.curved = FALSE, edge.width = E(statementGraphSub)$weight,
     vertex.label.dist = 0.5, vertex.label.cex = 0.8, vertex.label.font = 2)



## ------ compute communities ---------
## A nice blog entry that decribes various community approaches can be found here:
## browseURL("http://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph")

## We use the wakltrap algorthim
xc <- walktrap.community(statementGraph, weights = E(statementGraph)$weight, step = 6)   ## very quick
membership(xc)        ## gives the community membership for each node

edgeList <- tapply(seq_along(membership(xc)), membership(xc), function(xx) xx) ## communities without labels
edgeList                                ## communities as list (containing node index)
comList <- tapply(membership(xc), membership(xc), names)  
comList                                 ## communities as list (containing node names)
length(comList)                         ## number of communities
comsize <- sapply(comList, length)
table(comsize)          
barplot(table(comsize), main = "Community Size")   ## frequency distribution of community sizes

bigComIndex <- which(comsize > 20) ## big communities (having more than 20 nodes)
bigComIndex 
comList[bigComIndex]

coms3 <- edgeList[bigComIndex[[1]]]    ## pull out first community
coms4 <- edgeList[bigComIndex[[3]]]    ## pull out third community


## ------- produce a fancy plot: full network, big communities added
x11()
vcolo <- c(hcl(h=0,35,60),hcl(h=72,c=35,l=60),hcl(h=144,c=35,l=60),hcl(h=216,c=35,l=60),hcl(h=320,c=35,l=60))
op <- par(mar = c(1,2,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(statementGraph, layout = ly, vertex.size = 2, edge.color = "lightgray", 
     edge.arrow.size = 0.05, edge.curved = FALSE, 
     vertex.label = NA, 
     rescale = FALSE, 
     xlim = c(-530, 630), ylim = c(-440, 600),
     main= "Republican Statements (Large Communities)", asp = 0,
     margin = -1)
comGraphBig1 <- induced.subgraph(statementGraph, coms3[[1]])
ly1 <- ly[unlist(coms3[[1]]),] 
plot(comGraphBig1, layout = ly1, vertex.size = 3, edge.color = vcolo[2],
     edge.arrow.size = 0.2, edge.curved = FALSE,  vertex.label = NA, 
     vertex.label.dist = 0, vertex.label.cex = 0.8, vertex.label.font = 2,
     rescale = FALSE, add=TRUE, vertex.color=vcolo[2],vertex.frame.color=vcolo[2],vertex.label.color=vcolo[2])
comGraphBig2 <- induced.subgraph(statementGraph, coms4[[1]])
ly1 <- ly[unlist(coms4[[1]]),]
plot(comGraphBig2, layout = ly1, vertex.size = 3, edge.color = vcolo[3], 
     edge.arrow.size = 0.2, edge.curved = FALSE,  vertex.label = NA, 
     vertex.label.dist = 0, vertex.label.cex = 0.8, vertex.label.font = 2,
     rescale = FALSE, add=TRUE, vertex.color=vcolo[3],vertex.frame.color=vcolo[3],vertex.label.color=vcolo[3])
op2 <- par(mar = c(4,1,1,1))
colo <- c(hcl(72,35,40),hcl(72,35,80))
comGraphBig1 <- induced.subgraph(statementGraph, unlist(edgeList[bigComIndex[1]]))
ly1 <- ly[unlist(edgeList[bigComIndex[1]]),] 
plot(comGraphBig1, layout = ly1, vertex.size = 2, edge.color = colo[2],
     vertex.color=colo[1],vertex.label.color=colo[1],
     edge.arrow.size = 0.20, edge.curved = FALSE, edge.width = E(comGraphBig1)$weight/2,
     vertex.label.dist = 0, vertex.label.cex = 0.8, vertex.label.font = 2,
     rescale = FALSE, xlim=range(ly1[,1]), ylim=range(ly1[,2]))
colo <- c(hcl(144,35,40),hcl(144,35,80))
comGraphBig2 <- induced.subgraph(statementGraph, unlist(edgeList[bigComIndex[3]]))
ly1 <- ly[unlist(edgeList[bigComIndex[3]]),] 
plot(comGraphBig2, layout = ly1, vertex.size = 2, edge.color = colo[2],
     vertex.color=colo[1],vertex.label.color=colo[1],
     edge.arrow.size = 0.20, edge.curved = FALSE, edge.width = E(comGraphBig2)$weight/2,
     vertex.label.dist = 0, vertex.label.cex = 0.8, vertex.label.font = 2,
     rescale = FALSE, xlim=range(ly1[,1]), ylim=range(ly1[,2]))
par(op)
par(op2)

