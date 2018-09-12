proteinList_normal <- read.csv("Normal_pancrease.csv", header=TRUE, sep=",")
#proteinList <- read.table("ADEX_out_0.05.csv", header=TRUE, sep=",")
#proteinList <- read.table("Immunogenic_out_0.05.csv", header=TRUE, sep=",")
#proteinList <- read.table("Pancreatic_Progenitor_out_0.05.csv", header=TRUE, sep=",")
#proteinList <- read.table("Squamous_output_0.05.csv", header=TRUE, sep=",")
t_n <- as.data.frame(proteinList_normal)
#############Calculate Z-scare#############
#################https://stackoverflow.com/questions/6148050/creating-z-scores
#t$Correlation<-scale(t$Correlation)
t_n <- t_n[,1:3]
y_n <- graph.data.frame(t_n, directed=FALSE, vertices=(NULL))
