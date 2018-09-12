source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

Sys.putenv("http\_proxy" = "http://my.proxy.org:9999")

library(biomaRt)
listMarts()
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

attributes = listAttributes(ensembl)
attributes
#mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

genes<- c(read.table("Desktop/p2.txt"))

g <- getBM( attributes= c("hgnc_symbol"), values=genes,mart= ensembl, uniqueRows=FALSE)

#dupRows <- union(which(duplicated(g[,1])), which(duplicated(g[,2])))
#ensemblIDs <- g[-dupRows, 2]
#names(ensemblIDs) <- g[-dupRows, 1]

write.csv(g, "prot2.csv") 
 

