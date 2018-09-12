

##### from GPL
getGEOdataObjects <- function(x, getGSEobject=FALSE){
  # Make sure the GEOquery package is installed
  require("GEOquery")
  # Use the getGEO() function to download the GEO data for the id stored in x
  GSEDATA <- getGEO(x, GSEMatrix=T, AnnotGPL=FALSE)
  # Inspect the object by printing a summary of the expression values for the first 2 columns
  print(summary(exprs(GSEDATA[[1]])[, 1:2]))
  
  # Get the eset object
  eset <- GSEDATA[[1]]
  # Save the objects generated for future use in the current working directory
  save(GSEDATA, eset, file=paste(x, ".RData", sep=""))
  
  # check whether we want to return the list object we downloaded on GEO or
  # just the eset object with the getGSEobject argument
  if(getGSEobject) return(GSEDATA) else return(eset)
}
# Store the dataset ids in a vector GEO_DATASETS just in case you want to loop through several GEO ids
GEO_DATASETS <- c("GSE73072")

# Use the function we created to return the eset object
eset <- getGEOdataObjects(GEO_DATASETS[1])
# Inspect the eset object to get the annotation GPL id
eset 

# Get the annotation GPL id (see Annotation: GPL10558)
gpl <- getGEO('GPL14604', destdir=".")
Meta(gpl)$title

# Inspect the table of the gpl annotation object
colnames(Table(gpl))

# Get the gene symbol and entrez ids to be used for annotations
Table(gpl)[1:4, c(1, 2, 3, 4)]
dim(Table(gpl))

# Get the gene expression data for all the probes with a gene symbol
geneProbes <- which(!is.na(Table(gpl)$ID))
probeids <- as.character(Table(gpl)$ID[geneProbes])

probes <- intersect(probeids, rownames(exprs(eset)))
length(probes)

geneMatrix <- exprs(eset)[probes, ]

inds <- which(Table(gpl)$ID %in% probes)
# Check you get the same probes
head(probes)
head(as.character(Table(gpl)$ID[inds]))

# Create the expression matrix with gene ids
geneMatTable <- cbind(geneMatrix, Table(gpl)[inds, c(1, 2, 3, 4)])
head(geneMatTable)

# Save a copy of the expression matrix as a csv file
#write.csv(geneMatTable, paste(GEO_DATASETS[1], "_DataMatrix.csv", sep=""), row.names=T)

virus.genes <-read.table(file = 'RF_train_104_sub_0_AvgOnset_VarIMP_2.tsv', sep = '\t', header = TRUE)
virus.genes1<-virus.genes[1:17,]
res1 <-match(virus.genes1$Gene_ID,geneMatTable$ID)
res1[!is.na(res1)]
new <- merge(geneMatTable,virus.genes1, by.x='ID', by.y='Gene_ID',) 


######add gene symbols 
library(org.Hs.eg.db)
library(annotate)
##Testing
e2s = toTable(org.Hs.egSYMBOL)
lookUp('11104', 'org.Hs.eg', 'SYMBOL') 
getSYMBOL(c('3815', '3816', '2341'), data='org.Hs.eg')
final_Gene <- merge(new,e2s, by.x='GENE_ID', by.y='gene_id')
write.table(final_Gene,file = "final-gene.txt",sep='\t')
