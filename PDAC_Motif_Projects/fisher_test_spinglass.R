#Test the hypothesis that the geneset in the previously published genesets in the pathways are independent of the enriched pathways identified in my results.
#If we obtain a low p-value (<<.05), we can reject the null hypothesis.

setwd("/Users/BC/Documents/spinglass/genesets/")

path <- "../genesets/"
file.names <- list.files(path, pattern=".txt")

my_gene_sets <- list.files(path, pattern = "_pathway.txt")
#my_func <- function(file.names {
  df <- as.data.frame(read.table(file.names, header= T, sep="\t"))
  my_df <- as.data.frame(read.table(my_gene_sets, header= T, sep="\t"))
  for (i in 1:length(file.names)){
    for (j in 1:length(my_gene_sets)){
      matches<- merge(df,my_df, by.y="ID", by.x=df[,1], all.x = T)
      common <- nrow(matches)
      total <- nrow(my_df) + nrow(df) - common
      my_set <- nrow(my_df) - common
      df_count <- nrow(df) - common
      nrows <- sapply(file.names, function(f) nrow(file.names))
      not_in_either <- sum(nrows)-nrow(df)
      c <- matrix(c(common,my_set,df_count,not_in_either), nrow=2)
      all_tables <- rbind(c[i])
    }
  }
#}

#(lapply(file.names, my_func))


my_notch_pway <- as.data.frame(read.table("comm5_notch_matching_pathway.txt", header=T, sep="\t"))
my_basal_pway <- as.data.frame(read.table("comm13_basal_pathway.txt", header=T, sep="\t"))

publ_notch <- as.data.frame(read.table("notch_geneset.txt", header=T,sep="\t"))
publ_basal <- as.data.frame(read.table("basal_geneset.txt", header=T, sep="\t"))


#Number of genes in both sets
notch_matches <- merge(publ_notch, my_notch_pway, by.y = "ID", by.x ="Notch.Geneset", all.x = T)#my_notch_pway$ID[my_notch_pway$ID %in% publ_notch$Notch.Geneset,]
common_genes <- nrow(notch_matches)

#Combine the data frames
total_genes <- nrow(my_notch_pway) + nrow(publ_notch) - common_genes

#Number only in my geneset
my_count <- nrow(my_notch_pway)- common_genes 

#Number in ref but not in my set
ref_set <- nrow(publ_notch) - common_genes #

nrows <- sapply(file.names, function(f) nrow(read.table(f, header=T, sep="\t")) )
not_in_either <- sum(nrows)- my_count - ref_set #297


c <- matrix(c(common_genes,my_count,ref_set,not_in_either), nrow=2)#,97,202,439,696584)
f <- fisher.test(c)

####################################Basal Pathway#############################
my_notch_pway <- as.data.frame(read.table("comm5_notch_matching_pathway.txt", header=T, sep="\t"))
my_basal_pway <- as.data.frame(read.table("comm13_basal_pathway.txt", header=T, sep="\t"))

publ_notch <- as.data.frame(read.table("notch_geneset.txt", header=T,sep="\t"))
publ_basal <- as.data.frame(read.table("basal_geneset.txt", header=T, sep="\t"))


#Number of genes in both sets
basal_matches <- merge(publ_basal, my_basal_pway, by.y = "ID", by.x ="Basal.Geneset", all.x = T)#my_basal_pway$ID[my_basal_pway$ID %in% publ_basal$basal.Geneset,]
common_genes <- nrow(basal_matches)

#Combine the data frames
total_genes <- nrow(my_basal_pway) + nrow(publ_basal) - common_genes

#Number only in my geneset
my_count <- nrow(my_basal_pway)- common_genes 

#Number in ref but not in my set
ref_set <- nrow(publ_basal) - common_genes #

nrows <- sapply(file.names, function(f) nrow(read.table(f, header=T, sep="\t")) )
not_in_either <- sum(nrows)- my_count - ref_set #297


c <- matrix(c(common_genes,my_count,ref_set,not_in_either), nrow=2)
f <- fisher.test(c)

#Perform Fisher Exact test
apply(all_tables,1, function(x) fisher.test(matrix(x,nr=2))$p.value)
c <- matrix(c(7,22,0,297), nrow=2)
f <- fisher.test(c)
