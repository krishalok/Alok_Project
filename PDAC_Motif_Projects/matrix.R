
####################################################################################
setwd("C:\\Users\\alokum\\Desktop\\Menuscript\\PNAS\\Menuscript\\ADEX")
dat1 <- read.table("ids.txt", sep = "\t", header = T, stringsAsFactors = F)
dat2 <- read.table("nature16965-s2.txt", sep = "\t", header = T, stringsAsFactors = F)

library(dplyr)
dat2 <- dat2 %>%
  select(noquote(order(colnames(dat2))))

dat1 <- dat1 %>%
  arrange(icgc_id)

#colnames(dat2) <- c(dat1$membership.ordered, "X")

#names(dat2) <- dat1$icgc_id[match(names(dat2), dat1$icgc_id)]
names(dat2) <- plyr::mapvalues(names(dat2), from = dat1$icgc_id, to = dat1$membership.ordered)

dat3 <-read.table("ENSB2HUGO.txt", sep = "\t", header = T, stringsAsFactors = F)

dat2$X <- plyr::mapvalues(dat2$X, from = dat3$ensembl_gene_id, to = dat3$hgnc_symbol)

which(colnames(dat2)=="1")
which(colnames(dat2)=="2")
which(colnames(dat2)=="3")
which(colnames(dat2)=="4")

varnames_1 <- c(97,1,7,9,11,13,22,26,40,49,50,73,74,76,80,85,88)
varnames_2 <- c(97,2,4,18,19,25,28,29,31,32,34,39,41,42,43,52,54,63,64,65,68,79,81,82,93,95)
varnames_3 <- c(97,3,5,6,10,15,21,30,33,37,44,46,48,51,58,62,69,70,71,72,77,78,91,92,94,96)
varnames_4 <- c(97,8,12,14,16,17,20,23,24,27,35,36,38,45,47,53,55,56,57,59,60,61,66,67,75,83,84,86,87,89,90)

dat_1 <- dat2[,varnames_1]
dat_2 <- dat2[,varnames_2]
dat_3 <- dat2[,varnames_3]
dat_4 <- dat2[,varnames_4]
#1=ADEX
#2=Immunogenic
#3=Squamous
#4=Pancreatic Progenitor
ADEX<-dat_1
Immunogenic<-dat_2
Squamous<-dat_3
Pancreatic_Progenitor<-dat_4
library(data.table) 
fwrite(ADEX, "ADEX.csv")
fwrite(Immunogenic, "Immunogenic.csv")
fwrite(Squamous, "Squamous.csv")
fwrite(Pancreatic_Progenitor, "Pancreatic_Progenitor.csv")
































































