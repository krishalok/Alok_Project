library(tidyr)
setwd("C:/Users/alokum/Desktop/Menuscript/Oncotarget/PRAD_Predict/temp/heatmap")
var=c("IL1RAPL2","PHOX2A","EIF2B2","FOXS1","PIGW","ALPL","GPR4","FHOD1","EVI2B","CCL27","ALPK1","TOMM34","B3GNT9","RAX2","S100A3","AZI1","PDLIM4","MEIS3","NCRNA00115","YY2","CCDC39","RACGAP1P","FAM76A","C18orf10","KLHL30","PRPF39","FBXW8","SLC15A3","SLC25A12","ALG9","TEX101","EFCAB10","TMEM175","SLC16A14","CARNS1","CMPK2","RNGTT","TMEM26","PIK3CD","ARL11","LOC283663","ANXA10","ZNF611","PLAG1","FAT4","PSPC1","CHIC2","GUCY2D","SLC27A3","LPAR1","TBX15","TOP3A","FOXR1","C5orf25","KLHDC8B","BCL10","RBMY2EP","OVCA2","KIAA0355","FUT1","BBS12","ST7","ARHGAP11A","MTERF","THBS3","CIITA","SLC35F4","ABCB4","GJD3","CD86","EFNA4","ANTXRL","ANKRD2","FOLR3","LRRC70","MAL","SIKE1","OPN3","DAPP1","FAM115A","SAMD11","RAB11FIP2","KIF18B","ZNF564","ZNF404","NFATC4","IL16","CBR4","MS4A7","ROPN1","NEDD1","TWIST2","GDF7","SLC41A2","ENTPD2","ASPM","APBB1","FTO","LPA","MRAS","BTLA","FOXD3","ICAM4","LOC152217","DCAF4L2","PRLHR","UGT1A6","UCK2","SRPK3","IL1F8","CHI3L1","ORMDL1","FAM55B","PYROXD1","METTL13","HOXB9","ZNF507","HOXB4","KLF1","KIAA0415","TLL1","BTBD16","DUSP8","AP3M2","ALOX5AP","ACER1","MYF6","SPRED1","FGD1","PXMP2","KDELR3","JRK","DAPK2","RCAN1","OR8A1","RPL10L","NCEH1","PTPLB","CHRNB1","CYP7A1","RIPPLY2")
matrix=read.table("PRAD_genomicMatrix",sep = "",header = TRUE)
clinical=read.table("PRAD_Gleason.txt",sep = "\t",header=TRUE)
#data<-data[rowSums(data=="")!=ncol(data), ]
data=matrix[var,]
ind <- match(names(data),clinical$sampleID)
names(data) <- clinical$gleason_score[ind]
unique(names(data))
table(names(data))
library(som)
#esetSel=normalize(data)
library(gplots)
pdf("PRAD_141_HM.pdf", height=25, width=30)
heatmap.2(as.matrix(data),ColSideColors=c(rep("green",49),rep("blue",277),rep("red",195),rep("blue",29)),col=redgreen(20),scale="none",key=TRUE,symkey=FALSE,density.info="none",trace="none",cexRow=1.4,dendrogram="row",Colv = F, main="TCGA-PRAD 141 Genes expression profile GS6-GS7-GS[8-10]")
#heatmap.2(as.matrix(esetSel),ColSideColors=c(rep("green",49),rep("blue",277),rep("red",195),rep("blue",29)),col=redgreen(20),breaks = seq(-1,1,0.1),scale="none",key=TRUE,symkey=FALSE,density.info="none",trace="none",cexRow=1.4,dendrogram="row",Colv = F, main="TCGA-PRAD 141 Genes expression profile GS6-GS7-GS[8-10]")
dev.off()
#boxplot for each data sets
write.csv(data,file="boxplot_data.csv")

