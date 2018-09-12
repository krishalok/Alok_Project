# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Wed Oct 18 14:21:24 EDT 2017

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE15471", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00000000000000000000000000000000000000011111111111",
               "1111111111111111111111111111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=54678)
var=c("HCLS1","GP2","PRSS3P2","PRSS1","CPA1","CLPS","COL5A2","PTPRC","CTRC","HLA-DRA","CDH11","ENTPD1","PLA2G1B","PNLIPRP1","ERP27","CELA3B","TCF4","SAMSN1","SULF1","CPA2","CUZD1","SYCN","SERPINI2","HLA-DPA1","GATM","AP1S2","FBN1","LAPTM5","CELP","WIPF1","LCP2","CELF2","ATP8B4","RAB34","CELA3A","PRSS3","ARHGDIB","HLA-DMB","C3","DACT1","PNLIP","CEL","KLK1","CASP1","COL5A1","Sep-06","CLIC2","DPYSL3","MS4A6A","CPB1","CTRL","CELA2A","A2M","CD48","RAB31","GNB4","FCER1G","VCAN","MSRB3","GIMAP6","RASSF2","PCDHAC1","CELA2B","CYP1B1","B2M","LY96","EVI2A","CD300A","RNASE6","CD37","REG1A","TMEM97","NINJ2","PDIA2","SCN10A","RAC2","CMAHP","CCDC80","LOXL1","CD86","COL6A2","SULF2","EFEMP2","PRKCB","HLA-A","RNF215","PMP22","ASAP1","KCTD12","C1S","IMPA2","SPARC","TMEM52","ROBO1","COL6A3","ACRV1","HPN")
tT <- subset(tT, select=c("Gene.symbol","logFC"))

z=subset(tT[tT$Gene.symbol%in%var,])

write.table(z, file="C:/Users/alokum/Desktop/Study/GSE15471.txt", row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE15471", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("00000000000000000000000000000000000000011111111111",
               "1111111111111111111111111111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("normal","tumor")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE15471", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")
