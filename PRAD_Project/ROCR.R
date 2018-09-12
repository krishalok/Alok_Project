setwd("C:\\Users\\alokum\\Desktop\\Data-20170612T233518Z\\Data")

dd=read.table("3-pca.txt",sep="\t", header=TRUE)
##var=c("Samples","MYBPC1","DHFR","EYA1","RPE65","FUCA1","ALOX15B","CSE1L","DRAP1","BCAS1","MCM4","SMC4","NEK2","TROAP","SPAG5","MKI67","CENPE","PYGM","ASPN","MELK","UBFD1","PGM5","ALDH1A2","SELE","TGFB3","GNE","GRIA3","FLNA","ARG2","FMOD","SCN7A")
##var=c("Samples","ALDH1A2","ASPN","EYA1","FUCA1","MELK","MKI67","MYBPC1","NEK2","PGM5","SCN7A","SPAG5")
##var=c("Samples","MYL9","S100A4","CNN1","MYLK","MYH11","MT1F","DBT","LUM","EYA1","PKM2","C20orf24","SORBS1","MYBPC1","ALOX15B","LDHA","ATP5B","ARG2","GARS","KLF5","GMDS","HMGN3","USP9X","STX3","FGF7","UBAP2L","GNE","SPG20","GSTP1","FLNA","ZNF3","SSPN","FMOD","CSE1L","ENAH","DKK3","SMC4","GNG4","MEIS2","CYP27A1","TGFB3")
#var=c("Samples","MYL9","S100A4","CNN1","MYLK","MYH11","MT1F","DBT","LUM","EYA1","PKM2","C20orf24","SORBS1","MYBPC1","ALOX15B","LDHA","ATP5B","ARG2","GARS","KLF5","GMDS","HMGN3","USP9X","STX3","FGF7","UBAP2L","GNE","SPG20","GSTP1","FLNA","ZNF3","SSPN","FMOD","CSE1L","ENAH","DKK3","SMC4","GNG4","MEIS2","CYP27A1","TGFB3","FUCA1","DPT","VIPR1","DRAP1","ASPN","XRCC2","PTTG1","SLC15A2","SLC39A8","BCAS1","SCN7A","MAOA","PTGER2","CTHRC1","FMO5","DHFR","TYMS","MT1H","CHRNA2","MS4A2","LMNB1","GSTM3","MCM4","GRIA3","NEK2","TPX2","C2","PGM5","KIF14","TROAP","TCF21","ESRRG","CEACAM1","ALDH2","BMP6","ALDH1A2","NUSAP1","MT1E","UBFD1","TOP2A","SPAG5","CENPI","UBE2C","TRIP13","ESM1","ASPA","RPE65","PLAGL2","CCNA2","CENPF","ANGPTL2","CKMT2","NCAPH","STMN1","ARID5A","GTSE1","ATP2A1","PKMYT1","MMP11","MKI67","MELK","PRR4","ZNF212","SELE","CENPA","RFPL1","PYGM","CCK","ESPL1","LGALS3","CENPE","PYGL","IRF5","RAD54L","CDC25C","GCGR","PAH")
var=c("Samples","ATP5B","NUSAP1","XRCC2","ESPL1","CENPI","RPE65","GCGR","STMN1","LDHA","MMP11","MYH11","ZNF212","C20orf24","PKM2","CSE1L","SPAG5","NCAPH","DRAP1","GARS","CENPE","ATP2A1","CENPA","GTSE1","RFPL1","CENPF","PKMYT1","RAD54L","EYA1","FUCA1","FLNA","MYL9","ARG2","PTTG1","CEACAM1","DHFR","KIF14","CTHRC1","MT1E","UBFD1","CKMT2","TROAP","MAOA","MYLK")
####LOOCV
library(caret)
library(e1071)
library(pROC)
##23-genes

set.seed(123)


data=dd[var]
ctrl <- trainControl(method = "loocv", savePred=T, classProb=T)
svm.fit <- train(Samples~., data,method='svmPoly',cost=10, trControl=ctrl)
gbm.probs <- predict(svm.fit, data,type="prob")
gbm.ROC <- roc(predictor=gbm.probs$normal,
               response=data$Samples, 
               levels=rev(levels(data$Samples)))

pdf("3.pdf")
plot(gbm.ROC,print.auc=TRUE,main="3-LOOCV",  print.auc.y = 0.8, print.auc.x = 0.65, xlim=c(0, 1), ylim=c(0, 1))
dev.off()

sink("3.txt")
tune.out <- tune(svm, Samples~., data=data, kernel='polynomial',
                 ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100)))
yhat <- predict(tune.out$best.model, data)
confusionMatrix(yhat,data$Samples)
sink()

