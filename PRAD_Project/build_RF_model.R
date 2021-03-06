#Install necessary packages
library(forestFloor)
install.packages("randomForest")
install.packages("ROCR")
install.packages("Hmisc")
source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
library(rgl)

source("https://bioconductor.org/biocLite.R")
biocLite("genefilter",dependencies=TRUE)

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
library(mlbench)
library(caret)
library(forestFloor)

#Set working directory and filenames for Input/output
setwd("C:\\Users\\alokum\\Desktop\\PRAD\\Code_Working\\RF_Predictor_v2\\PRAD_Final")

datafile="PRAD_DATA.txt" 
clindatafile="trainset_clindetails.txt"

outfile="trainset_RFoutput_2.txt"
varimp_pdffile="trainset_varImps_2.pdf"
MDS_pdffile="trainset_MDS_3.pdf"
ROC_pdffile="trainset_ROC_2.pdf"
case_pred_outfile="trainset_CasePredictions_2.txt"
vote_dist_pdffile="trainset_vote_dist_2.pdf"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep=" ")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GEO.asscession.number"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[2:length(colnames(data_import))])+1 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1,data_order)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata)

#If there are predictor variables that are constant/invariant, consider removing them
#Preliminary gene filtering
X=rawdata[,2:length(header)]
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
#At least 20% of samples should have raw intensity greater than 100 
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_Data=rawdata[filt,] 

#Get potential predictor variables
predictor_data=as.data.frame(t(filt_Data[,2:length(header)])) #Filtered
#predictor_data<-model.matrix(~ ., data = predictor_data1)
predictor_names=c(as.vector(filt_Data[,1])) #gene symbol
colnames(predictor_data)=predictor_names

#Get target variable and specify as factor/categorical
target= clindata[,"relapse..1.True."]
target[target==0]="NoRelapse"
target[target==1]="Relapse"
target=as.factor(target)

#Run RandomForests
#NOTE: use an ODD number for ntree. When the forest/ensembl is used on test data, ties are broken randomly.
#Having an odd number of trees avoids this issue and makes the model fully deterministic
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)
rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 25001, proximity=TRUE, sampsize=sampsizes)
#rf_output=lda(x=predictor_data, y=target, importance = TRUE, ntree = 25001, proximity=TRUE, sampsize=sampsizes)
varImpPlot(rf_output)
varImpPlot(rf_output,type=2)
importance(rf_output)
importanceOrder=order(-rf_output$importance)
#names=rownames(rf_output$importance)[importanceOrder][1:15]
#par(mfrow=c(5, 3), xpd=NA)
#for (name in names)
 # +   partialPlot(rf_output, df, eval(name), main=name, xlab=name,ylim=c(-.2,.9))
#######################################s
#Save RF classifier with save()
save(rf_output, file="RF_model")

#Load saved model - this will save time if re-running and you want to skip RF run
load("RF_model")

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy=1-overall_error
class1_error=paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error=paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy=100-overall_error

#Prepare stats for output to file
sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")
misclass_1=paste(confusion[1,2], rownames(confusion)[1],"misclassified as", colnames(confusion)[2], sep=" ")
misclass_2=paste(confusion[2,1], rownames(confusion)[2],"misclassified as", colnames(confusion)[1], sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Print results to file
write.table(rf_importances[,4],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2), file=outfile, append=TRUE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()
###### 10th sep 2017
plot(rf_output, main = "Variable Importance of Top 30", top = 30)

#Produce MDS plot
pdf(file=MDS_pdffile)
target_labels=as.vector(target)
target_labels[target_labels=="NoRelapse"]="N"
target_labels[target_labels=="Relapse"]="R"
target_labels[target_labels=="NA"]="NA"

MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")
#MDSplot(rf_output, target, k=3, xlab="", ylab="", zlab="", pch=target_labels, palette=c("red", "blue","green"), main="MDS plot")

dev.off()

#Create ROC curve plot and calculate AUC
#Can use Relapse/non-relapse vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling relapse vs non-relapse
predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#### added on 10-09-2017
plot(predictions)
plot(predictions$finalModel)

plot(varImp(predictions), top = 10)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()

#Produce back-to-back histogram of vote distributions for Relapse and NoRelapse
options(digits=2) 
pdf(file=vote_dist_pdffile)
out <- histbackback(split(rf_output$votes[,"Relapse"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for patients classified by RF', axes=TRUE, ylab="Fraction votes (Relapse)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
dev.off()

#Save case predictions to file:
case_predictions=cbind(clindata,target,rf_output$predicted,rf_output$votes)
write.table(case_predictions,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)



