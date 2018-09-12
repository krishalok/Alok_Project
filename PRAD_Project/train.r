library(randomForest)
library(mlbench)
library(caret)
library(data.table)
library(snow)
library(varSelRF)
library(forestFloor)
library(ROCR)
library(pROC)
library(caret)
library(rattle)
library(rpart.plot)
library(rpart)
library(limma)
#install.packages("randomForest")
#install.packages("ROCR")
#install.packages("Hmisc")
#source("http://bioconductor.org/biocLite.R")
#biocLite("genefilter")
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
library(snow)
library(varSelRF)

setwd("C:\\Users\\alokum\\Desktop\\PRAD")
datafile="PRAD_DATA.txt" 
clindatafile="trainset_clindetails.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GEO.asscession.number"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[2:length(colnames(data_import))])+1 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1,data_order)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata)
X=rawdata[,2:length(header)]
cvfun <- cv(14,14.1)
#cvfun <- cv(18,20)
ffun=filterfun(cvfun)
filt=genefilter(2^X,ffun)
filt_Data=rawdata[filt,] 
#Get potential predictor variables
predictor_data=as.data.frame(t(filt_Data[,2:length(header)])) #Filtered
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

#expression <- read.table(data_file, sep = "\t", header = T, row.names = 1)
#expression <- data.frame(t(expression))
#expression2 <- data.frame(predictor_data)

new_row_names <- c()
current_names <- row.names(predictor_data)
for (r in 1:length(current_names)){
  name <- current_names[r]
  new_row_names = c(new_row_names, paste("TCGA-", substr(name, 6,7), "-", 
                                         substr(name, 9, 12), sep = ""))
}
setDT(predictor_data, keep.rownames = TRUE)[]
predictor_data$rn<- new_row_names

predictor_data$rn = paste0(predictor_data$rn,'-01')
merged_data <- merge(clindata, predictor_data, by.x = "PID", by.y = "rn")
# Training  Data is ready here
setwd("C:\\Users\\alokum\\Desktop\\Training")

#Method-1 : specify that the resampling method is 
fit_control <- trainControl(## 10-fold CV
  method = "cv",
  number = 10)

rf_fit <- train(as.factor(relapse..1.True.) ~ ., 
                data = merged_data[,c(5,8:148)], 
                method="rf",
                trControl=fit_control,
                prox=TRUE,allowParallel=TRUE)
rf_fit
# Confusion Metrix 
print(rf_fit$finalModel)
plot(rf_fit)
# Method-2 : Random Search Method
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
#set.seed(seed)
mtry <- sqrt(ncol(merged_data))
rf_random <- train(as.factor(relapse..1.True.) ~ ., data=merged_data[,c(5,8:148)], method="rf", tuneLength=15, trControl=control)
print(rf_random)
plot(rf_random)
# Method-3 : Grid Search Method
metric <- "Accuracy"
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(.mtry=c(1:15))
rf_gridsearch <- train(as.factor(relapse..1.True.) ~ ., data=merged_data[,c(5,8:148)], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
print(rf_gridsearch)
plot(rf_gridsearch)


# Algorithm Tune (tuneRF)

set.seed(seed)
bestmtry <- tuneRF(merged_data, target, stepFactor=1.5, improve=1e-5, ntree=25001)
print(bestmtry)

# Manual Search
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(.mtry=c(sqrt(ncol(merged_data[,c(5,8:148)]))))
modellist <- list()
for (ntree in c(1001, 1501, 2001, 2501,25001)) {
  fit <- train(as.factor(relapse..1.True.) ~ ., data=merged_data[,c(5,8:148)], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}
# compare results
results <- resamples(modellist)
summary(results)
dotplot(results)


customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

# train model
control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(1:15), .ntree=c(1000, 1500, 2000, 2500))
#set.seed(seed)
custom <- train(as.factor(relapse..1.True.) ~ ., data=merged_data[,c(5,8:148)], method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
summary(custom)
plot(custom)

merged_data_gp<-merged_data
merged_data_gp$relapse..1.True. <- make.names(merged_data$relapse..1.True., unique = FALSE, allow_ = TRUE)

mod0 <- train(as.factor(relapse..1.True.) ~ ., data=merged_data_gp[,c(5,8:148)],
              method = "rf",
              metric = "ROC",
              tuneGrid = data.frame(mtry = 8),
              ntree = 1000,
              trControl = trainControl(method = "repeatedcv",
                                       repeats = 5,
                                       classProbs = TRUE,
                                       summaryFunction = twoClassSummary))


glmnetCM <- confusionMatrix(mod0, norm = "none")

library(pROC)

getTrainPerf(mod0)

coef(mod0$finalModel)



###plot decsion tree


## Testing the model_testing_GSE25136"

setwd("C:\\Users\\alokum\\Desktop\\PRAD\\testing_GSE25136")
data_file1="testset_gcrma_GSE25136_141.csv"
clin_datafile1="clinical_GSE251361.csv"
#data_file="PRAD_DATA.txt" 
#clindatafile="trainset_clindetails.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import1=read.table(data_file1, header = TRUE, na.strings = "NA", sep=",")
clin_data_import1=read.table(clin_datafile1, header = TRUE, na.strings = "NA", sep=",")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order1=order(clin_data_import1[,"geo_accn"])
clindata1=clin_data_import1[clin_data_order,]
data_order1=order(colnames(data_import1)[2:length(colnames(data_import1))])+1 #Order data without first three columns, then add 3 to get correct index in original file
rawdata1=data_import1[,c(1,data_order1)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata1)
X1=rawdata1[,2:length(header)]
#cvfun <- cv(14,14.1)
#cvfun <- cv(18,20)
#ffun=filterfun(cvfun)
#filt=genefilter(2^X,ffun)
#filt_Data=rawdata[filt,] 
filt_Data1<-rawdata1 
#Get potential predictor variables
predictor_data1=as.data.frame(t(filt_Data1[,2:length(header)])) #Filtered
predictor_names1=c(as.vector(filt_Data1[,1])) #gene symbol
colnames(predictor_data1)=predictor_names1
#Get target variable and specify as factor/categorical
target1= clindata1[,"event.rfs"]
target1[target1==0]="NoRelapse"
target1[target1==1]="Relapse"
target1=as.factor(target1)
#Run RandomForests
#NOTE: use an ODD number for ntree. When the forest/ensembl is used on test data, ties are broken randomly.
#Having an odd number of trees avoids this issue and makes the model fully deterministic
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target1))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

#new_row_names <- c()
#current_names <- row.names(predictor_data1)
#for (r in 1:length(current_names)){
#  name <- current_names[r]
#  new_row_names = c(new_row_names, paste("TCGA-", substr(name, 6,7), "-", 
#                                         substr(name, 9, 12), sep = ""))
#}
setDT(predictor_data1, keep.rownames = TRUE)[]
setnames(predictor_data1, "rn", "PID")
#predictor_data1$rn<- new_row_names

#predictor_data1$rn = paste0(predictor_data1$rn,'-01')
merged_data_test <- merge(clindata1, predictor_data1, by.x = "geo_accn", by.y = "PID")

#setDT(predictor_data1, keep.rownames = TRUE)[]
merged_data_test[is.na(merged_data_test)] <-0
#test_data<-subset(predictor_data1,clindata1,x.)
merged_data_test$event.rfs <- make.names(merged_data_test$event.rfs , unique = FALSE, allow_ = TRUE)

roc0 <- roc(as.factor(merged_data_test$event.rfs), 
            predict(mod0, merged_data_test[,4:144], type = "prob")[,1], 
            levels = rev(levels(as.factor(merged_data_test$event.rfs))))
roc0

plot(roc0, print.thres = c(.5), type = "S",
     print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)",
     print.thres.cex = .8, 
     legacy.axes = TRUE)

predictions_1<-  predict(mod0, merged_data_test[,4:144], type = "prob")
confusionMatrix(predictions_1, as.factor(merged_data_test$event.rfs))

## Testing the model_testing_GSE21032"
## Testing the model_testing_GSE21032"
setwd("C:\\Users\\alokum\\Desktop\\PRAD\\testing_GSE21032")
#setwd("C:\\Users\\alokum\\Desktop\\PRAD\\testing_GSE25136")
data_file2="testset_gcrma_MSKCC_141.csv"
clin_datafile2="testset_clindetails_MSKCCC.txt"
#data_file="PRAD_DATA.txt" 
#clindatafile="trainset_clindetails.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import2=read.table(data_file2, header = TRUE, na.strings = "NA", sep=",")
clin_data_import2=read.table(clin_datafile2, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order2=order(clin_data_import2[,"geo_accn"])
clindata2=clin_data_import2[clin_data_order2,]
data_order2=order(colnames(data_import2)[2:length(colnames(data_import2))])+1 #Order data without first three columns, then add 3 to get correct index in original file
rawdata2=data_import2[,c(1,data_order2)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata2)
X2=rawdata2[,2:length(header)]
#cvfun <- cv(14,14.1)
#cvfun <- cv(18,20)
#ffun=filterfun(cvfun)
#filt=genefilter(2^X,ffun)
#filt_Data=rawdata[filt,] 
filt_Data2<-rawdata2 
#Get potential predictor variables
predictor_data2=as.data.frame(t(filt_Data2[,2:length(header)])) #Filtered
predictor_names2=c(as.vector(filt_Data2[,1])) #gene symbol
colnames(predictor_data2)=predictor_names2
#Get target variable and specify as factor/categorical
target2= clindata2[,"event.rfs"]
target2[target2==0]="NoRelapse"
target2[target2==1]="Relapse"
target2=as.factor(target2)
#Run RandomForests
#NOTE: use an ODD number for ntree. When the forest/ensembl is used on test data, ties are broken randomly.
#Having an odd number of trees avoids this issue and makes the model fully deterministic
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target2))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

#new_row_names <- c()
#current_names <- row.names(predictor_data1)
#for (r in 1:length(current_names)){
#  name <- current_names[r]
#  new_row_names = c(new_row_names, paste("TCGA-", substr(name, 6,7), "-", 
#                                         substr(name, 9, 12), sep = ""))
#}
setDT(predictor_data2, keep.rownames = TRUE)[]
setnames(predictor_data2, "rn", "PID")
#predictor_data1$rn<- new_row_names

#predictor_data1$rn = paste0(predictor_data1$rn,'-01')
merged_data_test2 <- merge(clindata2, predictor_data2, by.x = "geo_accn", by.y = "PID")

#setDT(predictor_data1, keep.rownames = TRUE)[]
merged_data_test2[is.na(merged_data_test2)] <-0
#test_data<-subset(predictor_data1,clindata1,x.)
merged_data_test2$event.rfs <- make.names(merged_data_test2$event.rfs , unique = FALSE, allow_ = TRUE)

roc0_test2 <- roc(as.factor(merged_data_test2$event.rfs), 
            predict(mod0, merged_data_test2[,4:144], type = "prob")[,1], 
            levels = rev(levels(as.factor(merged_data_test2$event.rfs))))
roc0_test2

plot(roc0_test2, print.thres = c(.5), type = "S",
     print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)",
     print.thres.cex = .8, 
     legacy.axes = TRUE)

probs <- predict(mod0, merged_data_test2[,4:144], type = "prob") 


#predictions_1<-  predict(mod0, merged_data_test[,4:144], type = "prob")
confusionMatrix(probs, as.factor(merged_data_test2$event.rfs))
## Testing the model_testing_TCGA"
## Testing the model_testing_TCGA"
setwd("C:\\Users\\alokum\\Desktop\\PRAD\\\\testing_tcga")
data_file3="GS7_141_127.txt"
clin_datafile3="clin_127.txt"
#data_file="PRAD_DATA.txt" 
#clindatafile="trainset_clindetails.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import3=read.table(data_file3, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import3=read.table(clin_datafile3, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order3=order(clin_data_import3[,"geo_accn"])
clindata3=clin_data_import3[clin_data_order,]
data_order3=order(colnames(data_import3)[2:length(colnames(data_import3))])+1 #Order data without first three columns, then add 3 to get correct index in original file
rawdata3=data_import3[,c(1,data_order3)] #grab first three columns, and then remaining columns in order determined above
header=colnames(rawdata3)
X3=rawdata3[,2:length(header)]
#cvfun <- cv(14,14.1)
#cvfun <- cv(18,20)
#ffun=filterfun(cvfun)
#filt=genefilter(2^X,ffun)
#filt_Data=rawdata[filt,] 
filt_Data3<-rawdata3 
#Get potential predictor variables
predictor_data3=as.data.frame(t(filt_Data3[,2:length(header)])) #Filtered
predictor_names3=c(as.vector(filt_Data3[,1])) #gene symbol
colnames(predictor_data3)=predictor_names3
#Get target variable and specify as factor/categorical
target3= clindata3[,"event.rfs"]
target3[target3==0]="NoRelapse"
target3[target3==1]="Relapse"
target3=as.factor(target3)
#Run RandomForests
#NOTE: use an ODD number for ntree. When the forest/ensembl is used on test data, ties are broken randomly.
#Having an odd number of trees avoids this issue and makes the model fully deterministic
#Use down-sampling to attempt to compensate for unequal class-sizes
tmp = as.vector(table(target3))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)

new_row_names <- c()
current_names <- row.names(predictor_data3)
for (r in 1:length(current_names)){
  name <- current_names[r]
  new_row_names = c(new_row_names, paste("TCGA-", substr(name, 6,7), "-", 
                                         substr(name, 9, 12), sep = ""))
}
setDT(predictor_data3, keep.rownames = TRUE)[]
setnames(predictor_data3, "rn", "PID")
#predictor_data1$rn<- new_row_names
predictor_data3$rn<- new_row_names

predictor_data3$PID = paste0(predictor_data3$PID,'-01')
#predictor_data1$rn = paste0(predictor_data1$rn,'-01')
merged_data_test3 <- merge(clindata3, predictor_data3, by.x = "geo_accn", by.y = "PID")

#setDT(predictor_data1, keep.rownames = TRUE)[]
merged_data_test3[is.na(merged_data_test3)] <-0
#test_data<-subset(predictor_data1,clindata1,x.)
merged_data_test3$event.rfs <- make.names(merged_data_test3$event.rfs , unique = FALSE, allow_ = TRUE)

roc0_test3 <- roc(as.factor(merged_data_test3$event.rfs), 
                  predict(mod0, merged_data_test3[,14:154], type = "prob")[,1], 
                  levels = rev(levels(as.factor(merged_data_test3$event.rfs))))
roc0_test3
str(roc0_test3)
plot(roc0_test3, print.thres = c(.5), type = "S",
     print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)",
     print.thres.cex = .8, 
     legacy.axes = TRUE)

plot(roc0, col = 1, lty = 2, main = "ROC")
plot(roc0_test2, col = 4, lty = 3, add = TRUE)
plot(roc0_test3, col = 4, lty = 3, add = TRUE)



probs <- predict(mod0, merged_data_test3[,14:154], type = "prob") 
################Plotting of data################

vi <- varImp(mod0)
plot(mod0$finalModel)
plot(varImp(mod0), top = 40)

summary(mod0)
roc0



################Logistic Regression################
################Plotting of data################
predict <- predict(mod0, type = 'prob')[,2]
table(merged_data_gp$relapse..1.True.)
ROCRpred <- prediction(predict, merged_data_gp$relapse..1.True.)
ROCRperf <- performance(ROCRpred, 'tpr','fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2,1.7))

################Plotting of data################

mod0 <- train(as.factor(relapse..1.True.) ~ ., data=merged_data_gp[,c(5,8:148)],
              method = "rf",
              metric = "ROC",
              tuneGrid = data.frame(mtry = 8),
              ntree = 1000,
              trControl = trainControl(method = "repeatedcv",
                                       repeats = 5,
                                       classProbs = TRUE,
                                       summaryFunction = twoClassSummary))

#logistic regression model
model <- glm (as.factor(relapse..1.True.) ~ ., data = merged_data_gp[,c(5,8:148)], family = binomial)
summary(model)
predict <- predict(model, type = 'response')
table(merged_data_gp$relapse..1.True., predict > 0.5)


library(ROCR)
ROCRpred <- prediction(predict, merged_data_gp$relapse..1.True.)
ROCRperf <- performance(ROCRpred, 'tpr','fpr')
plot(ROCRperf, colorize = TRUE, text.adj = c(-0.2,1.7))

predTst <- predict(model, merged_data_test2[,4:144], type="response")


#plot glm
library(ggplot2)
ggplot(merged_data_gp[,c(5,8:148)], aes(x=Rating, y=Recommended)) + geom_point() + 
  stat_smooth(method="glm", family="binomial", se=FALSE)




###########################Finding The best Model#######################


# load the library
library(mlbench)
library(caret)
# load the dataset
data(PimaIndiansDiabetes)
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the LVQ model
set.seed(7)
modelLvq <- train(as.factor(relapse..1.True.) ~ ., data = merged_data_gp[,c(5,8:148)], method="lvq", trControl=control)
# train the GBM model
set.seed(7)
modelGbm <- train(as.factor(relapse..1.True.) ~ ., data = merged_data_gp[,c(5,8:148)], trControl=control, verbose=FALSE)
# train the SVM model
set.seed(7)
modelSvm <- train(as.factor(relapse..1.True.) ~ ., data = merged_data_gp[,c(5,8:148)], trControl=control)
# collect resamples
results <- resamples(list(LVQ=modelLvq, GBM=modelGbm, SVM=modelSvm))
# summarize the distributions
summary(results)
# boxplots of results
bwplot(results)
# dot plots of results
dotplot(results)

# load the library
library(mlbench)
library(caret)
# load the dataset
data(PimaIndiansDiabetes)
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the LVQ model
set.seed(7)
modelLvq <- train(diabetes~., data=PimaIndiansDiabetes, method="lvq", trControl=control)
# train the GBM model
set.seed(7)
modelGbm <- train(diabetes~., data=PimaIndiansDiabetes, method="gbm", trControl=control, verbose=FALSE)
# train the SVM model
set.seed(7)
modelSvm <- train(diabetes~., data=PimaIndiansDiabetes, method="svmRadial", trControl=control)
# collect resamples
results <- resamples(list(LVQ=modelLvq, GBM=modelGbm, SVM=modelSvm))
# summarize the distributions
summary(results)
# boxplots of results
bwplot(results)
# dot plots of results
dotplot(results)



##########Condusion Matrix#######

#mod0 <- train(as.factor(relapse..1.True.) ~ ., data=merged_data_gp[,c(5,8:148)],
mod0 <- train(as.factor(relapse..1.True.) ~ ., data=merged_data_gp_S,
              method = "rf",
              metric = "ROC",
              tuneGrid = data.frame(mtry = 8),
              ntree = 1001,
              trControl = trainControl(method = "repeatedcv",
                                       repeats = 5,
                                       classProbs = TRUE,
                                       summaryFunction = twoClassSummary,savePredictions = T))

selectedIndices <- mod0$pred$mtry == 8

plot(roc(mod0$pred$obs[selectedIndices],mod0$pred$X0[selectedIndices]))

plot.roc(mod0$pred$obs[selectedIndices],mod0$pred$X0[selectedIndices])
varImp(mod0)$importance
merged_data_test_N<- as.data.frame(scale(merged_data_test[,4:144]))



rocN0 <- roc(as.factor(merged_data_test$event.rfs), 
            predict(mod0, merged_data_test_N, type = "prob")[,1], 
            levels = rev(levels(as.factor(merged_data_test$event.rfs))))
#colnames(merged_data_test)[which(names(merged_data_test) == "time.rfs")] <- "time.to.relapse.or.last.follow.up..months."
rocN0

merged_data_test_N2<- as.data.frame(scale(merged_data_test2[,4:144]))



rocN1 <- roc(as.factor(merged_data_test2$event.rfs), 
            predict(mod0, merged_data_test_N2, type = "prob")[,1], 
            levels = rev(levels(as.factor(merged_data_test2$event.rfs))))
#colnames(merged_data_test)[which(names(merged_data_test) == "time.rfs")] <- "time.to.relapse.or.last.follow.up..months."
rocN1

merged_data_test_N2<- as.data.frame(scale(merged_data_test3[,14:154]))


rocN2 <- roc(as.factor(merged_data_test3$event.rfs), 
             predict(mod0, merged_data_test3[,14:154], type = "prob")[,1], 
             levels = rev(levels(as.factor(merged_data_test3$event.rfs))))
#colnames(merged_data_test)[which(names(merged_data_test) == "time.rfs")] <- "time.to.relapse.or.last.follow.up..months."
rocN2

modelroc<-roc(mod0$pred$obs[selectedIndices],mod0$pred$X0[selectedIndices])
plot(modelroc, print.auc=TRUE,location=c(0.2,1),auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
plot(modelroc,
     print.auc = TRUE,
     auc.polygon = TRUE,
     grid=c(0.1, 0.2),
     grid.col = c("green", "red"),
     max.auc.polygon = TRUE,
     auc.polygon.col = "skyblue",
     print.thres = TRUE,
     print.auc.x = 0.3,
     print.auc.y = 0.2)
plot(roc(mod0$pred$obs[selectedIndices],mod0$pred$X0[selectedIndices]))
plot(rocN0, col = 1, lty = 2, main = "ROC")
plot(rocN1, col = 4, lty = 3, add = TRUE)
plot(rocN2, col = 4, lty = 3, add = TRUE)


library(dplyr)
col_index <- varImp(mod0)$importance %>% 
  mutate(names=row.names(.)) %>%
  arrange(-Overall)

imp_names <- col_index$names[1:31]


merged_data_gp_S<-merged_data_gp[,imp_names]


merged_data_gp_S<-merged_data_gp[,c("relapse..1.True.","FOXS1","KIF18B","MS4A7","ASPM","SLC41A2","ARHGAP11A","HOXB9","DUSP8","SLC15A3",  
                                  "CD86","ALG9","ARL11","CMPK2","TWIST2" ,"THBS3","MAL","SLC16A14","RACGAP1P", 
                                  "TEX101","ACER1","FAT4","EFCAB10","ICAM4","JRK" ,"AP3M2","CCDC39","SLC35F4",  
                                  "RIPPLY2","FBXW8","KLHL30","HOXB4")]

merged_data_test_S1<-merged_data_test[,c("event.rfs","FOXS1","KIF18B","MS4A7","ASPM","SLC41A2","ARHGAP11A","HOXB9","DUSP8","SLC15A3",  
                                    "CD86","ALG9","ARL11","CMPK2","TWIST2" ,"THBS3","MAL","SLC16A14","RACGAP1P", 
                                    "TEX101","ACER1","FAT4","EFCAB10","ICAM4","JRK" ,"AP3M2","CCDC39","SLC35F4",  
                                    "RIPPLY2","FBXW8","KLHL30","HOXB4")]

merged_data_test_N<- as.data.frame(scale(merged_data_test_S1[,2:32]))
rocN1 <- roc(as.factor(merged_data_test_S1$event.rfs), 
             predict(mod0, merged_data_test_N, type = "prob")[,1], 
             levels = rev(levels(as.factor(merged_data_test_S1$event.rfs))))
#colnames(merged_data_test)[which(names(merged_data_test) == "time.rfs")] <- "time.to.relapse.or.last.follow.up..months."
rocN1


merged_data_test_S2<-merged_data_test2[,c("event.rfs","FOXS1","KIF18B","MS4A7","ASPM","SLC41A2","ARHGAP11A","HOXB9","DUSP8","SLC15A3",  
                                         "CD86","ALG9","ARL11","CMPK2","TWIST2" ,"THBS3","MAL","SLC16A14","RACGAP1P", 
                                         "TEX101","ACER1","FAT4","EFCAB10","ICAM4","JRK" ,"AP3M2","CCDC39","SLC35F4",  
                                         "RIPPLY2","FBXW8","KLHL30","HOXB4")]
merged_data_test_N2<- as.data.frame(scale(merged_data_test_S2[,2:32]))
rocN2 <- roc(as.factor(merged_data_test_S2$event.rfs), 
             predict(mod0, merged_data_test_N2, type = "prob")[,1], 
             levels = rev(levels(as.factor(merged_data_test_S2$event.rfs))))
#colnames(merged_data_test)[which(names(merged_data_test) == "time.rfs")] <- "time.to.relapse.or.last.follow.up..months."
rocN2


merged_data_test_S3<-merged_data_test3[,c("event.rfs","FOXS1","KIF18B","MS4A7","ASPM","SLC41A2","ARHGAP11A","HOXB9","DUSP8","SLC15A3",  
                                         "CD86","ALG9","ARL11","CMPK2","TWIST2" ,"THBS3","MAL","SLC16A14","RACGAP1P", 
                                         "TEX101","ACER1","FAT4","EFCAB10","ICAM4","JRK" ,"AP3M2","CCDC39","SLC35F4",  
                                         "RIPPLY2","FBXW8","KLHL30","HOXB4")]

rocN3 <- roc(as.factor(merged_data_test_S3$event.rfs), 
             predict(mod0, merged_data_test_S3[2:32], type = "prob")[,1], 
             levels = rev(levels(as.factor(merged_data_test_S3$event.rfs))))
#colnames(merged_data_test)[which(names(merged_data_test) == "time.rfs")] <- "time.to.relapse.or.last.follow.up..months."
rocN3

##############################################################################################
##################################plotting ROC############################################################
##############################################################################################
plot.ROC = function(p,cutoff=0.5,filename=NA,...) {
  require(ROCR)
  require(RJSONIO)
  require(ggplot2)
  require(gridSVG)
  require(grid)
  
  # linear interpolation of x 
  interp = function(x,xvalues,yvalues) {
    if (x <= min(xvalues)) {
      return(yvalues[length(yvalues)])
    }
    if (x >= max(xvalues)) {
      return(yvalues[1])
    }
    i = 1
    while(x < xvalues[i]) {
      i = i + 1
    }
    return(yvalues[i-1]+ (yvalues[i]-yvalues[i-1]) *  
             (xvalues[i-1]-x) / (xvalues[i-1]-xvalues[i]))
  }
  
  auc.color.fill = rgb(216,216,255,128,maxColorValue=255)
  auc.color.text = rgb(0,0,0.6)
  bar.color.neg = rgb(179,179,179,179,maxColorValue=255)
  bar.color.pos = rgb(0,0,0.6,0.7)
  
  # grab some performance metrics
  auc = performance(p,"auc")@y.values[[1]]
  acc = performance(p,"acc")
  roc = performance(p,"tpr","fpr")
  roc@alpha.values[[1]][1]=1
  roc.df = data.frame(fpr=roc@x.values[[1]],tpr=roc@y.values[[1]],cutoff=roc@alpha.values[[1]])
  
  # draw the plots
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(3, 3)))
  v1 = viewport(layout.pos.col=c(1,2), layout.pos.row=c(1,2))
  p1= ggplot(data=roc.df) +
    xlab(roc@x.name) + ylab(roc@y.name) + ggtitle("ROC Curve") +
    geom_ribbon(aes(x=fpr,ymin=0,ymax=tpr),fill=auc.color.fill) +
    geom_line(aes(x=fpr,y=tpr),color=rgb(0,0,1,maxColorValue=255)) +
    annotate("segment",x=0,xend=1,y=0,yend=1,linetype=2,color="grey50") +
    annotate("text",x=0.96,y=0.05,size=4,hjust=1,label=paste("AUC =",format(auc,digits=2)),col=auc.color.text) +
    annotate("point",x=interp(cutoff,roc@alpha.values[[1]],roc@x.values[[1]]),y=interp(cutoff,roc@alpha.values[[1]],roc@y.values[[1]]),colour="red",size=3) +
    coord_cartesian(xlim = c(0, 1),ylim=c(0,1))
  
  v2 = viewport(layout.pos.col=c(3), layout.pos.row=c(1))
  baseline = max(c(p@n.pos[[1]],p@n.neg[[1]]))/(p@n.pos[[1]] + p@n.neg[[1]])
  p2 = ggplot(aes(xmin=0,ymin=0,xmax=1,ymax=1),data=data.frame(x=acc@x.values[[1]],y=acc@y.values[[1]])) +
    xlab(acc@x.name) + ylab(acc@y.name) +
    geom_line(aes(x=x,y=y),color=rgb(0,0,2,maxColorValue=255)) +
    annotate("segment",x=0,xend=1,y=baseline,yend=baseline,linetype=2,color=auc.color.text) +
    annotate("text",x=0.5,y=baseline,vjust=1.5,size=4,label="baseline",color=auc.color.text) +
    geom_vline(xintercept=cutoff,color="red",size=0.4) +
    coord_cartesian(xlim = c(0, 1))
  p2 = p2 + theme_classic()
  
  v3 = viewport(layout.pos.col=c(3), layout.pos.row=c(2))
  p3 = ggplot(aes(xmin=0,ymin=0,xmax=1,ymax=1),data=data.frame(x=roc@alpha.values[[1]],y=roc@y.values[[1]])) +
    xlab(roc@alpha.name) + ylab(roc@y.name) + 
    geom_line(aes(x=x,y=y),color=rgb(0,0,3,maxColorValue=255)) +
    geom_vline(xintercept=cutoff,color="red",size=0.4) +
    coord_cartesian(xlim = c(0, 1))
  p3 = p3 + theme_classic()
  
  v4 = viewport(layout.pos.col=c(3), layout.pos.row=c(3))
  p4 = ggplot(aes(xmin=0,ymin=0,xmax=1,ymax=1),data=data.frame(x=roc@alpha.values[[1]],y=roc@x.values[[1]])) +
    xlab(roc@alpha.name) + ylab(roc@x.name) + 
    geom_line(aes(x=x,y=y),color=rgb(0,0,4,maxColorValue=255)) +
    geom_vline(xintercept=cutoff,color="red",size=0.4) +
    coord_cartesian(xlim = c(0, 1))
  p4 = p4 + theme_classic()
  
  sp = p@predictions[[1]]
  sp = sort(sp,decreasing=TRUE)
  v5 = viewport(layout.pos.col=c(1,2), layout.pos.row=c(3))
  p5 = ggplot(aes(ymin=0,ymax=1),data=data.frame(x=1:length(sp),y=sp,col=p@labels[[1]])) +
    xlab(paste("Predictions  (N = ",length(sp),")",sep="")) +
    ylab("Probability") + 
    geom_segment(aes(x=x,y=0,xend=x,yend=y,colour=col),size=max(c(100/length(sp),0.15))) +
    scale_colour_manual(breaks=c("1","0"),labels=c("Positive Outcome","Negative Outcome"),
                        values = c(bar.color.neg,bar.color.pos)) +
    guides(colour=guide_legend(title=NULL)) +
    geom_hline(yintercept=cutoff,color="red",size=0.4) +
    coord_cartesian(ylim = c(0, 1)) + 
    theme(legend.position="top",legend.background=element_rect(fill=rgb(1,1,1,0.2)),
          legend.key.size=unit(0.7,"char"),
          legend.direction="horizontal",
          plot.margin = unit(c(0,4,2,3),"mm"),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  print(p1,vp=v1)
  print(p2,vp=v2)
  print(p3,vp=v3)
  print(p4,vp=v4)
  print(p5,vp=v5)
  
  if (!is.na(filename)) {
    # add animation
    jroc.df = toJSON(roc.df)
    grid.script(script=paste("roc =", jroc.df),inline=TRUE)
    code = "
  function interp(t,e,n){if(e[0]<e[e.length-1])return interp(t,e.reverse(),n.reverse());if(t<=e[e.length-1])return n[n.length-1];if(t>=e[0])return n[0];for(i=0;t<e[i];)i+=1;return n[i-1]+(n[i]-n[i-1])*(e[i-1]-t)/(e[i-1]-e[i])}function animAttr(t,e,i){t.setAttribute('to',i),t.setAttribute('from',e),t.beginElement()}function animLine(t,e,i,n){'x'==n?(x=viewportConvertPos(grobViewport(t),i,1,from='npc',to='svg').x,y1=viewportConvertPos(grobViewport(t),0,0,from='npc',to='svg').y,y2=viewportConvertPos(grobViewport(t),1,1,from='npc',to='svg').y,newpoints=x+','+y1+' '+x+','+y2,animAttr(e,document.getElementById(t).getAttribute('points'),newpoints)):(y=viewportConvertPos(grobViewport(t),1,i,from='npc',to='svg').y,x1=viewportConvertPos(grobViewport(t),0,0,from='npc',to='svg').x,x2=viewportConvertPos(grobViewport(t),1,1,from='npc',to='svg').x,newpoints=x1+','+y+' '+x2+','+y,animAttr(e,document.getElementById(t).getAttribute('points'),newpoints))}function moveCutoff(t){var e=document.getElementById(dot_id),i=document.getElementById('line1a'),n=document.getElementById('line2a'),r=document.getElementById('line3a'),o=document.getElementById('line4a'),l=document.getElementById('dota'),d=document.getElementById('dotb');animLine(line1_id,i,t,'x'),animLine(line2_id,n,t,'x'),animLine(line3_id,r,t,'x'),animLine(line4_id,o,t,'y');var u={x:interp(t,roc.cutoff,roc.fpr),y:interp(t,roc.cutoff,roc.tpr)},m=viewportConvertPos(grobViewport(dot_id),u.x,u.y,from='npc',to='svg');animAttr(l,e.getAttribute('x'),m.x),animAttr(d,e.getAttribute('y'),m.y)}function toSVG(t){var e=t.target.ownerSVGElement.createSVGPoint();e.x=t.clientX,e.y=t.clientY;var i=t.target.getScreenCTM(),n=e.matrixTransform(i.inverse());return n}function click2(t){var e=toSVG(t);moveCutoff(viewportConvertPos(grobViewport(t.target.id),e.x,e.y,from='svg',to='npc').x)}function click3(t){var e=toSVG(t);moveCutoff(viewportConvertPos(grobViewport(t.target.id),e.x,e.y,from='svg',to='npc').y)}function click(t){var e=toSVG(t),i=viewportConvertPos(grobViewport(t.target.id),e.x,e.y,from='svg',to='npc');moveCutoff(interp(i.x,roc.fpr,roc.cutoff))}function showNode(t){if(1==t.nodeType){if(t.hasAttribute('fill')){var e=t.getAttribute('fill');if(e&&'rgb(255,0,0)'==e){dot_id=t.getAttribute('id');var i=t.parentElement.parentElement;i.onclick=click,i.setAttribute('cursor','crosshair'),i=i.parentElement;var n=i.getAttribute('clip-path');n=n.slice(5,-1),n=document.getElementById(n),n=n.firstElementChild;var r=n.getAttribute('x');r=Number(r)-5,n.setAttribute('x',r),r=n.getAttribute('y'),r=Number(r)-5,n.setAttribute('y',r),r=n.getAttribute('width'),r=Number(r)+10,n.setAttribute('width',r),r=n.getAttribute('height'),r=Number(r)+10,n.setAttribute('height',r)}else e&&'rgb(216,216,255)'==e&&(dline5_id=t.getAttribute('id'))}if(t.hasAttribute('stroke')){if(e=t.getAttribute('stroke'),e&&'rgb(255,0,0)'==e)if(''==line1_id){line1_id=t.getAttribute('id');var i=t.parentElement.parentElement.parentElement;i.onclick=click2,i.setAttribute('cursor','crosshair')}else if(''==line2_id){line2_id=t.getAttribute('id');var i=t.parentElement.parentElement.parentElement;i.onclick=click2,i.setAttribute('cursor','crosshair')}else if(''==line3_id){line3_id=t.getAttribute('id');var i=t.parentElement.parentElement.parentElement;i.onclick=click2,i.setAttribute('cursor','crosshair')}else if(''==line4_id){line4_id=t.getAttribute('id');var i=t.parentElement.parentElement.parentElement;i.setAttribute('cursor','crosshair'),i.onclick=click3}e&&'rgb(0,0,1)'==e&&(dline1_id=t.getAttribute('id')),e&&'rgb(0,0,2)'==e&&(dline2_id=t.getAttribute('id')),e&&'rgb(0,0,3)'==e&&(dline3_id=t.getAttribute('id')),e&&'rgb(0,0,4)'==e&&(dline4_id=t.getAttribute('id')),e&&'rgb(179,179,179)'==e&&''==dline6_id&&(dline6_id=t.parentElement.parentElement.getAttribute('id'))}}}function walkDOM(t,e){for(e(t),t=t.firstChild;t;)walkDOM(t,e),t=t.nextSibling}function makeAnim(t,e,i)
{var n=document.createElementNS('http://www.w3.org/2000/svg','animate');n.setAttributeNS(null,'attributeName',i),n.setAttributeNS(null,'fill','freeze'),n.setAttributeNS(null,'begin','indefinite'),n.setAttribute('dur','200ms'),n.setAttributeNS('http://www.w3.org/1999/xlink','href','#'+e),n.setAttributeNS(null,'id',t),document.documentElement.appendChild(n)}function initDOM(){walkDOM(document.documentElement,showNode),makeAnim('dota',dot_id,'x'),makeAnim('dotb',dot_id,'y'),makeAnim('line1a',line1_id,'points'),makeAnim('line2a',line2_id,'points'),makeAnim('line3a',line3_id,'points'),makeAnim('line4a',line4_id,'points'),initAnim(),document.documentElement.setAttribute('width','100%'),document.documentElement.setAttribute('height','100%'),document.documentElement.setAttribute('preserveAspectRatio','xMinYMid meet')}function animPolyline(t,e,i,n){var r=document.createElementNS('http://www.w3.org/2000/svg','animate');r.setAttributeNS(null,'attributeName',n),r.setAttributeNS(null,'fill','freeze'),r.setAttributeNS(null,'begin','400ms'),r.setAttribute('dur','500ms'),r.setAttributeNS(null,'from',e),r.setAttributeNS(null,'to',i),r.setAttributeNS('http://www.w3.org/1999/xlink','href','#'+t),document.documentElement.appendChild(r)}function initAnim(){var t=dline1_id,e=document.getElementById(t);p1=viewportConvertPos(grobViewport(t),0,0,from='npc',to='svg'),p2=viewportConvertPos(grobViewport(t),1,1,from='npc',to='svg');var i,n='',r=roc.fpr.length,o=(p2.x-p1.x)/r,l=(p2.y-p1.y)/r;for(i=0;r>i;i++)n+=(p1.x+i*o).toFixed(2)+','+(p1.y+i*l).toFixed(2)+' ';var d=e.getAttribute('points');for(e.setAttribute('points',n),animPolyline(t,n,d,'points'),t=dline5_id,e=document.getElementById(t),i=r-1;i>=0;i--)n+=(p1.x+i*o).toFixed(2)+','+p1.y.toFixed(2)+' ';for(d=e.getAttribute('points'),e.setAttribute('points',n),animPolyline(t,n,d,'points'),t=dline2_id,p1=viewportConvertPos(grobViewport(t),0,0,from='npc',to='svg'),p2=viewportConvertPos(grobViewport(t),1,1,from='npc',to='svg'),o=(p2.x-p1.x)/(r/2),l=(p2.y-p1.y)/(r/2),e=document.getElementById(t),n='',i=0;r>i;i++)n+=r/2>i?p1.x.toFixed(2)+','+(p2.y-i*l).toFixed(2)+' ':(p1.x+i*l).toFixed(2)+','+p1.y.toFixed(2)+' ';for(d=e.getAttribute('points'),e.setAttribute('points',n),animPolyline(t,n,d,'points'),t=dline3_id,p1=viewportConvertPos(grobViewport(t),0,0,from='npc',to='svg'),p2=viewportConvertPos(grobViewport(t),1,1,from='npc',to='svg'),o=(p2.x-p1.x)/(r/2),l=(p2.y-p1.y)/(r/2),e=document.getElementById(t),n='',i=0;r>i;i++)n+=r/2>i?p1.x.toFixed(2)+','+(p2.y-i*l).toFixed(2)+' ':(p1.x+i*l).toFixed(2)+','+p1.y.toFixed(2)+' ';for(d=e.getAttribute('points'),e.setAttribute('points',n),animPolyline(t,n,d,'points'),t=dline4_id,p1=viewportConvertPos(grobViewport(t),0,0,from='npc',to='svg'),p2=viewportConvertPos(grobViewport(t),1,1,from='npc',to='svg'),o=(p2.x-p1.x)/(r/2),l=(p2.y-p1.y)/(r/2),e=document.getElementById(t),n='',i=0;r>i;i++)n+=r/2>i?p1.x.toFixed(2)+','+(p2.y-i*l).toFixed(2)+' ':(p1.x+i*l).toFixed(2)+','+p1.y.toFixed(2)+' ';d=e.getAttribute('points'),e.setAttribute('points',n),animPolyline(t,n,d,'points'),t=dline6_id,e=document.getElementById(t);var u=document.createElementNS('http://www.w3.org/2000/svg','animateTransform');u.setAttributeNS(null,'attributeName','transform'),u.setAttributeNS(null,'type','scale'),u.setAttributeNS(null,'fill','freeze'),u.setAttributeNS(null,'begin','400ms'),u.setAttribute('dur','500ms'),u.setAttributeNS(null,'from','1 ,0'),u.setAttributeNS(null,'to','1 ,1'),u.setAttributeNS('http://www.w3.org/1999/xlink','href','#'+t),e.setAttributeNS(null,'transform','scale(1 0)'),document.documentElement.appendChild(u);var u=document.createElementNS('http://www.w3.org/2000/svg','animate');u.setAttributeNS(null,'attributeName','opacity'),u.setAttributeNS(null,'fill','freeze'),u.setAttributeNS(null,'begin','0ms'),u.setAttribute('dur','400ms'),u.setAttributeNS(null,'from',0),u.setAttributeNS(null,'to',1),document.documentElement.setAttributeNS(null,'opacity',0),document.documentElement.appendChild(u)}var dot_id='',line1_id='',line2_id='',line3_id='',line4_id='',dline1_id='',dline2_id='',dline3_id='',dline4_id='',dline5_id='',dline6_id='';document.documentElement.onload=initDOM;
  
    "
    grid.script(script=code,inline=TRUE)
    grid.export(name=filename,exportJS="inline",exportCoords="inline",progress=FALSE)
  }
}

p = predict(mod0,newdata=merged_data_test_S1[2:32],type="prob")[,2]

png("ROCPlot.png",width=8,height=8,units="in",res=72)
plot.ROC(prediction(p,merged_data_test_S1$event.rfs),0.25,"ROCplot.svg")
dev.off()



p = predict(mod0,newdata=merged_data_test_S2,type="prob")[,2]

png("ROCPlot_2.png",width=8,height=8,units="in",res=72)
plot.ROC(prediction(p,merged_data_test_S2$event.rfs),0.25,"ROCplot2.svg")
dev.off()

p = predict(mod0,newdata=merged_data_test_S3,type="prob")[,2]

png("ROCPlot3.png",width=8,height=8,units="in",res=72)
plot.ROC(prediction(p,merged_data_test_S3$event.rfs),0.25,"ROCplot3.svg")
dev.off()





