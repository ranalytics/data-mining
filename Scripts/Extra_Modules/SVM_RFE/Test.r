# testing the SVM-RFE algorithm

library(e1071)

#set working directory
setwd("C:\\Hernan\\TestFolder\\SVM-RFE")
#load database
db <- read.csv(file="databases\\colon-cancer\\colon-cancer.csv",head=FALSE,sep=",")
#db <- read.csv(file="databases\\leukemia\\leukemia.csv",head=FALSE,sep=",")

x = as.matrix(db[,1:(ncol(db)-1)])
y = as.factor(db[,ncol(db)])

x = log(x)
#scale samples with mean zero and standard deviation one
for(i in 1:nrow(x))x[i,] = (x[i,]-mean(x[i,]))/sd(x[i,])
#scale features with mean zero and standard deviation one
for(i in 1:ncol(x))x[,i] = (x[,i]-mean(x[,i]))/sd(x[,i])
x = 2*atan(x/2)


# Feature Ranking with SVM-RFE
featureRankedList = svmrfeFeatureRanking(x,y)


# Leave One Out test using the best n ranked features of the SVM-RFE
for(pow in 1:10){
  nfeatures = 2^pow
  truePredictions = 0
  for(i in 1:nrow(x)){
    svmModel = svm(x[-i, featureRankedList[1:nfeatures]], y[-i], cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" ) 
    prediction = predict(svmModel,matrix(x[i, featureRankedList[1:nfeatures]],nrow=1))
    if(prediction[[1]] == y[[i]]) truePredictions = truePredictions + 1
  }
  cat(nfeatures,":",truePredictions/nrow(x),"\n")
}
