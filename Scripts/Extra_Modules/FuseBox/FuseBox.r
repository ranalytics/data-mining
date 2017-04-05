 Combining classifiers into a single ensemble model
Description

fuse provides a platform to take existing predictive models and combine there predictions into a single outcome.
Usage

fuse(mods, ...)

## Default S3 method:
fuse(mods,
             classes = NULL, 
             probs = TRUE, 
             predict = NULL, 
             weights = rep(1, length(mods)), 
             method = "vote", 
             methodArgs = NULL, 
             ...)

Arguments
mods 	

a named list of models
classes 	

a character string of possible classes
probs 	

a logical: will the model predict class probabilities (as opposed to the discrete class)
predict 	

an optional list the same length as mods that contains prediction functions for the models. By default, when probs = FALSE, samples are predicted using codepredict(model, newdata). When class probabilities are produced, the default syntax is predict(model, newdata, type = "prob"). The argument can be used for models that do not fit this convention or cases where the predictors do not use all the columns of newdata.
weights 	

a numeric vector the same length as mods of weights when averaging probabilities. These values will be normalized via weights/sum(weights).
method 	

usually a single method for combining the classifiers. Possible values are 'vote' (for majority vote), 'meanProb' (for weighted and unweighted averages of the class probabilities), 'prod' (the product of the class probabilities across models). Alternatively, a function with minimum arguments x and levels. See the Details section below.
methodArgs 	

an optional named list of arguments if a custom function is used with the method argument.
... 	

not currently used
Value

a list of class "fuse"
Author(s)

Max Kuhn
References

insert
See Also

predict.fuse, ~~~
Examples

sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
       if(trace) cat(nm,":")
       source(file.path(path, nm), ...)
       if(trace) cat("\n")
    }
 }
sourceDir("FuseBox/R")
## Not run: 
library(QSARdata)
data(Mutagen)

## Split the data three times into training and testing sets
library(caret)

set.seed(1)
inTrain <- createDataPartition(Mutagen_Outcome, p = 3/4, list = FALSE)

trainData <- Mutagen_Dragon[ inTrain,]
testData <- Mutagen_Dragon[-inTrain,]
trainClass <- Mutagen_Outcome[ inTrain]
testClass <- Mutagen_Outcome[-inTrain]

## There are some predictors with degenerate distirbutions, so 
## remvoe these

isNZV <- nearZeroVar(trainData)
trainData <- trainData[, -isNZV]
testData <- testData[, -isNZV]
calData <- calData[, -isNZV]

## Make a copy opf the data in a single data frame

training <- trainData
training$Class <- trainClass

dim(training)

## Fit a random forest model
library(randomForest)
rfModel <- randomForest(Class ~ ., data = training)

## Now an SVM model
library(kernlab)
svmModel <- ksvm(as.matrix(trainData), trainClass, C = 4, prob.model = TRUE)

## Create a list of the models and associated prediction functions:
models <- list(rf = rfModel,
               svm = svmModel)

pred <- list(function(x, dat, ...) predict(x, dat, type = "prob")[,1],
             function(x, dat, ...) predict(x, dat, type = "probabilities")[,1])

fusedMods <- fuse(list(rf = rfModel, svm = svmModel), 
                  probs = TRUE, 
                  predict = pred, 
                  method = "meanProb", 
                  classes = levels(testClass))
fuse()
confusionMatrix(predict(fusedMods, testData), testClass)
confusionMatrix(predict(rfModel, testData), testClass)
confusionMatrix(predict(svmModel, testData), testClass)

## End(Not run)
