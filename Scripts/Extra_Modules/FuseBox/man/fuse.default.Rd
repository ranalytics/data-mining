\name{fuse.default}
\alias{fuse.default}
\alias{fuse}
\title{
Combining classifiers into a single ensemble model
}
\description{
\code{fuse} provides a platform to take existing predictive models and combine there predictions into a single outcome.
}
\usage{
fuse(mods, ...)

\method{fuse}{default}(mods,
             classes = NULL, 
             probs = TRUE, 
             predict = NULL, 
             weights = rep(1, length(mods)), 
             method = "vote", 
             methodArgs = NULL, 
             ...)
}
\arguments{
  \item{mods}{
a named list of models
}
  \item{classes}{
a character string of possible classes
}
  \item{probs}{
a logical: will the model predict class probabilities (as opposed to the discrete class) 
}
  \item{predict}{
an optional list the same length as \code{mods} that contains prediction functions for the models. By default, when \code{probs = FALSE}, samples are predicted using code{predict(model, newdata)}. When class probabilities are produced, the default syntax is \code{predict(model, newdata, type = "prob")}. The argument can be used for models that do not fit this convention or cases where the predictors do not use all the columns of \code{newdata}.
 }
  \item{weights}{
a numeric vector the same length as \code{mods} of weights when averaging probabilities. These values will be normalized via \code{weights/sum(weights)}.
}
  \item{method}{
usually a single method for combining the classifiers. Possible values are 'vote' (for majority vote), 'meanProb' (for weighted and unweighted averages of the class probabilities), 'prod' (the product of the class probabilities across models). Alternatively, a function with minimum arguments \code{x} and \code{levels}. See the Details section below. 
}
  \item{methodArgs}{
an optional named list of arguments if a custom function is used with the \code{method} argument.
}
  \item{\dots}{
not currently used
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list of class "fuse"
}
\references{
insert
}
\author{
Max Kuhn
}

\seealso{
 \code{\link{predict.fuse}}, ~~~
}
\examples{
\dontrun{
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

confusionMatrix(predict(fusedMods, testData), testClass)
confusionMatrix(predict(rfModel, testData), testClass)
confusionMatrix(predict(svmModel, testData), testClass)
}

}
\keyword{models }