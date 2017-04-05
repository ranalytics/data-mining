#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 6. БИНАРНЫЕ КЛАССИФИКАТОРЫ С РАЗЛИЧНЫМИ РАДЕЛЯЮЩИМИ ПОВЕРХНОСТЯМИ
#########################################################################

#  6.1. Дискриминантный анализ
#-----------------------------------------------------------------------

DGlass <- read.table(file="Glass.txt", sep=","
,header=TRUE,row.names=1)
DGlass$F <- as.factor(ifelse(DGlass$Class==2,2,1))
library(mvnormtest) 
mshapiro.test(t(DGlass[DGlass$F == 1, 1:9])) 
mshapiro.test(t(DGlass[DGlass$F == 2, 1:9]))
library(biotools) 
boxM(as.matrix(DGlass[, 1:9]), DGlass$F ) 

# Функция вывода результатов классифицирования
Out_CTab <- function (model, group, type="lda") {
# Таблица неточности "Факт/Прогноз" по обучающей выборке
classified<-predict(model)$class  
t1 <-  table(group,classified)  
# Точность классифицирования  и расстояние Махалонобиса
# Функция вывода результатов классифицирования
Out_CTab <- function (model, group, type="lda") {
# Таблица неточности "Факт/Прогноз" по обучающей выборке
   classified<-predict(model)$class  
   t1 <-  table(group,classified)  
# Точность классифицирования  и расстояние Махалонобиса
   Err_S <- mean(group != classified) ; mahDist <- NA
   if (type=="lda") 
       { mahDist <- dist(model$means %*% model$scaling) }
# Таблица "Факт/Прогноз" и ошибка при скользящем контроле
   t2 <-  table(group, update(model, CV=T)$class->LDA.cv) 
   Err_CV <- mean(group != LDA.cv) 
   Err_S.MahD  <-c(Err_S, mahDist) 
   Err_CV.N <-c(Err_CV, length(group)) 
   cbind(t1, Err_S.MahD, t2, Err_CV.N)
}
# --- Выполнение расчетов
library(MASS)
lda.all <- lda(F ~ ., data=DGlass[,-10])
Out_CTab(lda.all, DGlass$F) 

library ( klaR)
stepclass(F ~ ., data=DGlass[,-10],   method="lda")
lda.step <- lda(F ~ Mg + Al, data=DGlass[,-10]) 
partimat(F ~ Mg + Al, data=DGlass[,-10], main='',
              method="lda") 
Out_CTab(lda.step, DGlass$F) 

library (caret)
ldaProfile <- rfe(DGlass[,1:9], DGlass$F1,
sizes = 2:9,
rfeControl = rfeControl(functions = ldaFuncs, 
             method = "repeatedcv",repeats = 6))

# Модель на основе всех 9 предикторов
lda.full.pro <- train(DGlass[, 1:9], DGlass$F1, 
        data= DGlass,method="lda",
  trControl = trainControl(method="repeatedcv",repeats=5, 
  classProbs = TRUE), metric = "Accuracy")
# Модель на основе 2 предикторов  stepclass
lda.step.pro <- train(F1 ~ Mg + Al, data= DGlass,method="lda",
  trControl = trainControl(method="repeatedcv",repeats=5, 
  classProbs = TRUE), metric = "Accuracy")
# Модель на основе 3 предикторов rfe
lda.rfe.pro <- train(F1 ~ Al + K + Fe, 
       data= DGlass,method="lda",
    trControl = trainControl(method="repeatedcv",repeats=5, 
    classProbs = TRUE), metric = "Accuracy")
plot(varImp(lda.full.pro))

#-----------------------------------------------------------------------
#  6.2. Метод опорных векторов
#-----------------------------------------------------------------------

DGlass <- read.table(file="Glass.txt", sep=",",
                header=TRUE,row.names=1)
DGlass$F <- as.factor(ifelse(DGlass$Class==2,2,1))
svm.all <-  svm(formula = F ~ ., data=DGlass[,-10], 
       cross=10, kernel = "linear")  
table(Факт=DGlass$F, Прогноз=predict(svm.all))
Acc = mean(predict(svm.all) == DGlass$F )
paste("Точность=", round(100*Acc, 2), "%", sep="")

#   Функция вычисления ошибки кросс-проверки  
library(e1071)    
CVsvm <- function(x, y) {
   n <- nrow(x) ; Err_S  <- 0 
   for(i in 1:n) {
      svm.temp <- svm(x=x[-i,], y = y[-i],  kernel = "linear")
      if (predict(svm.temp, newdata=x[i,]) != y[i]) 
             Err_S <- Err_S +1 }  
   Err_S/n }
Acc <- 1-CVsvm(DGlass[,1:9],DGlass$F)
paste("Точность=", round(100*Acc, 2), "%", sep="")

library ( klaR)
stepclass(F ~ ., data=DGlass[,-10], method="svmlight", 
          pathsvm = "D:/R/SVMlight")

source("SVM-RFE.r")
featureRankedList = svmrfeFeatureRanking(DGlass[,1:9],
             DGlass$F)
ErrSvm <- sapply(1:9, function(nf)  {
     svmModel = svm(DGlass[,featureRankedList[1:nf]],
     DGlass$F, 	kernel="linear") 
     mean( predict(svmModel) != DGlass$F )  }  )
data.frame(Index=featureRankedList, NameFeat =
       names(DGlass[,featureRankedList]), ErrSvm = ErrSvm)

DGlass.sel <- DGlass[, c(featureRankedList[1:7], 11)]
(svm.sel <- svm(formula = F ~ ., data=DGlass.sel, cross=10,
       kernel = "linear", prob=TRUE))
table(Факт=DGlass.sel$F, Прогноз=predict(svm.sel))
Acc = mean(predict(svm.all) == DGlass$F )
paste("Точность=", round(100*Acc, 2), "%", sep="")

library (caret)
DGlass$F1 <- as.factor(ifelse(DGlass$Class==2,"No","Flash"))
svmProfile <- rfe(DGlass[,1:9], DGlass$F1, sizes = 2:9,
    rfeControl = rfeControl(functions = caretFuncs,
    method = "cv"), method = "svmLinear")

ctrl <- trainControl(method="repeatedcv",repeats=5,	
     summaryFunction=twoClassSummary, classProbs=TRUE)
print(" Модель на основе всех 9 предикторов ")
train(DGlass[,1:9], DGlass$F1, method = "svmLinear",
              metric="ROC",trControl=ctrl)	
print(" Модель на основе 7 предикторов ")
train(DGlass[,c(4,6,9,8,3,5,7)], DGlass$F1,
   method = "svmLinear", metric="ROC", trControl=ctrl)
print(" Модель на основе 5 предикторов ")
train(DGlass[,c(4,6,9,8,5)], DGlass$F1, method = "svmLinear",
              metric="ROC",trControl=ctrl)

#-----------------------------------------------------------------------
#  6.3. Классификаторы с использованием нелинейных разделяющих поверхностей
#-----------------------------------------------------------------------

DGlass <- read.table(file="Glass.txt", sep=",",
                header=TRUE,row.names=1)
DGlass$F <- as.factor(ifelse(DGlass$Class==2,2,1))
library(e1071)
tune.svm(F ~ ., data=DGlass[,-10], gamma = 2^(-1:1), 
                      cost = 2^(2:4))

svm.rbf <- svm(formula = F ~ ., data=DGlass[,-10], 
   		kernel = "radial",cost = 4, gamma = 0.5)
table(Факт=DGlass$F, Прогноз=predict(svm.rbf)) 
Acc <- mean(predict(svm.rbf) == DGlass$F )
paste("Точность=", round(100*Acc, 2), "%", sep="") 

svm2.rbf <- svm(formula = F ~ Al+Mg, data=DGlass[,-10], 
   kernel = "radial",cost = 4, gamma = 2.5)
Acc <- mean(predict(svm2.rbf) == DGlass$F )
paste("Точность=", round(100*Acc, 2), "%", sep="")
plot(svm2.rbf,DGlass[,c(3,4,11)],
svSymbol = 16, dataSymbol = 2,color.palette = terrain.colors)
legend("bottomleft",c("Опорные кл.1", "Опорные кл.2", 
      "Данные кл.1", "Данные кл.2"),col=c(1,2,1,2),
       pch=c(16,16,2,2))

library (caret)
ctrl <- trainControl(method="repeatedcv",repeats=5,	
   summaryFunction=twoClassSummary, classProbs=TRUE)
grid <- expand.grid(sigma = (3:6)/10, C = 2:5)
svmRad.tune <- train(DGlass[,1:9], DGlass$F1, 
   method = "svmRadial", metric="ROC",
   tuneGrid = grid, trControl=ctrl)

pred <- predict(svmRad.tune,DGlass[,1:9])
table(Факт=DGlass$F1, Прогноз=pred)
Acc <- mean(pred == DGlass$F1)
paste("Точность=", round(100*Acc, 2), "%", sep="")

library('kernlab') ; library('ggplot2')
set.seed(2335246L) ; data('spirals')
#  Спектральная функция выделяет две различные спирали
sc <- specc(spirals, centers = 2)
s <- data.frame(x=spirals[,1],y=spirals[,2],
class=as.factor(sc))
#  Данные делятся на обучающую и тестовую выборки
s$group <- sample.int(100,size=dim(s)[[1]],replace=T)
sTrain <- subset(s,group>10)
sTest <- subset(s,group<=10)
#  Снова используем Гауссово или радиальное ядро
mSVMG <- ksvm(class~x+y,data=sTrain,kernel='rbfdot')
sTest$predSVMG <- predict(mSVMG,newdata=sTest,type='response')
table(Факт=sTest$class, Прогноз=sTest$predSVMG) 

ggplot() +
   geom_text(data=sTest,aes(x=x,y=y,
      label=predSVMG),size=12) +
   geom_text(data=s,aes(x=x,y=y,
      label=class,color=class),alpha=0.7) +
   coord_fixed() +
   theme_bw() + theme(legend.position='none') 

#-----------------------------------------------------------------------
#  6.4. Деревья классификации, случайный лес и логистическая регрессия
#-----------------------------------------------------------------------

DGlass <- read.table(file="Glass.txt", sep=",",
                header=TRUE,row.names=1)
DGlass$F1 <- as.factor(ifelse(DGlass$Class==2,"No","Flash"))
library (caret)
control<-trainControl(method="repeatedcv",number=10,repeats=3)

# Деревья классификации
set.seed(7)
fit.cart <- train(DGlass[,1:9], DGlass$F1, 
                 method="rpart", trControl=control)
library(rpart.plot)
rpart.plot(fit.cart$finalModel)
pred <- predict( fit.cart,DGlass[,1:9])
table(Факт=DGlass$F1, Прогноз=pred)
Acc <- mean(pred == DGlass$F1)
paste("Точность=", round(100*Acc, 2), "%", sep="")

# Случайный лес (Random Forrest)
set.seed(7)
fit.rf <- train(DGlass[,1:9], DGlass$F1, 
                 method="rf", trControl=control)
fit.rf$finalModel
pred <- predict(fit.rf,DGlass[,1:9])
table(Факт=DGlass$F1, Прогноз=pred)
Acc <- mean(pred == DGlass$F1)
paste("Точность=", round(100*Acc, 2), "%", sep="")

# Логистическая регрессия
DGlass$F1 <- ifelse(DGlass$Class==2,1,0)
logit.all <- glm(F1 ~ .,data=DGlass[,-10], family=binomial)
mp.all <- predict(logit.all, type="response")
Acc <- mean(DGlass$F1 == ifelse(mp.all>0.5,1,0))
paste("Точность=", round(100*Acc, 2), "%", sep="")

library(boot)
cost <- function(r, pi = 0) mean(abs(r-pi) > 0.5)
Acc <- 1 - cv.glm(DGlass[,-10], logit.all, cost)$delta[1]
paste("Точность=", round(100*Acc, 2), "%", sep="")

# Регрессия на информативные переменные
logit.step <- step(logit.all)
summary(logit.step)

# Ошибка на обучающей выборке
mp.step <- predict(logit.step, type="response")
Acc <- mean(DGlass$F1== ifelse(mp.step>0.5,1,0))
paste("Точность=", round(100*Acc, 2), "%", sep="")
# Ошибка при скользящем контроле
Acc <- 1 - cv.glm(DGlass[,-10], logit.step, cost)$delta[1]
paste("Точность=", round(100*Acc, 2), "%", sep="")

#-----------------------------------------------------------------------
#  6.5. Процедуры сравнения эффективности моделей классификации
#-----------------------------------------------------------------------

set.seed(7)
train <- createDataPartition(DGlass$F1, p=0.7)
# определение схемы тестирования
control <- trainControl(method="repeatedcv", 
    number=10, repeats=3, classProbs = T)

#     Выполняем обучение шести моделей
# LDA - линейный дискриминантный анализ
set.seed(7)
fit.lda <- train(DGlass[train,3:4], DGlass$F1[train],
       method="lda", trControl=control)
# SVML - метод опорных векторов с линейным ядром
set.seed(7)
fit.svL <- train(DGlass[train,1:9], DGlass$F1[train],
      method="svmLinear", trControl=control)
# SVMR  - метод опорных векторов с радиальным ядром
set.seed(7)
fit.svR <- train(DGlass[train,1:9], DGlass$F1[train],
    method="svmRadial", trControl=control,
    tuneGrid = expand.grid(sigma = 0.4, C = 2))
# CART - дерево классификации
set.seed(7)
fit.cart <- train(DGlass[train,1:9], DGlass$F1[train],
    method="rpart", trControl=control)
# RF - случайный лес
set.seed(7)
fit.rf <- train(DGlass[train,1:9], DGlass$F1[train],
    method="rf", trControl=control)
# GLM - Логистическая регрессия
set.seed(7)
fit.glm <- train(DGlass[train,-c(1,4,10,11)],DGlass$F1[train],
 method="glm", family=binomial, trControl=control)
caret.models <-  list(LDA=fit.lda, SVML=fit.svL,
    SVMR=fit.svR, CART=fit.cart, RF=fit.rf, GLM=fit.glm)

# ресэмплинг коллекции моделей
results <- resamples(caret.models)
# обобщение различий между моделями
summary(results, metric = "Accuracy")
#  Оценка доверительных интервалов и построение графика
scales <- list(x=list(relation="free"),
                        y=list(relation="free"))
dotplot(results, scales=scales) 
diffs <- diff(results)
# summarize p-values for pair-wise comparisons
summary(diffs) 

#  прогноз на тестовой выборкей
pred.fm <- data.frame(
    LDA=predict(fit.lda,DGlass[-train,3:4]),
    SVML=predict(fit.svL, DGlass[-train,1:9]),
    SVMR=predict(fit.svR, DGlass[-train,1:9]),
    CART=predict(fit.cart, DGlass[-train,1:9]),
    RF=predict(fit.rf, DGlass[-train,1:9]),
    GLM=predict(fit.glm, DGlass[-train,-c(1,4,10,11)])
)
CombyPred <- apply(pred.fm, 1, function (voice) {
     voice2 <- c(voice, voice[5]) # У RF двойной голос
     ifelse(sum(voice2=="Flash") > 3,"Flash", "No") }
         )
pred.fm <- cbind(pred.fm, COMB=CombyPred)
head(pred.fm)

# Функция формирования строки критериев
ModCrit <- function (fact, pred) {
cM <- table(fact, pred)
c( Accur <- (cM[1,1]+cM[2,2])/sum(cM),
   Sens <- cM[1,1]/(cM[1,1]+cM[2,1]),
   Spec <- cM[2,2]/(cM[2,2]+cM[1,2]),
   F.score <-  2 * Accur * Sens / (Accur + Sens),
   MCC <- ((cM[1,1]*cM[2,2])-(cM[1,2]*cM[2,1])) / 
        sqrt((cM[1,1]+cM[1,2])*(cM[1,1]+cM[2,1])*
        (cM[2,2]+cM[1,2])*(cM[2,2]+cM[2,1])) )
}
Result <- t(apply(pred.fm,2, function (x) 
                   ModCrit(DGlass$F1[-train], x)))
colnames(Result) <- c("Точность","Чувствит.","Специфичн.",
            "F-мера","КК Мэтьюса")
round(Result,3)

# Извлекаем вероятности классов по всей выборке
pred.lda.roc <- predict(fit.lda,DGlass[,3:4],type = "prob")
pred.svmR.roc <- predict(fit.svR,DGlass[,1:9],type = "prob")
pred.rf.roc <- predict(fit.rf,DGlass[,1:9],type = "prob")
library(pROC)
# Строим три ROC-модели
m1.roc <- roc(DGlass$F1, pred.lda.roc[,1])
m2.roc <- roc(DGlass$F1, pred.svmR.roc[,1])
m3.roc <- roc(DGlass$F1, pred.rf.roc[,1])
plot(m1.roc, grid.col=c("green", "red"), grid=c(0.1, 0.2),
      print.auc=TRUE,print.thres=TRUE)
plot(m2.roc , add = T, col="green", print.auc=T,
      print.auc.y=0.45,print.thres=TRUE)
plot(m3.roc , add = T, col="blue", print.auc=T,
      print.auc.y=0.40,print.thres=TRUE)
legend("bottomright", c("LDA","SVM","RF"),lwd=2,
      col=c("black","green","blue"))
# Доверительные интервалы для параметров ROC-анализа:
ci.auc(m2.roc); ci.auc(m3.roc)
roc.test(m2.roc,m3.roc)

weight <- c(as.numeric(m1.roc$auc)^2,as.numeric(m2.roc$auc)^2,
         as.numeric(m3.roc$auc)^2)
weight <- weight/sum(weight)
Pred=data.frame(lda=pred.lda.roc[,1],svm=pred.svmR.roc[,1],
            rf=pred.rf.roc[,1])
CombyPred <- apply(Pred,1,function (x) sum(x*weight))
m.Comby <-roc(DGlass$F1,  CombyPred)
plot(m.Comby, add = T, col="red", print.auc=T,
      print.auc.y=0.35)
coords(m.Comby, "best",
 ret=c("threshold", "specificity", "sensitivity", "accuracy"))















#-----------------------------------------------------------------------
#  6.4. Деревья классификации, случайный лес и логистическая регрессия
#-----------------------------------------------------------------------


library (caret)
DGlass <- read.table(file="Glass.txt", sep=",",
                header=TRUE,row.names=1)
DGlass$F1 <- as.factor(ifelse(DGlass$Class==2,"No","Flash"))

library(MASS)
lda.step.pro <- train(DGlass[,3:4], DGlass$F1,metod="lda",
  trControl = trainControl(method = "cv", classProbs = TRUE), metric = "Accuracy")
lda.step.pro$finalModel

# Result
Random Forest 
163 samples
  2 predictor
  2 classes: 'Flash', 'No' 
No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 147, 147, 147, 146, 146, 147, ... 
Resampling results:
  Accuracy   Kappa    
  0.8047794  0.6122885
Tuning parameter 'mtry' was held constant at a value of 2
# Result
Call:
 randomForest(x = x, y = y, mtry = param$mtry, metod = "lda") 
               Type of random forest: classification
                     Number of trees: 500
No. of variables tried at each split: 2

        OOB estimate of  error rate: 15.95%
Confusion matrix:
      Flash No class.error
Flash    71 16   0.1839080
No       10 66   0.1315789
# Result

lm.step.pro <- train(F1 ~., data= DGlass[,-c(1,4,10,12)], metod="glm",
 family=binomial, trControl = trainControl(method = "cv", classProbs = T,
 summaryFunction = twoClassSummary ))

# Result
Random Forest 

163 samples
  7 predictor
  2 classes: 'Flash', 'No' 

No pre-processing
Resampling: Cross-Validated (10 fold) 
Summary of sample sizes: 147, 146, 147, 147, 146, 147, ... 
Resampling results across tuning parameters:

  mtry  ROC        Sens       Spec     
  2     0.9176339  0.8402778  0.8035714
  4     0.9178323  0.8527778  0.8142857
  7     0.9221726  0.8527778  0.8535714

ROC was used to select the optimal model using  the largest value.
The final value used for the model was mtry = 7. 

# Result
Call:
 randomForest(x = x, y = y, mtry = param$mtry, metod = "glm",      family = ..2) 
               Type of random forest: classification
                     Number of trees: 500
No. of variables tried at each split: 7

        OOB estimate of  error rate: 11.66%
Confusion matrix:
      Flash No class.error
Flash    77 10   0.1149425
No        9 67   0.1184211
# Result