#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 4. ПОСТРОЕНИЕ РЕГРЕССИОННЫХ МОДЕЛЕЙ РАЗЛИЧНОГО ТИПА 
#########################################################################

#  4.1. Селекция оптимального набора предикторов линейной модели
#-----------------------------------------------------------------------

library(DMwR)
data(algae)
library(caret)
# Заполним пропуски в данных на основе алгоритма бэггинга
pPbI <- preProcess(algae[,4:11], method='bagImpute')
algae[,4:11] <-  predict(pPbI, algae[,4:11])
# Сохраним таблицу для использования в дальнейшем
save(algae, file="algae.RData") 

Mcor <-cor(algae[,4:11])
library(corrplot)
corrplot(Mcor, method="color", addCoef.col="green", 
        addgrid.col = "gray33",tl.col = "black")

lm.a1 <- lm(a1 ~ .,data=algae[,1:12])
summary(lm.a1)

lm_step.a1 <- step(lm.a1, trace=0)
summary(lm_step.a1)

anova(lm.a1, lm_step.a1)

lm.a1.cv <- train(a1 ~ .,data=algae[,1:12], method='lm',
       trControl = trainControl(method = "cv"))
lm_step.a1.cv <- train( a1 ~ size + mxPH + mnO2 + NO3 + NH4 +
            PO4,data=algae[,1:12], method='lm',
            trControl = trainControl(method = "cv"))

x <- model.matrix(a1 ~ .,data=algae[,1:12])[,-1]
set.seed(10)
#--Рекурсивное исключение переменных
ctrl <- rfeControl(functions = lmFuncs, method = "cv", 
verbose = FALSE,    returnResamp = "final")
lmProfileF <- rfe(as.data.frame(x), algae$a1, 
sizes = 1:10, rfeControl = ctrl)

ggplot(lmProfileF, metric = "Rsquared")
predictors(lmProfileF)
summary(lmProfileF$fit) 

#--Генетический алгоритм
set.seed(10)
ctrl <- gafsControl(functions = rfGA, method = "cv", 
verbose = FALSE,     returnResamp = "final")
lmProfGafs <- gafs(algae[,1:11], algae$a1, 
    iters = 10, # 10 generations of algorithm
    gafsControl = ctrl)
lmProfGafs
lm_gafs.a1 <- lm(a1 ~ size + mxPH +Cl + NO3 + PO4,
 data=algae[,1:12])
train(a1 ~ size + mxPH +Cl + NO3 + PO4,
data=algae[,1:12], method='lm',
trControl = trainControl(method = "cv"))

#--Тестирование моделей 
Eval <- read.table('Eval.txt',header=F,dec='.',
   col.names=c('season','size','speed','mxPH','mnO2','Cl',
  'NO3','NH4','oPO4','PO4','Chla'),na.strings=c('XXXXXXX'))
Sols <- read.table('Sols.txt',header=F,dec='.',          
          col.names=c('a1','a2','a3','a4','a5','a6','a7'),
          na.strings=c('XXXXXXX'))
ImpEval <- preProcess(Eval[,4:11], method='bagImpute')
Eval[,4:11] <-  predict(ImpEval, Eval[,4:11])
# Cохраним данные для дальнейшего использования
save(Eval,Sols, file="algae_test.RData")
y <- Sols$a1
EvalF <- as.data.frame(model.matrix(y ~ .,Eval)[,-1])
# Функция, выводящая вектор критериев
ModCrit <- function (pred, fact) {
  mae <- mean(abs(pred-fact))
  rmse <- sqrt(mean((pred-fact)^2))
  Rsq <- 1-sum((fact-pred)^2)/sum((mean(fact)-fact)^2)
  c(MAE=mae, RSME=rmse,  Rsq = Rsq )
} 
y <- Sols$a1
EvalF <- as.data.frame(model.matrix(y ~ .,Eval)[,-1])
Result <- rbind(
   lm_full = ModCrit(predict(lm.a1,Eval),Sols[,1]),
   lm_step = ModCrit(predict(lm_step.a1,Eval),Sols[,1]),
   lm_rfe = ModCrit(predict(lmProfileF$fit,EvalF),Sols[,1]),
   lm_gafs = ModCrit(predict(lm_gafs.a1,Eval),Sols[,1]))
Result

#-----------------------------------------------------------------------
#  4.2. Регуляризация, частные наименьшие квадраты и kNN-регрессия 
#-----------------------------------------------------------------------
#  Регрессия Лассо
load(file="algae.RData") # Загрузка таблицы algae - раздел 4.1
grid=10^seq(10,-2,length=100)
library( glmnet)
lasso.a1 <- glmnet(x,algae$a1, alpha=1, lambda=grid)  
plot(lasso.a1, xvar =  "lambda", label = TRUE, lwd=2)

x <- model.matrix(a1 ~ .,data=algae[,1:12])[,-1]
grid.train = seq(0.5,4.5,length=15)
lasso.a1.train <- train(as.data.frame(x), algae$a1,
    method='glmnet',
    tuneGrid = expand.grid(.lambda = grid.train, .alpha = 1),
    trControl = trainControl(method = "cv"))

coef(lasso.a1.train$finalModel,
 lasso.a1.train$bestTune$lambda)

#  Метод частных наименьших квадратов (PLS)
library(pls)
M.pls <- plsr(algae$a1~x, scale=TRUE, 
                validation="CV", method="oscorespls")
 summary(M.pls)

set.seed(100)
ctrl <- trainControl(method = "cv", number = 10)
plsTune.a1 <- train(x, algae$a1, method = "pls",
              tuneLength = 14, trControl = ctrl,
              preProc = c("center", "scale"))

set.seed(100)
ctrl <- trainControl(method = "cv", number = 10)
pcrTune.a1 <- train(x, algae$a1, method = "pcr",
         tuneLength = 14, trControl = ctrl,
         preProc = c("center", "scale"))

#  Регрессия по методу k ближайших соседей
knnTune.a1 <- train(x, algae$a1, method = "knn",
         preProc = c("center", "scale"),
         trControl = ctrl, tuneGrid = data.frame(.k = 4:30))
plot( knnTune.a1)

#  Тестирование моделей 
load (file="algae_test.RData") # Загрузка таблиц Eval,Sols
y <- Sols$a1
EvalF <- as.data.frame(model.matrix(y ~ .,Eval)[,-1])
# Функция, выводящая вектор критериев
ModCrit <- function (pred, fact) {
  mae <- mean(abs(pred-fact))
  rmse <- sqrt(mean((pred-fact)^2))
  Rsq <- 1-sum((fact-pred)^2)/sum((mean(fact)-fact)^2)
  c(MAE=mae, RSME=rmse,  Rsq = Rsq )
} 
Result <- rbind(
lasso = ModCrit(predict(lasso.a1.train,EvalF),Sols[,1]),
 pls_1c = ModCrit(predict(plsTune.a1, EvalF),Sols[,1]),
 pcr_1c = ModCrit(predict(pcrTune.a1, EvalF),Sols[,1]),
 kNN_21 = ModCrit(predict(knnTune.a1, EvalF),Sols[,1]))
Result

#-----------------------------------------------------------------------
#  4.3. Построение деревьев регрессии
#-----------------------------------------------------------------------

#  Построение деревьев рекурсивного разделения
load(file="algae.RData") # Загрузка таблицы algae - раздел 4.1
(rt.a1 <- rpart(a1 ~ .,data=algae[,1:12]))
prettyTree(rt.a1) 
printcp(rt.a1) 

#   Снижаем порог стоимости сложности с шагом .005
rtp.a1 <- rpart(a1 ~ .,data=algae[,1:12], 
   control=rpart.control(cp=.005)) 
#  График изменения относительных ошибок от числа узлов дерева
plotcp(rtp.a1) 
with(rtp.a1, {lines(cptable[,2]+1,cptable[,3],
   type="b",col="red")
   legend(locator(1),c("Ошибка обучения",
   "Ошибка крос-проверки (CV)","min(CV ошибка)+SE"),
   lty=c(1,1,2),col=c("red","black","black"),bty="n") })
rtp.a1 <- prune(rtp.a1, cp=0.029)
prettyTree(rtp.a1) 

library(caret)
cvCtrl <- trainControl(method = "repeatedcv", repeats = 3)
rt.a1.train <- train(a1 ~ .,data=algae[,1:12], 
   method = "rpart", tuneLength = 30, trControl = cvCtrl)
plot(rt.a1.train)
rtt.a1 <- rt.a1.train$finalModel
prettyTree(rtt.a1) 

#   Построение деревьев с использованием алгортма условного вывода
library(party)  # Построение дерева методом "условного вывода"
(ctree.a1 <- ctree(a1 ~ .,data=algae[,1:12]))
plot(ctree.a1)

ctree.a1.train <- train(a1 ~ .,data=algae[,1:12], 
    method = "ctree", tuneLength = 10, trControl = cvCtrl)
ctreet.a1 <- ctree.a1.train$finalModel
plot(ctreet.a1)

#  Тестирование моделей 
load (file="algae_test.RData") # Загрузка таблиц Eval,Sols
# Функция, выводящая вектор критериев
ModCrit <- function (pred, fact) {
  mae <- mean(abs(pred-fact))
  rmse <- sqrt(mean((pred-fact)^2))
  Rsq <- 1-sum((fact-pred)^2)/sum((mean(fact)-fact)^2)
  c(MAE=mae, RSME=rmse, Rsq=Rsq) 
} 
Result <- rbind(
 rpart_prune = ModCrit(predict(rtp.a1,Eval),Sols[,1]),
 rpart_train = ModCrit(predict(rt.a1.train,Eval),Sols[,1]),
 ctree_party = ModCrit(predict(ctree.a1,Eval),Sols[,1]),
 ctree_train = ModCrit(predict(ctree.a1.train,Eval),Sols[,1])
)
Result

#-----------------------------------------------------------------------
#  4.4. Ансамбли моделей: бэггинг, случайные леса, бустинг
#-----------------------------------------------------------------------
load(file="algae.RData") # Загрузка таблицы algae - раздел 4.1
x <- as.data.frame(model.matrix(a1~.,data=algae[,1:12])[,-1])
library(randomForest)
randomForest(x, algae$a1, mtry= ncol(x))

bag.a1 <- train(x, algae$a1,
     preProc=c('center', 'scale'),
     method='rf',trControl = trainControl(method = "cv"), 
     tuneGrid =expand.grid(.mtry=ncol(x)))

ranfor.a1 <- train(x, algae$a1,
     preProc=c('center', 'scale'),
     method='rf',trControl = trainControl(method = "cv"), 
     tuneGrid =expand.grid(.mtry=2:10),
            importance=TRUE)

varImpPlot(ranfor.a1$finalModel)

plot(ranfor.a1$finalModel, col="blue", lwd=2)
plot(bag.a1$finalModel, col="green", lwd=2, add=TRUE)
legend("topright",c("Bagging", "RandomForrest"),
          col=c("green","blue"), lwd=2)

#  Бустинг
library(gbm)
set.seed(1)
xd <- cbind(a1 = algae$a1, x)
boost.a1=gbm(a1 ~ ., data= xd, distribution="gaussian",
n.trees=1000,interaction.depth=3)
summary(boost.a1)

pred=predict(boost.a1,x,n.trees=1000)
mean((pred-algae$a1)^2)

modelLookup("gbm")
gbmFit.a1 <- train(a1 ~ ., data= xd, 
    method = "gbm", trControl = trainControl(method = "cv"), 
    tuneGrid = expand.grid(.shrinkage = c(0.1,0.05,0.02),
            .interaction.depth=2:5, .n.trees = 50),
    verbose = FALSE)

modelLookup("bstTree")
library(bst)                           
boostFit.a1 <- train(a1 ~ ., data= xd, 
      method='bstTree', trControl=trainControl(method = "cv"), 
      preProc=c('center','scale'))
plot(boostFit)

#  Тестирование моделей 
load (file="algae_test.RData") # Загрузка таблиц Eval,Sols
y <- Sols$a1
EvalF <- as.data.frame(model.matrix(y ~ .,Eval)[,-1])
# Функция, выводящая вектор критериев
ModCrit <- function (pred, fact) {
  mae <- mean(abs(pred-fact))
  rmse <- sqrt(mean((pred-fact)^2))
  Rsq <- 1-sum((fact-pred)^2)/sum((mean(fact)-fact)^2)
  c(MAE=mae, RSME=rmse, Rsq=Rsq) 
} 
Result <- rbind(
   bagging = ModCrit(predict(bag.a1,EvalF),Sols[,1]),
   ranfor = ModCrit(predict(ranfor.a1,EvalF),Sols[,1])
   bst.gbm = ModCrit(predict(gbmFit.a1,EvalF),Sols[,1]),
   bst.bst = ModCrit(predict(boostFit.a1,EvalF),Sols[,1]))
Result

#-----------------------------------------------------------------------
#  4.5. Сравнение построенных моделей и оценка информативности предикторов
#-----------------------------------------------------------------------
Models <- read.delim('Models.txt',header=T)
plot(Models$Rsq,Models$Rsquared,pch=CIRCLE<-16, 
     col=8-Models$col, cex=2.5,
     xlab="Rsquared на дополнительной выборке", 
     ylab="Rsquared при кросс-проверке")
text(Models$Rsq,Models$Rsquared,rownames(Models), 
     pos=4, font=4,cex=0.8)
legend('bottomright',c('Бэггинг/бустинг','Деревья',
     'Регрессия kNN','PLS/PCR','Лассо','Линейные модели'),
     col=2:7, pch=CIRCLE<-16, cex=1)

load(file="algae.RData") # Загрузка таблицы algae - раздел 4.1
library(Boruta)
algae.mod <- as.data.frame(model.matrix(a1 ~ .,
                    data=algae[,1:12])[,-1])
algae.mod <- cbind(algae.mod, a1=algae$a1)
algae.Boruta <- Boruta(a1 ~ ., data = algae.mod, 
            doTrace = 2, ntree = 500
getConfirmedFormula(algae.Boruta) 
attStats(algae.Boruta)

plot(algae.Boruta, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(algae.Boruta$ImpHistory),function(i)
   algae.Boruta$ImpHistory[is.finite(algae.Boruta$ImpHistory[,i])
                  ,i])
names(lz) <- colnames(algae.Boruta$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(algae.Boruta$ImpHistory), cex.axis = 0.7

#-----------------------------------------------------------------------
#  4.6. Деревья регрессии с многомерным откликомв
#-----------------------------------------------------------------------
load(file="algae.RData") # Загрузка таблицы algae - раздел 4.1
library(caret)
Transal = preProcess(algae[,12:18], 
method = c("BoxCox", "scale"))
Species = predict(Transal, algae[,12:18])

library(mvpart) 
spe.mvpart <- mvpart(data.matrix(Species) ~ ., algae[,1:11], 
xv="pick", xval=nrow(Species), xvmult = 1,
         margin=0.08, which=4, bars = TRUE)
plot(spe.mvpart) ; text(spe.mvpart)  
summary(spe.mvpart)

# Относительная доля групп водорослей в кластерах
groups.mrt <- levels(as.factor(spe.mvpart$where))
leaf.sum <- matrix(0, length(groups.mrt), ncol(Species))  
 colnames(leaf.sum) <- colnames(Species)
 rownames(leaf.sum) <- groups.mrt
for(i in 1:length(groups.mrt)){
     leaf.sum[i,] <- apply(Species[which
    (spe.mvpart$where==groups.mrt[i]),], 2, sum)
   }  
leaf.sum
opar <- par() 
#  Вывод диаграммы типа "разрезанный пирог"
par(mfrow=c(2,2)) ; for(i in 1:length(groups.mrt)){
pie(leaf.sum[i,which(leaf.sum[i,]>0)], radius=1, 
main = paste("Кл. №", groups.mrt[i])) }
par(opar)
#  Вывод диаграммы РСА
rpart.pca(spe.mvpart)















#--







