#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 7. МОДЕЛИ КЛАССИФИКАЦИИ ПРИ КАТЕГОРИАЛЬНОМ ОТКЛИКЕ
#########################################################################

#  7.1. Ирисы Фишера и метод k-ближайших соседей
#-----------------------------------------------------------------------
data(iris) ; library(ggplot2)
qplot(Sepal.Length, Sepal.Width, data = iris) +
      facet_grid(facets = ~ Species)+
      geom_smooth(color="red", se = FALSE) 
qplot(Petal.Length, Petal.Width,data = iris) +
      facet_grid(facets = ~ Species)+
      geom_smooth(color="red", se = FALSE) 

library(vegan) 
mod.pca <- rda(iris[,-5], scale = TRUE)
scores <- as.data.frame(scores(mod.pca, display ="sites",
 scaling = 3))
scores$Species <- iris$Species
# Включаем в названия осей доли объясненных дисперсий
axX <- paste("PC1 (",
   as.integer(100*mod.pca$CA$eig[1]/sum(mod.pca$CA$eig)),"%)")
axY <- paste("PC2 (",
   as.integer(100*mod.pca$CA$eig[2]/sum(mod.pca$CA$eig)),"%)")
# Составляем таблицу для hull "каркаса точек на графике"
l <- lapply(unique(scores$Species), function(c) 
         { f <- subset(scores,Species==c); f[chull(f),]})
hull <- do.call(rbind, l)
# Выводим ординационную диаграмму
ggplot() + 
  geom_polygon(data=hull,aes(x=PC1,y=PC2, fill=Species),
alpha=0.4, linetype=0) +  
  geom_point(data=scores,aes(x=PC1,y=PC2,shape=Species,
colour=Species),size=3) + 
  scale_colour_manual( values = c('purple', 'green', 'blue'))+
  xlab(axX) + ylab(axY) + coord_equal() + theme_bw() 

library(kknn)   #   -------- Метод к-ближайших соседей  kNN
train.kknn(Species ~ ., iris, kmax = 50, kernel="rectangular") 
max_K=20 ; gen.err.kknn <- numeric(max_K)
mycv.err.kknn <- numeric(max_K) ; n <- nrow(iris)
#  Рассматриваем число возможных соседей от 1 до 20
for (k.val in 1:max_K)  {  
  pred.train <- kknn(Species ~ .,
     iris,train=iris,test=iris,k=k.val,kernel="rectangular")
  gen.err.kknn[k.val] <- mean(pred.train$fit != iris$Species)
  for (i in 1:n)   {
   pred.mycv <- kknn(Species~., train=iris[-i,],test=iris[i,],
         k=k.val,kernel="rectangular")
   mycv.err.kknn[k.val] <- mycv.err.kknn[k.val] +
        (pred.mycv$fit != iris$Species[i])   }
}  ;  mycv.err.kknn <- mycv.err.kknn/n
plot(1:20,gen.err.kknn,type="l",xlab='k', ylim=c(0,0.07),
    ylab='Ошибка классификации', col="limegreen", lwd=2) 
points(1:max_K,mycv.err.kknn,type="l",col="red", lwd=2)
legend("bottomright",c("При обучении",
   "Скользящий контроль"),lwd=2,col=c("limegreen", "red")) 

library(caret) ; set.seed(123)
contrl <- trainControl(method="repeatedcv",repeats = 3)
train(Species ~ ., data = iris, method = "knn", 
       trControl = contrl, preProcess = c("center","scale"),
       tuneLength = 20)

set.seed(123) ;  (samp.size <- floor(nrow(iris) * .75))
train.ind <- sample(seq_len(nrow(iris)), size = samp.size)
train <- iris[train.ind,]  ;  test <- iris[-train.ind, ] 
knn.iris <- knn(train = train[,-5], test = test[,-5], 
cl = train[,"Species"], k = 13, prob = T)
table(Факт=test$Species,Прогноз=knn.iris)
Acc = mean(knn.iris == test$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")

#-----------------------------------------------------------------------
#  7.2. "Наивный" классификатор Байеса
#-----------------------------------------------------------------------

library(klaR)
naive_iris <- NaiveBayes(iris$Species ~ ., data = iris)
naive_iris$tables$Petal.Width 
plot(naive_iris,lwd=2)

pred <- predict(naive_iris,iris[,-5])$class
table(Факт=iris$Species, Прогноз=pred)
Acc = mean(pred == iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")

library(caret)
#  Определим условия кросс-проверки в объекте  
train_control <- trainControl(method='cv',number=10)
Test <- train(Species~., data = iris, trControl=train_control, method="nb")
print(Test)
Acc = mean(predict(Test$finalModel,iris[,-5])$class
                          == iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")

#-----------------------------------------------------------------------
#  7.3. Классификация в линейном дискриминантном пространстве
#-----------------------------------------------------------------------

require(MASS)
LDA.iris <- lda(formula = Species ~ .,data = iris)
LDA.iris$scaling # Коэффициенты линейных лискриминантов
LDA.iris$svd
(prop = LDA.iris$svd^2/sum(LDA.iris$svd^2)) 
prop =percent(prop)
pred <-predict(LDA.iris,newdata = iris)
scores = data.frame(Species = iris$Species,pred$x)
# Выводим ординационную диаграмму
require(ggplot2)
ggplot() +   geom_point(data=scores,
  aes(x=LD1,y=LD2,shape=Species,colour=Species),size=3) + 
  scale_colour_manual( values = c('purple', 'green', 'blue'))+
  labs(x = paste("LD1 (", prop[1], ")", sep=""),
       y = paste("LD2 (", prop[2], ")", sep="")) + theme_bw()

lda(scale(iris[,1:4]), gr = iris$Species)$scaling

train <- sample(1:150, 140)
LDA.iris3 <- lda(Species ~ .,iris, subset = train)
plda = predict(LDA.iris3,newdata = iris[-train, ])
data.frame(Species=iris[-train,5], plda$class, plda$posterior)

Acc = mean(pred$class == iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")

LDA.irisCV <- lda(Species ~ ., data = iris,  CV = TRUE)
table(Факт=iris$Species,Прогноз=LDA.irisCV$class)
Acc = mean(LDA.irisCV$class==iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")

#-----------------------------------------------------------------------
#  7.4. Нелинейные классификаторы в R
#-----------------------------------------------------------------------

# Квадратичный дискриминантный анализ 
require(MASS)
QDA.iris = qda(Species~ Petal.Length+Petal.Width, data = iris)
pred = predict(QDA.iris)$class
table(Факт=iris$Species,Прогноз=pred)
Acc = mean(pred ==iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")
library(klaR)
partimat(Species ~ Petal.Length + Petal.Width,
         data = iris, method="qda")

library(caret)
set.seed(123)
# Используем только размеры лепестка
train(Species ~ Petal.Length + Petal.Width, data = iris, method = "qda", trControl = trainControl(method = "cv"))

# Используем весь набор признаков
train(Species ~ ., data = iris, method = "qda", 
trControl = trainControl(method = "cv"))

# Регуляризованный дискриминантный анализ 
library(rda)
set.seed(123)
#  1- этап грубой оптимизации
train(Species ~ ., data = iris, method = "rda", 
trControl = trainControl(method = "cv"))

#  2- этап с диапазоном lambda = 0.1:0.5 и gamma = 0.02:0.1
RDAGrid <- expand.grid(.lambda = (1:5)/10, .gamma = (1:5)/50)
train(Species ~ ., data = iris, method = "rda", 
tuneGrid =RDAGrid, trControl = trainControl(method = "cv"))

#  Построение классификатора и оценка его точности
RDA.iris <- rda(Species~., data=iris, gamma=0.02, lambda=0.5)
pred = predict(RDA.iris)$class
Acc = mean(pred ==iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")

# Машина опорных векторов 
library(e1071)
SVM.iris <- svm(Species~., data=iris)
pred <- predict(SVM.iris, iris[,1:4], type="response")
table(Факт=iris$Species,Прогноз=pred)
Acc = mean(pred ==iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")
pred.DV <- predict(SVM.iris, iris[,1:4], 
decision.values = TRUE)
ind <- sort(sample(150,9))
data.frame(attr(pred.DV, "decision.values")[ind,],
                iris$Species[ind])

plot(cmdscale(dist(iris[,-5])), col = as.integer(iris[,5])+1,
     pch = c("o","+")[1:150 %in% SVM.iris$index + 1], font=2,
xlab="Шкала 1",ylab="Шкала 1" ) 
legend (0,1.2, c("setosa","versicolor","virginica"),pch = "o", 
        col =2:4)

#-----------------------------------------------------------------------
# 7.5. Мультиномиальная логистическая регрессия
#-----------------------------------------------------------------------

library(nnet)
MN.iris <- multinom(Species~., data=iris)
summary(MN.iris)
Probs <- fitted(MN.iris) 
pred=apply(Probs,1,function(x)
        colnames(Probs)[which(x==max(x))])
table(Факт=iris$Species,Прогноз=pred)
Acc = mean(pred ==iris$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")
z <- summary(MN.iris)$coefficients/
       summary(MN.iris)$standard.errors
# p-значения на основе теста Вальда 
(1 - pnorm(abs(z), 0, 1))*2 
# p-значения на основе  t-статистики
pt(z, df = nrow(iris) - 5, lower=FALSE) 

#-----------------------------------------------------------------------
# 7.6. Классификаторы на основе искусственных нейронных сетей
#-----------------------------------------------------------------------

data(iris)
ind = sample(2, nrow(iris), replace = TRUE, prob=c(0.7, 0.3))
trainset = iris[ind == 1,]
testset = iris[ind == 2,]
trainset$setosa = trainset$Species == "setosa"
trainset$virginica = trainset$Species == "virginica"
trainset$versicolor = trainset$Species == "versicolor

library(neuralnet)
net.iris = neuralnet(versicolor + virginica + setosa ~
Sepal.Length + Sepal.Width + Petal.Length + Petal.Width,
trainset, hidden=3)
net.iris$result.matrix
plot(net.iris)

net.prob = compute(net.iris, testset[-5])$net.result
pred = c("versicolor",   "virginica",   "setosa")
 [apply(net.prob,   1,   which.max)]
table(Факт=testset$Species, Прогноз= pred)
Acc = mean(pred == testset$Species)
paste("Точность=", round(100*Acc, 2), "%", sep="")























