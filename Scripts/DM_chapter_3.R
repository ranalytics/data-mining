#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 3. ПАКЕТ СARET  - ИНСТРУМЕНТ ПОСТРОЕНИЯ СТАТИСТИЧЕСКИХ МОДЕЛЕЙ В R 
#########################################################################

#  3.2. Обнаружение и удаление "ненужных" предикторов
#-----------------------------------------------------------------------
library(caret)
data(GermanCredit)
(u <-unique(GermanCredit$ResidenceDuration))
# Доля уникальных значений:
length(u)/nrow(GermanCredit)
(t<- sort(table(GermanCredit$ResidenceDuration), 
                 decreasing = TRUE))
t[1]/t[2]
# Создадим копию данных без столбца с откликом Class:
gcred = GermanCredit[, -10]
# Функция nearZeroVar() возращает вектор номеров переменных,
# обладающих околонулевой дисперсией:
(nz = nearZeroVar(gcred))
print("Имена этих переменных:") ; names(gcred)[nz]
# Удаляем предикторы с околонулевой дисперсией:
gcred.clean = gcred[, -nz]

a <- apply(gcred.clean, 2, max)
ab <- names(a[a==1])
gcred.bin <- gcred.clean[,ab]
head(gcred.bin)
m <-glm(y~., gcred.bin, family=binomial)
summary(m)
library("MASS") 
Fac.lda <- lda(gcred.bin, grouping = y) 


# Наибольшие значения триангулярной матрицы
top.mat <- function(X, level=0.45, N=10, values=TRUE) {
X.nam <- row.names(X)
X.tri <- as.vector(lower.tri(X))
X.rep.g <- rep(X.nam, length(X.nam))[X.tri]
X.rep.e <- rep(X.nam, each=length(X.nam))[X.tri]
X.vec <- as.vector(X)[X.tri]
X.df <- data.frame(Var1=X.rep.g, Var2=X.rep.e, Value=X.vec)
{if (values)
    {X.df <- X.df[abs(X.df$Value) >= level, ]
    X.df <- X.df[order(-abs(X.df$Value)), ]}
else
    {X.df <- X.df[order(-abs(X.df$Value)), ]
    X.df <- X.df[1:N, ]}}
row.names(X.df) <- seq(1, along=X.df$Value)
return(X.df)
}
top.mat(cor(gcred.clean))
# Функция findCorrelation() возращает вектор 
# номеров переменных с высокой корреляцией:
(highCor = findCorrelation(cor(gcred.clean), cutoff = 0.75))
print("Имена этих переменных:")
names( gcred.clean)[highCor]
# Удаляем эти переменные:
gcred.clean =  gcred.clean[, -highCor]
(linCombo <- findLinearCombos(gcred.clean))
# Удаляем эти переменные:
gcred.clean =  gcred.clean[, -linCombo$remove]
dim(gcred.clean)


#-----------------------------------------------------------------------
# 3.3. Препроцессинг: преобразование и групповая трансформация переменных
#-----------------------------------------------------------------------
data(GermanCredit)
TransPred <- c("Duration","Amount","Age")
preVar <- preProcess(GermanCredit[,TransPred])
TransVar = predict(preVar, GermanCredit[,TransPred])
print("До преобразования:")
summary(GermanCredit[,TransPred]) 
print("После преобразования:")
summary(TransVar) 
# при необходимости, обновление переменных
GermanCredit[,TransPred] <- TransVar

preBox <- preProcess(GermanCredit[,TransPred],
             method="BoxCox")
BoxVar <- predict(preBox,GermanCredit[,TransPred]) 
y <- factor(GermanCredit$Class)
trellis.par.set(theme = col.whitebg(), warn = FALSE)
featurePlot(GermanCredit$Amount, y, "density", 
          labels = c("Amount",""))
featurePlot(BoxVar[,2], y, "density", 
          labels = c("Box-Amount",""))

y <- factor(GermanCredit$Class)
gcred = GermanCredit[, -10]
nz = nearZeroVar(gcred)
gcred.clean = gcred[, -nz]
highCor = findCorrelation(cor(gcred.clean), cutoff = 0.75)
gcred.clean =  gcred.clean[, -highCor]
linCombo <- findLinearCombos(gcred.clean)
gcred.clean =  gcred.clean[, -linCombo$remove]
dim(gcred.clean)

library(vegan)
mod.pca <- rda(gcred.clean, scale=TRUE) 
ev <- mod.pca$CA$eig 
# Иллюстрация Критерия Кайзера-Гуттмана
barplot(ev, col="bisque", las=2)
abline(h=mean(ev), col="red") 
legend("topright", "Средние собственные значения", 
             lwd=1, col=2, bty="n")
c(length(ev[ev > mean(ev)]),sum(ev[ev > mean(ev)])/sum(ev)) 

prePCA <- preProcess(gcred.clean, 
          method= c("center", "scale", "pca"), pcaComp=19)
gcred.pca <- predict(prePCA,gcred.clean)

train.index <- createDataPartition(y, p = .8, times = 100)
trControl = trainControl(method="LGOCV", index = train.index)
print('Модель по натуральным предикторам')
modNat <- train(gcred.clean, y, "glm", 
        family=binomial, trControl = trControl)
print('Модель по главным компонентам')
modPCA <- train(gcred.pca, y, "glm", 
        family=binomial, trControl = trControl) 

plot(varImp(modNat, scale = FALSE))

#-----------------------------------------------------------------------
#  3.4. Заполнение пропущенных значений в данных
#-----------------------------------------------------------------------
library(DMwR)
library(ggplot2)
summary(algae) # вывод не приводится
ggplot(data=algae[!is.na(algae$mnO2),], aes(speed , mnO2)) + 
       geom_violin(aes(fill = speed), trim = FALSE, 
       alpha = 0.3) + 
       geom_boxplot(aes(fill = speed), width = 0.2,
       outlier.colour = NA) + 
       theme(legend.position = "NA")

qplot(PO4, a1, data = algae[!is.na(algae$PO4), ]) +
      facet_grid(facets = ~ season)+
      geom_smooth(color="red", se = FALSE) 

# Число строк с пропущенными значениями
nrow(algae[!complete.cases(algae),])
# Их удаление
algae <- na.omit(algae)

data(algae)
manyNAs(algae, 0.2)
algae <- algae[-manyNAs(algae, 0.2), ]

data(algae)
ind<- apply(algae, 1, function(x) sum(is.na(x))) > 0
algae[ind, 4:11]

# Восстановление пропущенных значений медианами
library(caret)
pPmI <- preProcess(algae[, 4:11], method = 'medianImpute')
algae[, 4:11] <- predict(pPmI, algae[, 4:11])
(Imp.Med <- algae[ind, 4:11])

data(algae)
lm(PO4 ~ oPO4, data = algae)
# Функция  вывода значений PO4 в зависимости от оPO4 
fillPO4 <- function(oP) {if (is.na(oP)) return(NA)
      else return(42.897 + 1.293 * oP)
}
# Восстановление пропущенных значений PO4
algae[is.na(algae$PO4), 'PO4'] <- 
   sapply(algae[is.na(algae$PO4), 'oPO4'], fillPO4)
algae[ind, 10]

# Восстановление пропущенных значений бэггингом
data(algae)
pPbI <- preProcess(algae[, 4:11], method = 'bagImpute')
algae[, 4:11] <- predict(pPbI, algae[, 4:11])
Imp.Bag <- algae[ind, 4:11]

# Восстановление пропущенных значений k-ближайшими соседями
data(algae)
pPkI <- preProcess(algae[, 4:11], method = 'knnImpute')
alg.stand <- predict(pPkI, algae[, 4:11])

m <- pPkI$mean
sd <- pPkI$std
algae[, 4:11] <- t(apply(alg.stand, 1, 
                    function (r) m + r * sd))
(Imp.Knn <- algae[ind, 4:11])

ImpVal <- rbind(Imp.Med, Imp.Knn)
ImpVal <- rbind(ImpVal, Imp.Bag)
Imp.Metod <- as.factor(c(rep("Med", 16), rep("Knn", 16), rep("Bag", 16)))

library(vegan)
Imp.M <- rda(ImpVal ~ Imp.Metod, ImpVal)
plot(Imp.M, display = "sites", type = "p")
ordihull(Imp.M, Imp.Metod, draw = "polygon", alpha = 67, 
         lty = 2, col = c(1, 2, 3), label = TRUE)

#-----------------------------------------------------------------------
#  3.5. Функция train() пакета caret
#-----------------------------------------------------------------------
library(caret)
ls(getModelInfo(model = "lm"))
modelLookup("lm")
modelLookup("rpart")

library(DAAG)
data("fruitohms")
set.seed(123) 
max.poly <- 7
degree <- 1:max.poly
RSquared <- rep(0,max.poly)
RMSE <- rep(0,max.poly)
# Выполним 10-кратную кросспроверку с 10 повторностями
fitControl <- trainControl(method = "repeatedcv",
          number = 10,repeats = 10)
# Тестируем модель для различных степеней
for ( d in degree)  {
  LinearRegressor <- train( ohms ~ poly(juice, d),
     data = fruitohms,
     method = "lm",trControl =fitControl)
     RSquared[d] <- LinearRegressor$results$Rsquared
     RMSE[d]<- LinearRegressor$results$RMSE
}
library(ggplot2)
Degree.RegParams = data.frame(degree,RSquared,RMSE)
ggplot(aes(x = degree,y = RSquared),
      data = Degree.RegParams) + geom_line()
ggplot(aes(x = degree,y = RMSE),
      data = Degree.RegParams) + geom_line()
Poly5 <- train(ohms ~ poly(juice,5), data = fruitohms,
      method = "lm")
summary(Poly5$finalModel)
summary(lm(ohms ~ poly(juice,5), data = fruitohms))
Poly5



