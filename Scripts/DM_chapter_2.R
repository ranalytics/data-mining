#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 2. СТАТИСТИЧЕСКИЕ МОДЕЛИ: КРИТЕРИИ И МЕТОДЫ ИХ ОЦЕНИВАНИЯ 
#########################################################################

#  2.1. Основные шаги построения и верификации моделей
#-----------------------------------------------------------------------

# Загрузка и визуализация данных в виде традиционного биплота
library(DAAG) 
data("fruitohms")
library('ggplot2')
ggplot(fruitohms, aes(x=juice, y=ohms)) + geom_point() +
    stat_smooth(method="lm") + 
    xlab("Содержание сока, %") + 
    ylab("Сопротивление, Ом")
    
library(hexbin)
ggplot(fruitohms, aes(x=juice, y=ohms)) +
    geom_hex(binwidth=c(3, 500)) +
    geom_smooth(color="red", se=F) + 
    xlab("Содержание сока, %") + 
    ylab("Сопротивление, Ом")

#  Строим модель lm() с одним предиктором
x <- fruitohms$juice
y <- fruitohms$ohms/1000
n <- dim(as.matrix(x))[1] ; m <- dim(as.matrix(x))[2] 
M_reg <- lm(y~x) ; pred <- predict(M_reg)
#  Показатели модели
RSS <- sum((y-pred)*(y-pred))
RMSE <- sqrt(RSS/n)
RSE <- sqrt(RSS/(n - m -1)) 
c(RSS, RMSE, RSE)
Rsquared <- 1-RSS/sum((mean(y)-y)^2)
F <- (sum((mean(y)-pred)^2)/m)/(RSS/(n - m -1))
p <- pf(q=F, df1=m, df2=(n-m-1), lower.tail=FALSE)
c(Rsquared, F, p)
summary(M_reg)
anova(M_reg)
ncvTest (M_reg)

#  Строим модель glm() с одним предиктором
M_glm <- glm(y~x)
lgLik <- logLik(M_glm)
D.null <- M_glm$null.deviance 
D <- M_glm$deviance
df <- with(M_glm, df.null - df.residual)
p <- pchisq(D.null - D, df, lower.tail = FALSE)
PRsquar = 1 - D/D.null
c(lgLik, D , D.null, p, PRsquar)
summary(M_glm)
with(M_glm, null.deviance - deviance)
anova(M_glm, glm(x ~ 1) , test ="Chisq")
k <- extractAIC(M_glm)[1] ; AIC <- extractAIC(M_glm)[2]
AIC <- AIC(M_glm) ; AICc <- AIC(M_glm) + 2*k*(k+1)/(n-k-1)
BIC <- BIC(M_glm) ;  BIC <- AIC(M_glm, k=log(n))
c(AIC, AICc, BIC)

# Построение моделей со степенью полинома от 1 до 7
max.poly <- 7
# Создание пустой таблицы для хранения значений AIC и BIC 
# рассчитанных для всех моделей и ее заполнение 
AIC.BIC <- data.frame(criterion=c(rep("AIC",max.poly),
rep("BIC",max.poly)), value=numeric(max.poly*2),
degree=rep(1:max.poly, times=2))
for(i in 1:max.poly)  {
     AIC.BIC[i,2] <- AIC(lm(y~poly(x,i)))
     AIC.BIC[i+max.poly,2] <- BIC(lm(y~poly(x,i)))
}
# График AIC и BIC для разных степеней полинома
library('ggplot2')
qplot(degree, value, data=AIC.BIC,
geom="line", linetype=criterion) +
xlab("Степень полинома") + ylab("Значение критерия")

#-----------------------------------------------------------------------
#2.2. Использование алгоритмов ресэмплинга для тестирования и оптимизации параметров моделей
#-----------------------------------------------------------------------
# Функция скользящего контроля для модели y~poly(x, degree) 
crossvalidate <- function(x, y, degree) {
   preds <- numeric(length(x))
   for(i in 1:length(x)) {
      x.in <- x[-i] ; x.out <- x[i]
      y.in <- y[-i] ; y.out <- x[i]
      m <- lm(y.in ~ poly(x.in, degree=degree) )
      new <- data.frame(x.in = seq(-3, 3, by=0.1))
      preds[i]<- predict(m, newdata=data.frame(x.in=x.out))
   }
  # Тестовая статистика - сумма квадратов отклонений:
  return(sum((y-preds)^2))
} 
# Заполнение таблицы результатами кросс-проверки 
# и сохранение квадрата ошибки в таблице "a"
a <- data.frame(cross=numeric(max.poly))
for(i in 1:max.poly)
{
  a[i,1] <- crossvalidate(x, y, degree=i)
}
# График суммы квадратов ошибки при кросспроверке
qplot(1:max.poly,cross, data=a, geom=c("line"))+
xlab("Степень полинома ") + ylab("Квадратичная ошибка")

(M_poly <- glm(y~poly(x, 4)))
anova(M_glm, M_poly , test ="Chisq")

set.seed(123)
Nboot = 1000
BootMat <- sapply(1:max.poly, function(k) 
   replicate(Nboot, {
       ind <- sample(n,replace=T)
       BIC(glm(y[ind]~poly(x[ind],k)))
     }  )
)
apply(BootMat, 2, mean)

library(reshape) # для функции melt()
BootMat.df <- data.frame( melt(BootMat))
ggplot(data = BootMat.df, aes(X2, value)) +
   geom_jitter(position = position_jitter(width = 0.1), 
                alpha = 0.2) + 
# Добавляем средние значения в виде точек красного цвета:
  stat_summary(fun.y = mean, geom = "point", color = "red",
                size = 5)  + 
# Добавляем отрезки, символизирующие 0.25 и 0.75 квантили:
  stat_summary(fun.y = mean,
     fun.ymin = function(x){quantile(x, p = 0.25)},
     fun.ymax = function(x){quantile(x, p = 0.75)},
     geom = "errorbar", color = "magenta", width = 0.5,
             size =1.5) +
     xlab("Степень полинома") + 
     ylab("Информационный критерий Байеса")

Nboot = 1000  ; n=100
mean(replicate(Nboot, length(unique(sample(n,replace=T)))))

#-----------------------------------------------------------------------
# 2.3. Модели предсказания класса объектов 
#-----------------------------------------------------------------------
# spamD <- read.table(
#        'https://raw.github.com/WinVector/zmPDSwR/master/Spambase/spamD.tsv',
#            header=T,sep='\t')
spamD <- read.table('spamD.tsv',header=T,sep='\t')
dim(spamD)
spamTrain <- subset(spamD,spamD$rgroup>=10)
spamTest <- subset(spamD,spamD$rgroup<10)
с(nrow(spamTrain), nrow(spamTest))
# Составляем список переменных и объект типа формула
spamVars <- setdiff(colnames(spamD),list('rgroup','spam'))
spamFormula <- as.formula(paste('spam=="spam"',
   paste(spamVars,collapse=' + '),sep=' ~ '))
spamModel <- glm(spamFormula,family=binomial(link='logit'),
   data=spamTrain)
# Добавляем столбец с прогнозируемыми вероятности спама
spamTrain$pred <- predict(spamModel,newdata=spamTrain,
             type='response')
spamTest$pred  <- predict(spamModel,newdata=spamTest,
             type='response')
#  На обучающей выборке
cM.train <- table(Факт=spamTrain$spam,
   Прогноз=spamTrain$pred>0.5)
#  При экзамене
cM <- table(Факт=spamTest$spam,
   Прогноз=spamTest$pred>0.5)
c( Точность <- (cM[1,1]+cM[2,2])/sum(cM),
   Чувствительность <- cM[1,1]/(cM[1,1]+cM[2,1]),
   Специфичность    <- cM[2,2]/(cM[2,2]+cM[1,2]))
library(caret)
pred <- ifelse(spamTest$pred>0.5,"spam","non-spam") 
confusionMatrix(spamTest$spam, pred) 

#  Информационные показатели
entropy <- function(x) { xpos <- x[x>0]
   scaled <- xpos/sum(xpos) ; sum(-scaled*log(scaled,2))
}
print(entropy(table(spamTest$spam))) 
conditionalEntropy <- function(t) {
(sum(t[,1])*entropy(t[,1]) + sum(t[,2])*entropy(t[,2]))/sum(t)
}
print(conditionalEntropy(cM))

# ROC кривая
library(pROC)
m_ROC.roc <- roc(spamTest$spam,spamTest$pred)
plot(m_ROC.roc, grid.col=c("green", "red"), grid=c(0.1, 0.2),
      print.auc=TRUE, print.thres=TRUE)
plot(smooth(m_ROC.roc), col="blue", add = T, print.auc=F)

# Анализ статзначимости
(LL <- logLik(spamModel))
df <- with(spamModel, df.null - df.residual)
c(D.null <- spamModel$null.deviance,
  D <- spamModel$deviance,
  Rsquared = 1-D/D.null,
  pchisq(D.null - D, df, lower.tail = FALSE))
Null_Model <- glm(spam ~ 1,family=binomial(link='logit'),
   data=spamTrain)
anova(spamModel, Null_Model , test ="Chisq")

ggplot (data=spamTest) +
        geom_density (aes (x=pred, color=spam, linetype=spam))

#-----------------------------------------------------------------------
#  2.4. Проецирование многомерных данных на плоскости
#-----------------------------------------------------------------------
DGlass <- read.table(file="Glass.txt", sep=","
,header=TRUE,row.names=1)
print(t(apply(DGlass[,-10],2,function (x) {
c(Минимум=min(x),Максимум=max(x),
 		Среднее=mean(x), Отклонение=sd(x),
 		Корреляция=cor(x,DClass$Class))  # признаков с Class
})),3)
library(vegan)
Y <- as.data.frame(DGlass[,1:9])
mod.pca <- rda(Y ~ 1)
summary(mod.pca)
F <- as.factor(ifelse(DGlass$Class==2,2,1))
pca.scores <- as.data.frame(summary(mod.pca)$sites [,1:2])
pca.scores <- cbind(pca.scores,F)
# Составляем таблицу для hull "каркаса точек на графике"
l <- lapply(unique(pca.scores$F), function(c) 
         { f <- subset(pca.scores,F==c); f[chull(f),]})
hull <- do.call(rbind, l)
# Включаем в названия осей доли объясненных дисперсий
axX <- paste("PC1 (",
   as.integer(100*mod.pca$CA$eig[1]/sum(mod.pca$CA$eig)),"%)")
axY <- paste("PC2 (",
   as.integer(100*mod.pca$CA$eig[2]/sum(mod.pca$CA$eig)),"%)")
# Выводим ординационную диаграмму
ggplot() + 
  geom_polygon(data=hull,aes(x=PC1,y=PC2, fill=F),
           alpha=0.4, linetype=0) +  
  geom_point(data=pca.scores,aes(x=PC1,y=PC2,shape=F,
           colour=F),size=3) + 
  scale_colour_manual( values = c('purple', 'blue'))+
  xlab(axX) + ylab(axY) + coord_equal() + theme_bw()

#-----------------------------------------------------------------------
#  2.5. Многомерный статистический анализ данных
#-----------------------------------------------------------------------
F <- as.factor(ifelse(DGlass$Class==2,2,1)) 
Y <- as.data.frame(DGlass[,1:9])
mod.rda <- rda(Y ~ F)
summary(mod.rda)
rda.scores <- as.data.frame(scores(mod.rda, display="sites",
                          scales=3))
rda.scores <- cbind(F,rda.scores)
centroids <- aggregate(cbind(RDA1,PC1)~F,rda.scores, mean)
f <- function(z)sd(z)/sqrt(length(z)) # функция для std.err
se <- aggregate(cbind(RDA1,PC1)~F, rda.scores,f)
names(se) <- c("F","RDA1.se","PC1.se")
# объединяем  координаты центроидов и стандартные ошибки 
centroids <- merge(centroids,se, by="F")    

#  Формируем диаграмму
ggplot() +  geom_point(data=rda.scores,
    aes(x=RDA1,y=PC1,shape=F,colour=F),size=2) + 
  xlim(c(-2, 5)) + ylim(c(-0.75, 1)) +
  scale_colour_manual( values = c('purple', 'blue')) +
  stat_ellipse(data=rda.scores,aes(x=RDA1,y=PC1,fill=F),
                    geom="polygon",level=0.8,alpha=0.2)+
  geom_point(data=centroids, aes(x=RDA1,y=PC1,colour=F),
             shape=1,size=4)+
#  geom_errorbar(data=centroids,  aes(ymin = (PC1 - PC1.se),
#         ymax = (PC1 + PC1.se)),width=0.1)+
# Error in eval(expr, envir, enclos) : object 'x' not found
#  geom_errorbarh(data=centroids, aes(xmin=RDA1-RDA1.se,
          xmax=RDA1+RDA1.se),height=0.1) + theme_bw()

#  Многомерный дисперсионный анализ
mod.an <- manova(as.matrix(Y) ~ F)
anova(mod.an)
library(Hotelling) ; split.data <- split(Y,F)  
summ.hot <- hotelling.test(split.data[[1]], split.data[[2]],
                 perm = T, B = 1000)

#-----------------------------------------------------------------------
#  2.6. Методы кластеризации без учителя 
##-----------------------------------------------------------------------

set.seed(32297) 
d <- data.frame(x=runif(100),y=runif(100))
clus <- kmeans(d,centers=5) ;  d$cluster <- clus$cluster
table(d$cluster)
library('ggplot2'); library('grDevices')
h <- do.call(rbind, lapply(unique(clus$cluster),
function(c) { f <- subset(d,cluster==c); f[chull(f),]}))
ggplot() +
geom_text(data=d,aes(label=cluster,x=x,y=y,
                           color=cluster),size=3) +
geom_polygon(data=h,aes(x=x,y=y,group=cluster,
   fill=as.factor(cluster)),alpha=0.4,linetype=0) +
theme(legend.position = "none")

library('reshape2')
n <- dim(d)[[1]] 
pairs <- data.frame(
  ca = as.vector(outer(1:n,1:n,function(a,b) d[a,'cluster'])),
  cb = as.vector(outer(1:n,1:n,function(a,b) d[b,'cluster'])),
  dist = as.vector(outer(1:n,1:n,function(a,b)
  sqrt((d[a,'x']-d[b,'x'])^2 + (d[a,'y']-d[b,'y'])^2)))
)
dcast(pairs,ca~cb,value.var='dist',mean)


