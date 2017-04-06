#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 8. МОДЕЛИРОВАНИЕ ПОРЯДКОВЫХ И СЧЕТНЫХ ПЕРЕМЕННЫХ
#########################################################################

#  8.1. Модель логита для порядковой переменной
#-----------------------------------------------------------------------
# abalone <- read.csv("http://archive.ics.uci.edu/ml/
# 		machine-learning-databases/abalone/abalone.data",
# 		header=FALSE)
abalone <- read.csv("abalone.data", header = FALSE)
names(abalone) <- c("пол", "длина", "диаметр", "высота", "вес.общ", "вес.тела", "вес.внут", "вес.рак", "rings")
summary(abalone[,c(1,2,9)]) 

require(ggplot2) ;  library(ggcorrplot)
# Отображение корреляционной матрицы (рис. 8.2)
M <- cor(abalone[,2:8])
ggcorrplot(M, hc.order = TRUE, type = "lower",
  colors = c("white","yellow","purple" ),    lab = TRUE)
# Отображение линейной зависимости (рис. 8.3)

ggplot(abalone) + aes(длина, rings, color = пол) +
    geom_point() + labs(x = "Длина раковины",
              y = "Число колец", color = "Пол") + 
    stat_smooth(method = "lm", se = FALSE, size = 2)

ggplot(abalone) + aes(rings, fill = пол) + 
     geom_density(position = "stack")+ 
     geom_vline(xintercept = quantile(abalone$rings, 
     p = c(0.25, 0.5, 0.75)),colour = "blue", 
     linetype = 5, size = 1.5)
table (cut(abalone$rings,breaks=quantile(abalone$rings,
     c(0,.25,.50,.75,1),include.lowest=TRUE)))

abalone$Возраст <- cut(abalone$rings,breaks=c(0,7,9,11,29),
  labels=c("Q1","Q2","Q3","Q4"),include.lowest=TRUE)
save(abalone,file="abalone.RData")
table (abalone$Возраст )

library(MASS)
CL.aq <- polr(Возраст ~ .,  data = abalone[,-9])
summary(CL.aq, digits = 3)
# Оценка доверительных интервалов для коэффициентов
confint.default(CL.aq)

CLs.aq <- stepAIC (CL.aq)
summary(CLs.aq, digits = 3) 
confint.default(CLs.aq) 
CL0.aq <- polr(Возраст ~ 1,  data = abalone[,-9])
anova(CL0.aq, CLs.aq)

Probs <- fitted(CLs.aq)
head(Probs,3)

#  Построение таблицы сопряженности "Факт-Прогноз"
pred=apply(Probs,1,function(x) colnames(Probs)[which(x==max(x))])
table(Факт=abalone$Возраст,Прогноз=pred)
Acc = mean(pred == abalone$Возраст)
paste("Точность=", round(100*Acc, 2), "%", sep="")

#   Подготовка данных для графика
VM <- apply(abalone[,2:8],2, mean)
d.plot <- as.data.frame(matrix(VM, ncol=7, nrow=50,
 byrow=TRUE, dimnames = list(1:50,names(VM))))
d.plot <- cbind(пол=rep("F",50), d.plot)
d.plot$диаметр <- seq(min(abalone$диаметр),
         max(abalone$диаметр), len = 50)
d.pplot <- cbind(d.plot, predict(CL.aq, newdata = d.plot,
         type = "probs", se = TRUE))
#   Прорисовка компонент графика
plot (1,1, xlim=c(0,max(abalone$диаметр)),ylim=c(0,0.8),
 	type='n', xlab="Диаметр раковины", ylab="Вероятность Р")
lines(d.pplot[,c(3,9)],lwd=2, col="green")
lines(d.pplot[,c(3,10)],lwd=2, col="blue")
lines(d.pplot[,c(3,11)],lwd=2, col="black")
lines(d.pplot[,c(3,12)],lwd=2, col="red")
legend("topright", c("Q1","Q2","Q3","Q4"), lwd=2, 
       col=c(3,4,1,2))

#-----------------------------------------------------------------------
#  8.2. Настройка параметров нейронных сетей в пакте caret
#-----------------------------------------------------------------------
library(nnet); library(caret)
set.seed(123) ; load (file="abalone.RData")
train.aba <- train(Возраст ~ ., data = abalone[,c(3:8,10)],
method = "nnet", trace = F, linout = 1,
tuneGrid = expand.grid(.decay = c(0,0.05,0.2), .size = 4:9),
trControl = trainControl(method = "cv")) 

source("nnet_plot_update.r")
# plot.nnet(train.aba)  - можно все извлечь из train-объекта
# nn.aba <- train.aba$finalModel
nn.aba <- nnet(Возраст ~ .,  data = abalone[,c(3:8,10)],
decay =0,size = 7, niter=200)
plot.nnet(nn.aba)
summary(nn.aba) 
pred = predict(nn.aba,   abalone[,3:8],   type="class")
nn.table = table(abalone[,10],  pred)
confusionMatrix(nn.table) 

pcaNNet.Fit <- pcaNNet(abalone[,3:8], abalone[,10], size = 7, thresh = 0.975, linout = TRUE, trace = FALSE)
pred <- predict(pcaNNet.Fit ,abalone[,3:8],   type="class")
table(Факт=abalone$Возраст,Прогноз=pred)
Acc = mean(pred == abalone$Возраст)
paste("Точность=", round(100*Acc, 2), "%", sep="")

avNNet.Fit <- avNNet(Возраст ~ ., data = abalone[,c(3:8,10)],
    size = 7, repeats = 10, linout = TRUE, 
    trace = FALSE, bag = TRUE)
pred <- predict(avNNet.Fit ,abalone[,3:8],   type="class")
Acc = mean(pred == abalone$Возраст)
paste("Точность=", round(100*Acc, 2), "%", sep="")

#-----------------------------------------------------------------------
#  8.3. Методы комплексации модельных прогнозов
#-----------------------------------------------------------------------
load(file="abalone.RData")
# Формирование матрицы предикторов
X <- model.matrix(rings ~ .,data=abalone[,1:9])[,-1]
Y <- abalone$rings
# Разделение на обучающую и тестовую выборки
train <- runif(nrow(X)) <= .66

library(caret); library(party)
myControl <- trainControl(method='cv', number=10, 
  savePredictions=TRUE, returnData=FALSE, verboseIter=TRUE)
PP <- c('center', 'scale')
# Обучение избранных моделей
# Регрессия на k-ближайших соседей 
m.knn <- train(X[train,], Y[train], method='knn',
              trControl=myControl, preProcess=PP)
# Линейная модель
m.lm <- train(X[train,], Y[train], method='lm',
              trControl=myControl, preProcess=PP)
# Гребневая регрессия с регуляризацией
m.rlm <- train(X[train,], Y[train], method='glmnet', trControl=myControl, preProcess=PP)
# Модель опорных векторов
m.svm <- train(X[train,], Y[train], method='svmRadial', trControl=myControl, preProcess=PP)
#  Метод случайного леса
m.rf <- train(X[train,], Y[train], method='rf', trControl=myControl)
# Бэггинг деревьев условного вывода
m.ctr <- train(X[train,], Y[train], "bag",
    B = 10, bagControl = bagControl(fit = ctreeBag$fit,
    predict = ctreeBag$pred, aggregate = ctreeBag$aggregate)) 

# Создание списка всех моделей
all.models <- list(m.knn, m.lm, m.rlm, m.svm, m.rf, m.ctr)
names(all.models) <- sapply(all.models, function(x) x$method)
sort(sapply(all.models, function(x) min(x$results$RMSE)))

preds.all <- data.frame(sapply(all.models,
                   function(x){predict(x, X)}))
head(preds.all)
rmse <- function(x,y){sqrt(mean((x - y)^2))}
print("Среднеквадратичная ошибка на тестовой выборке")
print(sort(apply(preds.all[!train,], 2, rmse, y = Y[!train])), digits= 3) 

library(ForecastCombinations)
scheme= c("simple",  "variance based", "ols", "robust",
         "cls", "best")
combine_f <- list()
for (i in scheme) {
   combine_f[[i]] <- Forecast_comb(obs = Y[train],
        fhat = as.matrix(preds.all[train, ]),
        fhat_new = as.matrix(preds.all[!train, ]),
        Averaging_scheme = i)
    tmp <- round(sqrt(mean((combine_f[[i]]$pred -
           Y[!train])^2)), 3)
    cat(i, ":", tmp, "\n")
 }

w.list <- sapply(combine_f, function(sp.list) sp.list$weights)
weights <- as.data.frame.list(w.list)
rownames(weights) <- scheme  
wstan <- function (x) abs(x)/sum(abs(x)) 
weights[3,] <- wstan(weights[3,])
weights[4,] <- wstan(weights[4,])
print(round(weights,3) 

require(nnls)
m.nnls <- nnls(as.matrix(preds.all[train, ]), Y[train])
coef(m.nnls)

combine_f_all <- Forecast_comb_all(Y[train],
       fhat = as.matrix(preds.all[train, ]),
      fhat_new = as.matrix(preds.all[!train, ]))
# Усредняем комбинированные прогнозы по всем регрессиям
Combined_f_all_simple <- apply(combine_f_all$pred, 1, mean)
print(sqrt(mean((Combined_f_all_simple - Y[!train])^2)),
      digits= 3 )

# Комбинирование прогнозов с использованием Cp Mallow:
Combined_f_all_mal <- 
         t( combine_f_all$mal %*% t(combine_f_all$pred) )
print(sqrt(mean((Combined_f_all_mal - Y[!train])^2)),
      digits= 3)

#-----------------------------------------------------------------------
#  8.4. Обобщенные линейные модели на основе регрессии Пуассона
#-----------------------------------------------------------------------
RK <- read.table(file = "RoadKills.txt", 
                      header = TRUE, dec = ".")
y <- RK$TOT.N
library(vcd)   ##  Visualizing Categorical Data
gf<-goodfit(table(y),type= "poisson",method= "ML") 
#  Визуализация подогнанных данных гистограммой 
plot(gf,ylab="Частоты", xlab="Число классов")

M1 <-  glm(TOT.N ~ D.PARK, family =  poisson, data =  RK)
summary(M1)
require(ggplot2)
MyData=data.frame(D.PARK=seq(from=0,to=25000,by=1000))
G<-predict(M5,newdata=MyData,type="link",se=T)
MyData$F<-exp(G$fit)
MyData$FSEUP<-exp(G$fit+1.96*G$se.fit)
MyData$FSELOW<-exp(G$fit-1.96*G$se.fit)
ggplot(MyData, aes(D.PARK, F)) +
  geom_ribbon(aes(ymin = FSEUP, ymax = FSELOW), 
  fill = 3, alpha = .25) +   geom_line(colour = 3) +
  labs(x = "Расстояние до парка", y = "Убийства на дорогах")+
  geom_point(data=RK, aes(D.PARK, TOT.N) )

RK[, 7:14] <- sqrt(RK[, 7:14])
library(car)
sort(vif(glm(TOT.N ~ ., data=RK[,c(5, 7:23)], , family =  poisson)))
RK.11 <- cbind(RK$D.WAT.RES, RK$D.WAT.COUR, RK$WAT.RES,
     RK$D.PARK, RK$L.D.ROAD, RK$L.P.ROAD, RK$MONT.S, RK$SHRUB,
     RK$L.WAT.C, RK$POLIC, RK$SQ.URBAN)
Y <- RK$TOT.N
M2 <-  glm(Y ~ ., family =  poisson, data =  RK.11)
summary(M2)
M3 <- step(M2)
summary(M3)

anova(M2, M3, test = "Chi")
anova(M2, M3, test = "Chi")
M4 <- glm(Y ~ D.WAT.RES + WAT.RES + D.PARK + L.P.ROAD +
           MONT.S + SHRUB + L.WAT.C + POLIC + URBAN, 
           family =  quasipoisson, data =  RK.11)
summary(M4)
drop1(M4, test = "Chi")

library(caret)
glmFuncs <- lmFuncs
glmFuncs$fit <- function (x, y, first, last, ...) { 
  tmp <- as.data.frame(x)
  tmp$y <- y 
  glm(y ~ ., data = tmp, family=quasipoisson(link='log'))
}
set.seed(13)
ctrl <- rfeControl(functions = glmFuncs, method = "cv", 
           verbose = FALSE, returnResamp = "final")
lmProfileF <- rfe(x = RK.11, y = RK$TOT.N, sizes = 1:10,
       rfeControl = ctrl) 
train(TOT.N ~ D.PARK, data=RK, method='glm', 
   family = quasipoisson, trControl = 
   trainControl(method = "cv")) 
   
train(TOT.N ~ D.WAT.RES + WAT.RES + D.PARK + L.P.ROAD +
    MONT.S + SHRUB + L.WAT.C + POLIC + URBAN, data=RK,
    method='glm', family = quasipoisson,
   trControl = trainControl(method = "cv"))

library(MASS)
M5 <-  glm.nb(Y ~ .,  data =  RK.11)
summary(M5)
M6 <- step(M5)
summary(M6)
train(TOT.N ~ D.PARK, data=RK, method='glm.nb', 
   trControl = trainControl(method = "cv"))

train(TOT.N ~ D.WAT.RES + WAT.RES + D.PARK + L.P.ROAD + 
    MONT.S + SHRUB + L.WAT.C, data = RK, method='glm.nb', 
   trControl = trainControl(method = "cv"))

#-----------------------------------------------------------------------
# 8.5. ZIP- и барьерные модели счетных данных 
#-----------------------------------------------------------------------
library(AER)
data("NMES1988")        # Загружаем данные и отбираем столбцы
nmes <- NMES1988[, c(1, 6:8, 13, 15, 18)]  # переменных
plot(table(nmes$visits))  
sum(nmes$visits < 1)   # наблюдаемое число нулей
mod1 <- glm(visits ~ ., data = nmes, family = "poisson")
mu <- predict(mod1, type = "response") # среднее Пуассона
exp <- sum(dpois(x = 0, lambda = mu))  # теоретическая частота
round(exp)                             # нулевых значений

library(pscl)
M.ZIP <- zeroinfl(visits ~ ., data = nmes, 
               dist = "poisson", link="logit")
summary(M.ZIP)
library(countreg)
rootogram(M.ZIP)

M.ZINB <- zeroinfl(visits ~ ., data = nmes, 
                 dist = "negbin", link="logit")
summary(M.ZINB)
rootogram(M.ZINB)

M.ZAP <-  hurdle(visits ~ ., data = nmes, 
    dist = "poisson", zero.dist = "binomial", link = "logit")
summary(M.ZAP)
rootogram(M.ZAP)

M.ZINB <-  hurdle (visits ~ ., data = nmes, 
                 dist = "negbin", link="logit")
summary(M.ZINB)
rootogram(M.ZINB)

fm <- list("ZIP" = M.ZIP, "ZINB" = M.ZINB, 
      "Hurdle-Pois" = M.ZAP,  "Hurdle-NB" = M.ZANB)
rbind(logLik = sapply(fm, function(x) round(logLik(x), 
                       digits = 0)),
     AIC = sapply(fm, function(x) round(AIC(x), digits = 0)))
