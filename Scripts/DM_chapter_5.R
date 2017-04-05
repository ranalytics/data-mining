#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 5. БИНАРНЫЕ МАТРИЦЫ И АССОЦИАТИВНЫЕ ПРАВИЛА
#########################################################################

#  5.1. Классификация в бинарных пространствах с использованием классических моделей
#-----------------------------------------------------------------------
DFace <- read.delim(file="Faces.txt", header=TRUE, row.names=1)
head(DFace,8)

# Логистическая регрессия
logit <- glm((Class-1) ~ ., data=DFas, family=binomial)
logit.step <- step(logit)
summary(logit.step)
mp <- predict(logit.step, type="response") 
pred = ifelse(mp>0.5,2,1)
CTab <- table(Факт=DFace$Class, Прогноз=pred))
Acc = mean(pred == DGlass$F)
paste("Точность=", round(100*Acc, 2), "%", sep="")
barplot(mp-0.5, col = "steelblue",xlab = "Члены выборки",
       ylab = "Вероятности P - 0.5"))

library(boot)
cv.glm(DFace, logit.step)$delta

# Дискриминантный анализ
library("MASS") 
Fac.lda <- lda(DFace[,1:16], grouping = DFace$Class) 
pred <- predict(Fac.lda, DFace[,1:16])
table(Факт=DFace$Class, Прогноз=pred$class)
Acc = mean(DFace$Class == pred$class)
paste("Точность=", round(100*Acc, 2), "%", sep="") 

# Бинарный Дискриминантный анализ
library(binda)
Xtrain <- as.matrix(DFace[,1:16])
is.binaryMatrix(Xtrain) # Проверяем бинарность матрицы
Class = as.factor(DFace$Class)
binda.fit = binda(Xtrain,Class,lambda.freq=1)
pred <- predict(binda.fit, Xtrain)
table(Факт=DFace$Class, Прогноз=pred$class)
Acc = mean(DFace$Class!= pred$class)
paste("Точность=", round(100*Acc, 2), "%", sep="")

# ранжируем предикторы
binda.x <-  binda.ranking(Xtrain, Class)
plot(binda.x, top=40, arrow.col="blue", zeroaxis.col="red",
ylab="Предикторы", main="")

library(caret) 
binda.train <- train(Xtrain,Class, method = "binda")
binda.train$finalModel 
pred.train <- predict(binda.train, Xtrain)
CTab <- table(Факт=DFace$Class, Прогноз=pred.train)
Acc = mean(DFace$Class == pred.train)
paste("Точность=", round(100*Acc, 2), "%", sep="")

#-----------------------------------------------------------------------
#  5.2. Бинарные деревья решений
#-----------------------------------------------------------------------

#  Набор пользовательских функций для выполнения расчетов
Entropy <- function( vls ) {
  res <- vls/sum(vls)*log2(vls/sum(vls)) ; res[vls == 0] <- 0
  -sum(res)
}
InformationGain <- function( tble ) {
  entropyBefore <- Entropy(colSums(tble))
  s <- rowSums(tble)
  entropyAfter <- sum(s/sum(s)*
         apply(tble, MARGIN = 1, FUN = Entropy ))
  informationGain <- entropyBefore - entropyAfter
  return (informationGain)
}
IsPure <- function(data) {
              length(unique(data[,ncol(data)])) == 1
}
TrainID3 <- function(node, data) {
      node$obsCount <- nrow(data)
  #  если текущий набор данных принадлежит к одному классу, то
  if (IsPure(data)) {
    #создается лист дерева с экземплярами этого класса
    child <- node$AddChild(unique(data[,ncol(data)]))
    node$feature <- tail(names(data), 1)
    child$obsCount <- nrow(data)
    child$feature <- ''
  } else {
    # рассчитывается вектор информационный выигрышей IG
    ig <- sapply(colnames(data)[-ncol(data)], 
            function(x) InformationGain(
              table(data[,x], data[,ncol(data)])))
    #  выбирается значение признака с наибольшей величиной IG
    feature <- names(which.max(ig))
    node$feature <- feature
    # создается подмножества данных на основе этого значения 
    childObs <- split(data[ ,names(data) != feature, 
		drop = FALSE], data[ ,feature], drop = TRUE)
    #  создаются дочерние узлы дерева c именем признака
    for(i in 1:length(childObs)) {
      child <- node$AddChild(names(childObs)[i])
    # осуществляется рекурсия алгоритма для дочернего узла
      TrainID3(child, childObs[[i]])
    }
  }
}
Predict <- function(tree, features) {
  if (tree$children[[1]]$isLeaf) 
             return (tree$children[[1]]$name)
  child <- tree$children[[features[[tree$feature]]]]
  return ( Predict(child, features))
}

#  Выполнение расчетов
DFace <- read.delim(file="Faces.txt",header=TRUE, row.names=1)
DFaceN <-DFace[,-17]; Class <- DFace$Class 
DFaceN[DFaceN==1]<- "Да" ; DFaceN[DFaceN==0] <- "Нет"
Class[Class==1] <- "Патриот" ; Class[Class==2] <- "Демократ"
DFaceN <- cbind(DFaceN,Class) 
colnames(DFaceN) <- c("голова_круглая","уши_оттопырен",
"нос_круглый","глаза_круглые","лоб_морщины",
"носогубн_складка","губы_толстые","волосы","усы","борода",
"очки","родинка_щеке","бабочка","брови_подняты","серьга",
"курит_трубка","Группа")
library(data.tree)
tree <- Node$new("DFaceN")  # Создаем "пустое" дерево
TrainID3(tree, DFaceN)
print(tree, "feature", "obsCount")
pred <- apply(DFaceN[,-17],1, function(x) Predict(tree,x))
table(Факт=DFaceN$Группа, Прогноз=pred)  

#  Скользящий контроль
Nerr <- 0
for (i in 1:nrow(DFaceN)) {
    tree <- Node$new("DFaceN")
    TrainID3(tree, DFaceN[-i,])
    Nerr = Nerr + 
sum(DFaceN[i,17]!=Predict(tree,DFaceN[i,-17]))  }
 (Nerr/nrow(DFaceN)) 

#-----------------------------------------------------------------------
#  5.3. Поиск логических закономерностей в данных 
#-----------------------------------------------------------------------

#  Набор пользовательских функций для выполнения расчетов
#  Преобазование числа в последовательность битов (5 -> "0101")
number2binchar = function(number, nBits) {
   paste(tail(rev(as.numeric(intToBits(number))),nBits),collapse="")}
#  Поиск коньюнкций по набору битовых масок 
MaskCompare = function(Nclass, KSize, BitMask, vec_pos, vec_neg, ColCom) {
  nK <- sapply(BitMask, function(x) {
       if (sum(x==vec_neg)>0) return (0)
       countK = sum(x==vec_pos) ; if (minNum>countK) return (0)
      #  Cохранение конъюнкции  в трех объектах класса list
       Value.list[[length(Value.list)+1]] <<- list(Nclass=Nclass, KSize=KSize, 
               countK=countK, Bits=x)
       ColCom.list[[length(ColCom.list)+1]] <<- list(ColCom)
       RowList.list[[length(RowList.list)+1]] <<- list(which(vec_pos %in% x))
      return (countK) } )
# return (nK[nK>0])
}
 
DFace <- read.delim(file="Faces.txt", header=TRUE, row.names=1)
maxKSize <- 4
minNum <- 4
#  Списки для хранения результатов
Value.list <- list()  # Значения Nclass, KSize, BitMask, countK
ColCom.list <- list()  # Вектора наименований переменных ColCom
RowList.list <- list()  # Номера индексов строк RowList

for (KSize in 2:maxKSize) {
   BitMask <- sapply(0:(2^KSize-1),function(x) number2binchar(x,KSize))
   cols <-  combn(colnames(DFace[,-17]),KSize)
     for (i in 1:ncol(cols))  {
        SubArr <- DFace[,(names(DFace) %in% cols[,i])]
        vec1 <- apply(SubArr[DFace$Class==1,],1,function (x) paste(x,collapse=""))
        vec2 <- apply(SubArr[DFace$Class==2,],1,function (x) paste(x,collapse=""))
        MaskCompare(1, KSize, BitMask, vec1, vec2, cols[,i])
        MaskCompare(2, KSize, BitMask, vec2, vec1, cols[,i])
     }
}
#  Создание результирующей таблицы
DFval = do.call(rbind.data.frame, Value.list)
nrow = length(Value.list)
DFvar <- as.data.frame(matrix(NA, ncol=maxKSize+1, nrow=nrow,
              dimnames = list(1:nrow, c(
              paste("L", 1:maxKSize, sep=""),"Объекты:"))))
for (i in 1:nrow) {
    Varl <- unlist(ColCom.list[[i]])
    DFvar[i, 1:length( Varl)] <- Varl
    Objl <- unlist(RowList.list[[i]])
    DFvar[i, maxKSize+1] <- paste(Objl,collapse=" ")
                 }
DFvar[is.na(DFvar)] <- " "
DFout <- cbind(DFval, DFvar)

#  Вывод результатов
print ("Конъюнкции класса 1") ; DFout[DFout$Nclass==1,]
print ("Конъюнкции класса 2") ; DFout[DFout$Nclass==2,]

#-----------------------------------------------------------------------
#  5.4.  Алгоритмы выделения ассоциативных правил
#-----------------------------------------------------------------------
# Переформирование исходных данных
DFace <- read.delim(file="Faces.txt", 
                          header=TRUE, row.names=1)
Class <- DFace$Class ; DFaceN <-DFace[,-17]
colnames(DFaceN) <- c("голова_круглая","уши_оттопырен",
   "нос_круглый","глаза_круглые","лоб_морщины",
   "носогубн_складка","губы_толстые","волосы","усы","борода",
   "очки","родинка_щеке","бабочка", "брови_подняты","серьга",
   "курит_трубка")
Class[Class==1] <- "Патриот" ; Class[Class==2] <- "Демократ"
items_list <- sapply(1:nrow(DFaceN),function(i) 
paste(c(Class[i], colnames(DFaceN[i,DFaceN[i,]==1])),
            collapse=",",sep="\n"))
head(items_list)
write(items_list, file = "face_basket.txt")

library("arules")
trans = read.transactions("face_basket.txt", 
                    format = "basket", sep=",")
inspect(trans)   #  Выводимые данные не показаны
summary(trans)   #  Выводимые данные показаны частично
image(trans)
itemFrequencyPlot(trans, support = 0.1, cex.names=0.8)
rules <- apriori(trans,
    parameter = list(support = 0.01, confidence = 0.6))
summary(rules)  #  Выводимые данные показаны частично
library("arulesViz")
plot(rules, measure=c("support","lift"), shading="confidence")
# Выделение подмножеств правил
rulesPat <- subset(rules, subset = rhs %in% "Патриот" & 
                           lift > 1.8)
inspect(head(rulesPat, n = 10, by = "support"))
plot(head(sort(rulesPat, by="support"), 10),
         method="paracoord")
rulesDem <- subset(rules, subset = rhs %in% "Демократ" &
                    lift > 1.8)
inspect(head(rulesDem, n = 10, by = "support"))
plot (head(sort(rulesDem, by="support"), 10), method="graph",
        control = list(nodeCol = grey.colors(10), 
        edgeCol = grey(.7), alpha = 1))

#-----------------------------------------------------------------------
# 5.5. Анализ последовательностей знаков или событий
#-----------------------------------------------------------------------
library(TraMineR)
data(biofam)
summary(biofam) #  Результаты не приведены
biofam.lab <- c("Родит", "Отд", "Сем+Род", 
        "Сем.", "Один+Род", "Один+Отд", "Сем+Дет", "Развод") 
biofam.seq <- seqdef(biofam[, 10:25], labels = biofam.lab)
par(mfrow = c(2, 2)); seqiplot(biofam.seq, withlegend = FALSE)
seqdplot(biofam.seq, withlegend = FALSE)
seqfplot(biofam.seq, withlegend = FALSE)
seqlegend(biofam.seq)

seqdplot(biofam.seq,group=biofam$sex, withlegend = FALSE)
seqHtplot(biofam.seq)
Entropy <- seqient(biofam.seq)
ageg = cut(biofam$birthy, c(1909,1918,1928,1938,1948,1958),
label = c("1909-18","1919-28","1929-38","1939-48","1949-58"),
           include.lowest = TRUE)
boxplot(Entropy ~ ageg, data = biofam, 
    xlab = "Диапазон годов рождения", 
    ylab = "Энтропия последовательностей", col = "cyan")
Turbulence <- seqST(biofam.seq)
plot(Turbulence, Entropy, xlab = "Турбулентность", 
                          ylab = "Энтропия")
m.turb <- lm(Turbulence ~ sex + birthyr, data = biofam)  
summary(m.turb)
couts <- seqsubm(biofam.seq, method = "TRATE")
biofam.om <- seqdist(biofam.seq, method = "OM", 
                     indel = 3, sm = couts)
round(biofam.om[1:10, 1:10], 1)

library(cluster)
clusterward <- agnes(biofam.om, diss = TRUE, method = "ward")
plot(clusterward, which.plots = 2)

# Распилим дерево на три части 
cluster3 <- cutree(clusterward, k = 3) 
cluster3 <- factor(cluster3, 
        labels = c("Кластер 1", "Кластер 2", "Кластер 3"))
table(cluster3)
par(mfrow = c(2,2))
seqmtplot(biofam.seq, group = cluster3)
dissvar(biofam.om) # Общая дисперсия попарных расстояний
da <- dissassoc(biofam.om, group=biofam$sex, R = 1000) 
print(da)
dlm <- dissmfac(biofam.om ~ sex + ageg, data =biofam, R = 100)
disstree (biofam.om ~ sex +  birthyr, data =biofam, R = 100)

