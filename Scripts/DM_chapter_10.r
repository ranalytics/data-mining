#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 10. КЛАСТЕРНЫЙ АНАЛИЗ
#########################################################################

#  10.1. Алгоритмы кластеризации, основанные на разделении
#-----------------------------------------------------------------------
library(cluster)
data("USArrests")
df.stand <- as.data.frame(scale(USArrests))
C(kmeans(df.stand,centers=5, nstart=1)$tot.withinss,
  kmeans(df.stand,centers=5, nstart=25)$tot.withinss) 

k.max <- 15 # максимальное число кластеров
wss <- sapply(1:k.max, function(k)
     {kmeans(df.stand, k, nstart=10 )$tot.withinss})
plot(1:k.max, wss, type="b", pch = 19, frame = FALSE, 
       xlab="Число кластеров K", 
       ylab="Общая внутригрупповая сумма квадратов")
# Формируем график с помощью fviz_nbclust()
library(factoextra)
fviz_nbclust(df.stand, kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)

set.seed(123)
gap_stat <- clusGap(df.stand, FUN = kmeans, nstart = 10,
                    K.max = 10, B = 50)
# Печать и визуализация результатов
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

set.seed(123)
gap_stat <- clusGap(df.stand, FUN =  pam, K.max = 7, B = 100)
print(gap_stat, method = "firstmax")
k.pam <- pam(df.stand, k=4)

fviz_silhouette(silhouette(k.pam))

fviz_nbclust(df.stand, pam, method = "silhouette")

ShState <- read.delim("St_USA.txt")
rownames(df.stand) <- ShState$Short
fviz_cluster(pam(df.stand, 4), stand=FALSE)
fviz_cluster(clara(df.stand, 4), stand=FALSE, frame.type = "t", frame.level = 0.7)


#-----------------------------------------------------------------------
#  10.2. Иерархическая кластеризация 
#-----------------------------------------------------------------------

library(cluster)
data("USArrests")
d <- dist(scale(USArrests), method = "euclidean")
# Параметр hang=-1  выравнивает метки
plot(hclust(d, method = "average" ), cex = 0.7, hang=-1)

plot(hclust(d, method = "single" ), cex = 0.7)

res.hc <- hclust(d, method = "complete" )
grp <- cutree(res.hc, k = 4)  # Разрезание дерева на 4 группы
plot(res.hc, cex = 0.7)
rect.hclust(res.hc, k = 4, border = 2:5)

hcd <- as.dendrogram(hclust(d, method = "ward.D2" ))
nodePar <- list(lab.cex = 0.7, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
plot(hcd,  xlab = "Height", nodePar = nodePar, horiz = TRUE,
                edgePar = list(col = 2:3, lwd = 2:1))

library(dendextend)
# Вычисляем 2 иерархические кластеризации
hc1 <- hclust(d, method = "average")
hc2 <- hclust(d, method = "ward.D2")
# Создаем две дендрограммы
dend1 <- as.dendrogram (hc1); dend2 <- as.dendrogram (hc2)
# На диаграмме типа "клубок" показываем цветом общие ветви
tanglegram(dend1, dend2,
  common_subtrees_color_branches = TRUE)

c(cor_cophenetic(dend1, dend2), # кофенетическая корреляция
cor_bakers_gamma(dend1, dend2))  # корреляция Бейкера

# Создаем множество дендрограмм для сравнения
dend1 <- d %>% hclust("com") %>% as.dendrogram
dend2 <- d %>% hclust("single") %>% as.dendrogram
dend3 <- d %>% hclust("ave") %>% as.dendrogram
dend4 <- d %>% hclust("centroid") %>% as.dendrogram
dend5 <- d %>% hclust("ward.D2") %>% as.dendrogram
# Создаем список дендрограмм и корреляционную матрицу
dend_list <- dendlist("Complete" = dend1, "Single" = dend2,
     "Average" = dend3, "Centroid" = dend4, "Ward.D2" = dend5)
cors <- cor.dendlist(dend_list)
round(cors, 2)  
library(corrplot) # изображение корреляционной матрицы
corrplot(cors, "pie", "lower")

#-----------------------------------------------------------------------
#  10.3. Оценка качества кластеризации 
#-----------------------------------------------------------------------

library(cluster)  ; data("USArrests")
df.stand <- as.data.frame(scale(USArrests))
library("factoextra")
get_clust_tendency(df.stand, n = 30,
     gradient = list(low = "steelblue", high = "white"))

library(NbClust)
nb <- NbClust(df.stand, distance = "euclidean", min.nc = 2,
        max.nc = 8, method = "average", index ="all")
nb$Best.nc
fviz_nbclust(nb) + theme_minimal()

d <- dist(df.stand, method = "euclidean")
library(vegan)
hc_list <- list(hc1 <- hclust(d,"com"),
hc2 <-  hclust(d,"single"), hc3 <-  hclust(d,"ave"),
hc4 <-  hclust(d, "centroid"),hc5 <-   hclust(d, "ward.D2"))
Coph <- rbind(
MantelStat <- unlist(lapply(hc_list, 
     function (hc) mantel(d, cophenetic(hc))$statistic)),
MantelP <- unlist(lapply(hc_list, 
     function (hc) mantel(d, cophenetic(hc))$signif)))
colnames(Coph) <- c("Complete", "Single","Average",
     "Centroid","Ward.D2") 
rownames(Coph) <- c("W Мантеля","Р-значение")
round(Coph, 3)

data(Boston, package = "MASS")
library(pvclust)
set.seed(123)
#  Бутстреп деревьев и расчет BP- и AU- вероятностей для узлов
boston.pv <- pvclust(Boston, nboot=100, method.dist="cor", 
                  method.hclust="average")
plot(boston.pv)  # дендрограмма с p-значениями
pvrect(boston.pv) # выделение боксами достоверных фрагментов


#-----------------------------------------------------------------------
#  10.4. Другие алгоритмы кластеризации  
#-----------------------------------------------------------------------

#-  Иерархическая кластеризация на главные компоненты
library(FactoMineR)
data(Boston, package = "MASS")
df.scale <- scale(Boston)
res.pca <- PCA(df.scale, ncp = 5, graph=TRUE)
get_eig(res.pca)

res.hcpc <- HCPC(res.pca, graph = TRUE)
plot(res.hcpc, choice ="tree")
plot(res.hcpc, choice = "3D.map", ind.names=FALSE)

#-  Метод нечетких k-средних (fuzzy analysis clustering)
library(cluster) 
data("USArrests")
set.seed(123)
res.fanny <- fanny(USArrests, k = 4, memb.exp = 1.7, 
      metric = "euclidean", stand = TRUE, maxit = 500)
print(head(res.fanny$membership),3)
res.fanny$coeff

# Визуализация с использованием corrplot
library(corrplot)
Dunn <- res.fanny$membership^2
corrplot(Dunn[rev(order(rowSums(Dunn))),], is.corr = FALSE)
# Ординационная диаграмма
library(factoextra)
fviz_cluster(res.fanny, frame.type = "norm",
             frame.level = 0.7)

# Статистическая модель кластеризации
data("faithful")
head(faithful, 3)
library("ggplot2")
ggplot(faithful, aes(x=eruptions, y=waiting)) +
  geom_point() + geom_density2d() 

mc <- Mclust(faithful)
summary(mc) ; head(mc$z)

plot(mc, "classification")
plot(mc, "uncertainty")

plot(mc, "BIC")


#-----------------------------------------------------------------------
#  10.5. Самоорганизующиеся карты Кохонена   
#-----------------------------------------------------------------------

data(Boston, package = "MASS")
VarName = c("indus", "dis", "nox",  "medv", "lstat", "age", "rad") 
# отбор переменных для обучения SOM
data_train <- Boston[, VarName]
data_train_matrix <- as.matrix(scale(data_train))
set.seed(123)
som_grid <- somgrid(xdim = 9, ydim=6, topo="hexagonal") 
som_model <- som(data_train_matrix, 
     grid=som_grid, rlen=100, alpha=c(0.05,0.01), 
     keep.data = TRUE)
plot(som_model, type = "changes")

# Зададим палитру цветов
coolBlueHotRed <- function(n, alpha = 1) {
        rainbow(n, end=4/6, alpha=alpha)[n:1] }
# Сколько объектов связано с каждым узлом?
plot(som_model, type = "counts", palette.name=coolBlueHotRed)
# Каково среднее расстояние объектов узла до его прототипов?
plot(som_model, type = "quality", palette.name=coolBlueHotRed)

colB <- ifelse(Boston$black <= 100, "red", "gray70")
plot(som_model, type = "mapping", col =colB, pch = 16)
plot(som_model, type = "codes")

plot(som_model, type = "property", 
        property = som_model$codes[[1]][,1], 
        main = "indus - доля домов, продаваемых в розницу",
        palette.name=coolBlueHotRed)
var_unscaled <- aggregate(as.numeric(data_train[,3]),
 by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2]
plot(som_model, type = "property", property=var_unscaled,
        main="nox - содержание окислов азота",
        palette.name=coolBlueHotRed) 

## Формируем матрицу "узлы   переменные"
mydata <- as.matrix(som_model$codes[[1]])
# Используем иерархическую кластеризацию с порогом при k=5
som_cluster <- cutree(hclust(dist(mydata)), 5)
# Определяем палитру цветов
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c',
           '#d62728', '#9467bd', '#8c564b', '#e377c2')
# Показываем разными цветами кластеры узлов и переменные
plot(som_model, type="codes", 
                bgcol = pretty_palette[som_cluster])
add.cluster.boundaries(som_model, som_cluster) 
