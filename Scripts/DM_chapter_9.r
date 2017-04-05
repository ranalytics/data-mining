#########################################################################
# Шитиков В.К., Мастицкий С.Э. (2017) Классификация, регрессия и другие алгоритмы Data Mining 
# с использованием R. (Адрес доступа: http://www.ievbras.ru/ecostat/Kiril/R/DM )
#########################################################################
#########################################################################
# Глава 9. МЕТОДЫ МНОГОМЕРНОЙ ОРДИНАЦИИ
#########################################################################

# 9.2.  Непараметрический дисперсионный анализ матриц дистанции
#-----------------------------------------------------------------------

library(vegan)
data(dune)
data(dune.env)
results <- adonis(dune ~ Management*A1, data=dune.env, 
                   permutations=499) 

#  Предварительно рассчитываем матрицу расстояний Брея-Кёртиса D <- vegdist(dune)  
mod <- betadisper(D, dune.env$Management,  type = "centroid")
plot(mod, hull= FALSE, main="") ; legend("topright",
             levels(dune.env$Management), pch=1:4, col=1:4)
# Выполнение рандомизационного теста и парных сравнений
permutest(mod, pairwise = TRUE) 

#   Тест Тьюки HSD и график доверительных интервалов
(mod.HSD <- TukeyHSD(mod)) ;  plot(mod.HSD)

#-----------------------------------------------------------------------
#  9.3. Методы ординации объектов и переменных: сравнение и визуализация биплотов
#-----------------------------------------------------------------------

library(labdsv)
data(bryceveg) ;  data(brycesite)
#  Рассчитываем три матрицы расстояний по разным формулам
dis.mht <- dist(bryceveg,"manhattan")
dis.euc <- dist(bryceveg,'euclidean')
dis.bin <- dist(bryceveg,"binary")
#  Создаем объекты cmdscale
mht.pco <- cmdscale(dis.mht, k=10,eig=TRUE)
euc.pco <- cmdscale(dis.euc, k=10,eig=TRUE)
bin.pco <- cmdscale(dis.bin, k=10,eig=TRUE)
# График изменения относительных собственных чисел
plot(1:10, (mht.pco$eig/sum(mht.pco$eig))[1:10],
     type="b",xlab="Число координат",
     ylab="Относительные собственные значения")
lines((euc.pco$eig/sum(euc.pco$eig))[1:10],type="b",col=2)
lines((bin.pco$eig/sum(bin.pco$eig))[1:10],type="b",col=3)
legend ("topright",c("Манхеттен","Евклид", "Хемминг"),
      lwd=1, col=1:3)

plot(euc.pco$points[,1:2], col = 4, pch=17, xlab = "PCO1",
        ylab = "PCO2", main="Расстояние Евклида")
plot(0-bin.pco$points[,1:2], col = 3,pch=17, xlab = "PCO1",
        ylab = "PCO2", main="Расстояние Хемминга")
CCor <-c(cor(dis.euc, dist(euc.pco$points[,1:2],'euclidean')),
cor(dis.mht, dist(mht.pco$points[,1:2],'manhattan')))
names(CCor) <- c("Евклидово","Манхеттенское") 
print(CCor)

data(brycesite)
H <- brycesite$elev  
mht.pco <- pco(dis.mht,k=2) # аналог функции cmdscale()
source("surf_pso.r")
plot(mht.pco, pch=16, col="red")
surf(mht.pco,H, col="blue")
ordtest(mht.pco,H)$p

library(vegan)
ca <- cca(bryceveg)
pro.comp <- procrustes(mht.pco, ca)
plot(pro.comp)

library(FactoMineR)
library("factoextra")
data("decathlon")
colnames(decathlon) <- c("100м","Длина.прыжок","Ядро",
   "Высота.прыжок","400м","110м.барьер","Диск","Шест.прыжок",
   "Копье","1500м",   "Место","Очки", "Соревнование" ) 
head(res.pca$eig)

res.pca$var$contrib
fviz_pca_var(res.pca, col.var="contrib",
   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Раскрасим также текст

res.pca <- PCA(decathlon, quanti.sup = 11:12, quali.sup = 13,
                graph=FALSE)
fviz_pca_ind(res.pca, geom = "point", habillage = 13,
  addEllipses =TRUE, ellipse.level = 0.68) +
  scale_color_brewer(palette="Dark2") +  theme_minimal()


#-----------------------------------------------------------------------
# 9.4. Оценка связи ординации с внешними факторами 
#-----------------------------------------------------------------------

library(vegan)
data(varechem) ; data(varespec)
B <- rev(sort(colSums(varespec)))
BF <- as.matrix(t(round(B/min(B))))
rad <- rad.zipf(BF)
barplot(BF)
lines(1:ncol(BF), rad$fitted, col="green", ,lwd=2)

varespec.chi <- decostand(varespec, method="chi.square")
mod1.cca <- cca(varespec.chi ~ ., varechem)
spec <- as.data.frame(mod1.cca$CCA$v[,1:2])  # Виды
sites <- as.data.frame(mod1.cca$CCA$u[,1:2])  # Площадки
fact <-  as.data.frame(mod1.cca$CCA$biplot[,1:2])  # Факторы
library(ggplot2)
ggplot() + 
  geom_point(data=sites,aes(x=CCA1,y=CCA2),size=3,shape=17,
   color="red") + # добавляем точки площадок
  geom_text(data=spec,aes(x=CCA1,y=CCA2,
   label = abbreviate(rownames(spec),4)),colour="darkgreen",
   fontface = 3) +  # добавляем метки видов
  geom_segment(data=fact,aes(x=0,xend=CCA1,y=0,yend=CCA2),
      arrow = arrow(length = unit(0.5, "cm")),
      colour="darkblue") +  # добавляем стрелки факторов
  geom_text(data = fact, aes(x = CCA1, y = CCA2, 
      label = rownames(fact), size = 6), hjust = 1.2, 
       fontface = 2) + theme_bw()

library(car)
sort(vif(mod1.cca))
paste("AIC=", extractAIC(mod1.cca)[2])

anova(mod1.cca)
anova(mod1.cca, by = "term", step=200)

mod0.cca <- cca(varespec.chi ~ 1, varechem)
mod.cca <- step(mod0.cca, scope = formula(mod1.cca), 
               test = "perm", trace=0)
mod.cca$anova
anova(mod.cca) ; anova(mod.cca, by="axis")


#-----------------------------------------------------------------------
# 9.5. Неметрическое многомерное шкалирование и построение распределения чувствительности видов  
#-----------------------------------------------------------------------

library(vegan)
data(varechem) ; data(varespec)
varespec.chi <- decostand(varespec, method="chi.square")
# Коэффициенты корреляции Спирмена для различных метрик 
rankindex(varechem, varespec.chi,
                c("euc","man","bray","jac","kul"))

vare.mds <- metaMDS(varespec,distance = "kul", trace = FALSE) 
ord.mds <- MDSrotate(vare.mds, varechem$Al)
ef <- envfit(vare.mds, varechem, permu = 999)

plot(vare.mds,  type="n") 
 ordipointlabel(vare.mds, 
         display = "sp", font=4, col="darkgreen", add=TRUE)
abline(h = 0, lty = 3) ;   abline(v = 0, lty = 3)
with(varechem, ordisurf(vare.mds, Al, add = TRUE, col =2))
with(varechem, ordisurf(vare.mds, Ca, add = TRUE, col =4))

Al.log <-log(varechem$Al)
fit.Al = ordisurf(vare.mds, Al.log)
a <- as.data.frame(vare.mds$species)
Sp.Al <- cbind(a, 
         calibrate(fit.Al, newdata = a, se.fit = TRUE))
Sp.Al$val = exp(Sp.Al$fit)
tc <- qt(0.975, nrow(Sp.Al)-1)
Sp.Al$valC = exp(Sp.Al$fit)+ exp(Sp.Al$se.fit)*tc
head(Sp.Al)


# Сортировка по возрастанию val 
df <- Sp.Al[,5:6]
df <- df[order(df$val), ]
df$frac <- ppoints(df$val, 0.5)
head(df)


library(fitdistrplus)
fit <-fitdistr(df$valC, "lognormal")
ep1s=fit$estimate[1]; ep2s=fit$estimate[2]
hcs <- qlnorm(c(0.05, 0.1, 0.2), meanlog=ep1s,sdlog=ep2s)

# 1. Функция для нахождения p-квантили случайной выборки из 
#  логнормального распределения с параметрами  fit
myboot <- function(fit, p){
  xr <- rlnorm(fit$n, meanlog = fit$estimate[1], 
                      sdlog = fit$estimate[2])
  fitr <- fitdist(xr, 'lnorm')
    hc5r <- qlnorm(p, meanlog = fitr$estimate[1], 
                    sdlog = fitr$estimate[2])
  return(hc5r)
}
# 2. Функция, возвращающая значения вероятностей для случайной
#  выборки из логнормального распределения с параметрами  fit
myboot2 <- function(fit, newxs){
  xr <- rlnorm(fit$n, meanlog = fit$estimate[1], 
                      sdlog = fit$estimate[2])
  fitr <- fitdist(xr, 'lnorm')
  pyr <- plnorm(newxs, meanlog = fitr$estimate[1], 
                       sdlog = fitr$estimate[2])
  return(pyr)
}

set.seed(1234) # Установка генератора случайных чисел
require(reshape2)
# новые данные для построения плавной кривой
newxs <- 10^(seq(log10(0.01), log10(max(df$valC)), 
                              length.out = 1000))
# получение матрицы для построения 1000 кривых
boots <- replicate(1000, myboot2(fit, newxs))
bootdat <- data.frame(boots) ; bootdat$newxs <- newxs
bootdat <- melt(bootdat, id = 'newxs')
# извлечение доверительных интервалов
cis <- apply(boots, 1, quantile, c(0.025, 0.975))
rownames(cis) <- c('lwr', 'upr')
# Формирование итоговой таблицы подогнанных значений и  cis
pdat <- data.frame(newxs, py = plnorm(newxs, 
       meanlog = fit$estimate[1], sdlog = fit$estimate[2]))
pdat <- cbind(pdat, t(cis))
# координаты x для названия видов
df$fit <- 10^(log10(qlnorm(df$frac, 
   meanlog = fit$estimate[1], sdlog = fit$estimate[2])) -0.4)

library(ggplot2)
ggplot() +
geom_line(data = bootdat, aes(x = newxs, y = value, 
     group = variable), col = 'steelblue', alpha = 0.05) +
geom_point(data = df, aes(x = valC, y = frac)) +
geom_line(data = pdat, aes(x = newxs, y = py), col = 'red') +
geom_line(data = pdat, aes(x = newxs, y = lwr), 
             linetype = 'dashed') +
geom_line(data = pdat, aes(x = newxs, y = upr), 
             linetype = 'dashed') +
geom_text(data = df, aes(x = fit, y = frac, 
    label = rownames(df)), hjust = 1, size = 4) +
scale_x_log10(breaks = c(5,  10, 50, 100, 300), 
        limits = c(1, max(df$valC))) +
labs(x = 'Содержание алюминия в почве',
            y = 'Доля видов с неблагоприятными условиями') + theme_bw() 



