library(chemCal)
library(investr)

Data = "arsenic"
data("arsenic")

D = get(Data)
colnames(D)=c('x','y')
m <- lm(y ~ x, data = D)
slope = coefficients(m)[2]
chemCal::calplot(m,xlim = c(0,1),alpha=0.05)
abline(v=0)

# chemCal
print(chemCal::lod(m)$x)

# PP orig
Sxy   = summary(m)$sigma
Lod   = 3 * Sxy / slope
print(Lod)
abline(v=Lod,lty =2, col = 1)

# use predict
Q99 = diff(predict(
  m,
  newdata = data.frame(x=0),
  interval = "prediction",
  level = 0.99
)[c(1,3)])
Lod   = Q99 / slope
print(Lod)
abline(v=Lod, lty =2, col = 2)

# Yet another way
y3 =predict(
  m,
  newdata = data.frame(x=0),
  interval = "prediction",
  level = 0.99
)[3]
Lod = investr::calibrate(m,y0=y3)
print(Lod)
abline(h=y3,lty=2,col=4)
abline(v=unlist(Lod)[1:3],lty =c (2,3,3), col = 4)

xo =D[,1]; yo = D[,2]
plot(
  xo, yo,
  pch = 16,
  col = 4,
  xlab = 'cAA', xlim = c(0,2),
  ylab = 'aire_AA / aire_IS', ylim = c(0,2)
)
grid()
# Fits
x1 = c(0,xo)
p = predict(m,
            newdata = list(x = x1),
            interval = 'pred')
matlines(x1[!is.na(p[,1])],p,
         col = 2,
         lty=c(1,2,2))

# LOD
abline(v=unlist(Lod)[1:3] , col=3, lty = c(1,2,2))
mtext('LOD',side=3,col=3,at=Lod[1], cex=0.75)

box()

stop()

# Regression of the means
weights <- with(D, {
  yx <- split(y, x)
  ybar <- sapply(yx, mean)
  s <- round(sapply(yx, sd), digits = 2)
  w <- round(1 / (s^2), digits = 3)
})
D.means <- aggregate(y ~ x, D, mean)
m.means <- lm(y ~ x, w = weights, data = D.means)
print(lod(m.means,method = "din")$x)

