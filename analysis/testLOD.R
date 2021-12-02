library(chemCal)
library(investr)

# Data = "arsenic"
# data("arsenic")
# D = get(Data)
# colnames(D)=c('x','y')

file = 'data_6J_MMA.csv'
# file = 'data_6J_Glut.csv'
# file = 'data_6J_Adip.csv'

D = read.csv(file, sep = ';', header = FALSE)
x = D[, 1]
Y = D[, 2:ncol(D)]
Ns = ncol(Y)

regl = list()
tabc = matrix(NA, ncol = 2, nrow = Ns)
rownames(tabc) = paste0('J',1:Ns)
colnames(tabc) = c('Intercept', 'slope')
tabu = tabc
colnames(tabu) = c('u_Intercept', 'u_slope')
tabl = tablu = c()
for (i in 1:Ns) {
  y = Y[, i]
  regl[[i]]  = lm(y ~ x)
  tabc[i, ]  = summary(regl[[i]])$coefficients[, 1]
  tabu[i, ]  = summary(regl[[i]])$coefficients[, 2]
  Sxy        = summary(regl[[i]])$sigma
  slope      = tabc[i,2]
  uSlope     = tabu[i,2]
  tabl[i]    = 2 * qnorm(0.95) * Sxy / slope
  tablu[i]   = 2 * qnorm(0.95) * Sxy / slope^2 * uSlope

}
c.mean = apply(tabc, 2, mean)
c.sd   = apply(tabc, 2, sd)
LOD_mean = mean(tabl)
u.LOD_mean   = sd(tabl) / sqrt(length(tabl))

tabw    = 1 / tablu^2
tabw    = tabw / sum(tabw)
LOD_wmean = sum(tabl*tabw)
u.LOD_wmean = sqrt(sum(tabw*(tabl^2+tablu^2))-LOD_wmean^2)

# Weighted alternatives
# https://en.wikipedia.org/wiki/Mixture_distribution
tabw    = 1 / tabu^2
tabw    = apply(tabw,2,function(x) x/sum(x))
c.wmean = colSums(tabc*tabw)
c.wsd   = sqrt(colSums(tabw*(tabc^2+tabu^2))-c.wmean^2)

u.tot = sqrt(c.sd^2 + colMeans(tabu^2))

# MC check of u.tot
V = c()
for(i in 1:1e4) {
  z = c()
  for (j in 1:Ns) {
    z[j] = rnorm(1,tabc[j,1],tabu[j,1])
  }
  V[i] = var(z)
}
u.totMC = sqrt(mean(V))
print(u.totMC)

xs = 1:Ns
ys = tabc[,1]
uys = tabu[,1]
plot(xs,ys, xlab = 'Series', ylab = 'Intercept',
     ylim = c(0,0.07))
segments(xs,ys-2*uys,xs,ys+2*uys)
abline(h=c(c.mean[1],c.mean[1]-2*c.sd[1],c.mean[1]+2*c.sd[1]),
       col=2,lty=c(1,2,2))
abline(h=c(c.mean[1],c.mean[1]-2*u.totMC,c.mean[1]+2*u.totMC),
       col=3,lty=c(1,2,2))
abline(h=c(c.mean[1],c.mean[1]-2*u.tot[1],c.mean[1]+2*u.tot[1]),
       col=4,lty=c(1,2,2))


cat('\n --- Line params by series + Mean + SD ---')
print(
  knitr::kable(
    rbind(
      cbind(tabc,tabu,LOD=tabl,u_LOD=tablu),
      Mean=c(c.mean,NA,NA,LOD_mean,LOD_wmean),
      SD = c(c.sd, NA,NA,u.LOD_mean,u.LOD_wmean),
      u_tot = c(NA,NA,u.tot,NA,NA)
    ),
    digits = 5
  )
)

# LOD estimators

LOD_excel = 3 * c.sd[1] / c.mean[2]

# LOD_wgt   = 3 * c.wsd[1] / c.wmean[2]

LOD_stud  = 2 * qnorm(0.95) * u.tot[1] / c.wmean[2]

# Aggregated data
xo = rep(x, Ns)
yo = unlist(Y)
m = lm(yo ~ xo)

# IUPAC
LOD_chemcal = chemCal::lod(m)$x

# RepAndQuant
Sxy     = summary(m)$sigma
slope   = coefficients(m)[2]
LOD_RAQ = 3 * Sxy / slope

# # Regression of the means
# yx <- split(yo, xo)
# ybar <- sapply(yx, mean)
# s <- sapply(yx, sd)
# w <- 1 / (s^2)
# Y.means <- aggregate(y ~ x, cbind(xo,yo), mean)
# m.means <- lm(y ~ x, data = Y.means, weights = w)
#
# x.mu  = coef(m.means)
# x.cov = vcov(m.means)
# x.cov = x.cov /summary(m.means)$sigma^2
# x.u   = diag(x.cov)^0.5 # Standard error
#
# LOD_RAQ_m = 2 * qt(0.95, df = Ns - 1) * x.u[1] / x.mu[2]

matplot(
  x,
  Y,
  type = 'p',
  pch = 16,
  col = 4,
  xlim = c(0, 10),
  ylim = c(0, 0.75*mean(unlist(Y))),
  # xaxs = 'i',
  yaxs = 'i',
  main = file
)
grid()
for (i in 1:Ns)
  abline(regl[[i]], lty = 1, col = 'gray50')
abline(v = 0, h = 0)
abline(v = LOD_excel,
       col = 1,
       lwd = 2,
       lty = 3)
abline(v = tabl,
       col = 2,
       lwd = 2,
       lty = 3)
abline(v = LOD_mean,
       col = 2,
       lwd = 4,
       lty = 2)
segments(tabl-2*tablu,0.15,
         tabl+2*tablu,0.15,
         col = 2,lwd=4)
segments(LOD_mean-2*u.LOD_mean,0,
         LOD_mean+2*u.LOD_mean,0,
         col = 2,lwd=8)

abline(v = LOD_wmean,
       col = 6,
       lwd = 4,
       lty = 2)
segments(LOD_wmean-2*u.LOD_wmean,0,
         LOD_wmean+2*u.LOD_wmean,0,
         col = 6,lwd=8)




abline(v = LOD_stud,
       col = 3,
       lwd = 2,
       lty = 3)

p = predict(
  m,
  newdata = data.frame(xo = unique(xo)),
  interval = 'pred',,
  # se.fit = TRUE,
  level = 0.99
)
matlines(unique(xo),
         p[,1],
         col = 4,
         lwd = 2,
         lty = c(1, 2, 2))

abline(v = LOD_chemcal,
       col = 4,
       lwd = 2,
       lty = 3)

abline(v = LOD_RAQ,
       col = 5,
       lwd = 2,
       lty = 3)

legend(
  'bottomright',
  bty = 'n',
  title = 'LOD',
  legend = c(
    paste0(signif(LOD_excel, 4),   ' (Excel)'),
    paste0(signif(LOD_mean, 4),    ' (Mean)'),
    paste0(signif(LOD_stud, 4),    ' (Studentized)'),
    paste0(signif(LOD_chemcal, 4), ' (Chemcal)'),
    paste0(signif(LOD_RAQ, 4),     ' (RepAndQuant)')
  ),
  col = 1:5,
  lty = 3,
  lwd = 2,
  pch = NA
)

stop()

# Regression of the means
weights <- with(D, {
  yx <- split(y, x)
  ybar <- sapply(yx, mean)
  s <- sapply(yx, sd)
  w <- 1 / (s^2)
})
D.means <- aggregate(y ~ x, D, mean)
m.means <- lm(y ~ x, w = weights, data = D.means)

# chemCal::calplot(
#   m.means,
#   xlim = c(0,10),
#   ylim= c(0,0.25),
#   alpha=0.05)
# abline(v=0,h=0)
#
# print(lod(m.means,method = "din")$x)
#
