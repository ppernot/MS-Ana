library(chemCal)
library(investr)

# Data = "arsenic"
# data("arsenic")
# D = get(Data)
# colnames(D)=c('x','y')

file = 'data_6J_MMA.csv'
# file = 'data_6J_Glut.csv'
# file = 'data_6J_Adip.csv'

D = read.csv(file,sep =';',header=FALSE)
x = D[,1]
Y = D[,2:7]
Ns = ncol(Y)

regl = list()
tabc = tabu = matrix(NA,ncol=2,nrow=Ns)
for(i in 1:Ns){
  y = Y[,i]
  regl[[i]] = lm(y ~ x)
  tabc[i,]  = summary(regl[[i]])$coefficients[,1]
  tabu[i,]  = summary(regl[[i]])$coefficients[,2]
}
c.mean = apply(tabc,2,mean)
c.sd   = apply(tabc,2,sd)

# Add variance of fit
u.mean = apply(tabu,2,function(x) sqrt(mean(x^2)))
u.tot  = sqrt(c.sd^2 + u.mean^2)

# # Validate u.tot by simulation
# vv = c(); k=0
# for(i in 1:1000) {
#   for(j in 1:Ns) {
#     k = k+1
#     vv[k] = rnorm(1,tabc[j,1],tabu[j,1])
#   }
# }
# sd(vv)

LOD_orig = 3 * c.sd[1] / c.mean[2]

LOD_new  = qt(0.99,df = Ns-1) * u.tot[1] / c.mean[2]

matplot(
  x, Y,
  type='p', pch = 16, col = 4,
  xlim = c(0,10),
  ylim = c(0,0.2),
  # xaxs = 'i',
  yaxs = 'i',
  main = file
)
grid()
for(i in 1:Ns)
  abline(regl[[i]], lty=1, col='gray50')
abline(v=0,h=0)
abline(v = LOD_orig, col=1, lwd =2, lty = 3)
abline(v = LOD_new,  col=2, lwd =2, lty = 3)


xo = rep(x,Ns)
yo = unlist(Y)
m = lm(yo ~ xo)
# chemCal::calplot(
#   m,
#   xlim = c(0,10),
#   ylim= c(0,0.2),
#   alpha=0.05)
# abline(v=0,h=0)

p = predict(
  m,
  newdata = data.frame(xo = unique(xo)),
  interval = 'pred',
  level = 0.99)
matlines(unique(xo), p,
         col = 4, lwd=2,
         lty=c(1,2,2))

# chemCal
Lod1 = chemCal::lod(m)$x
abline(v=Lod1,col = 3, lwd =2, lty = 3)

# PP orig
Sxy   = summary(m)$sigma
slope = coefficients(m)[2]
Lod   = 3 * Sxy / slope
abline(v=Lod,col = 4, lwd =2, lty = 3)

legend(
  'bottomright', bty='n',
  title = 'LOD',
  legend = c(
    paste0(signif(LOD_orig,4),' (Excel)'),
    paste0(signif(LOD_new,4), ' (Excel_Stud)'),
    paste0(signif(Lod1,4),    ' (Chemcal)'),
    paste0(signif(Lod,4),     ' (RepAndQuant)')
  ),
  col = 1:4,
  lty = 3, lwd = 2,
  pch = NA
)


stop()

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
chemCal::calplot(
  m.means,
  xlim = c(0,60),
  ylim= c(0,5),
  alpha=0.05)
abline(v=0,h=0)

# print(lod(m.means,method = "din")$x)

