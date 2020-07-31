# Compare mean ratios for fit_dim = 0 and 2
source('functions.R')

D0 = read.csv(
  paste0(tabRepo,'2020_07_21_09_57_42_fit_dim_0_compilation.csv'),
  check.names = FALSE)
# sel = D0$tag == 'Mean'
sel = !is.na(D0$ratio)
D0 = D0[sel,]
tags = unique(D0$tag)
# icol = sapply(D0$tag, function(x) which(grepl(x,tags)) )
icol = 6
D2 = read.csv(
  paste0(tabRepo,'2020_07_21_11_57_09_fit_dim_2_compilation.csv'),
  check.names = FALSE)
# sel = D2$tag == 'Mean'
D2 = D2[sel,]

pdf(file = paste0(figRepo,'/compare02.pdf'),
    width = 7, height= 10)
par(mfrow=c(4,2),
    mar = mar,
    pty = 's',
    mgp = mgp,
    tcl = tcl)
for(tag in tags) {
  sel = D0$tag == tag
  x  = D2$ratio[sel]
  ux = 2*D2$u_ratio[sel]
  y  = D0$ratio[sel]
  uy = 2*D0$u_ratio[sel]

  xlim = range(c(x-ux,x+ux))
  ylim = range(c(y-uy,y+uy))
  lim  = range(xlim,ylim)
  plot(x,y,
       pch=16,
       xlim = lim, xlab = 'ratio_2',
       ylim = lim, ylab = 'ratio_0',
       col = cols[icol],
       main = tag, cex.main = 0.6)
  grid()
  abline(a=0,b=1)
  segments(x,y-uy,x,y+uy,col = cols[icol])
  segments(x-ux,y,x+ux,y,col = cols[icol])

  rel = (x-y)/x
  plot(x,rel,
       pch=16,
       xlim = lim, xlab = 'ratio_2',
       ylab = 'Rel. error',
       col = cols[icol])
  grid()
  abline(h=c(-0.2,0.2), col = cols[2], lty = 2)
  out = which(abs(rel) > 0.2)
  text(x[out], rel[out], D0$Name[sel][out], adj = 0)

}
dev.off()


# Check u_ratio
# X1 = rnorm(10000, 1, 0.1)
# X2 = rnorm(10000, 2, 0.2)
# R = X1/X2
# r0 = mean(R)
# ur = sd(R)
# ur_lup = r0*sqrt(sd(X1)^2/mean(X1)^2 + sd(X2)^2/mean(X2)^2)
# ur/ur_lup
