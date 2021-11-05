fit_1peak = function(
  mu1, sigma1, a1,
  x, y, weights,
  const_fwhm = NA
) {

  s2p = sqrt(2*pi)

  ## 1-peaks fit

  # First pass
  lower = upper = NULL
  start = c(
    mu1    = mu1,
    sigma1 = sigma1,
    a1     = a1
  )
  if(!is.na(const_fwhm)) {
    lower = c(
      mu1    = mu1 - const_fwhm/2,
      sigma1 = 0.8 * const_fwhm/2.355,
      a1     = 0.5 * a1
    )
    start = c(
      mu1    = mu1,
      sigma1 = sigma1,
      a1     = a1
    )
    upper = c(
      mu1    = mu1 + const_fwhm/2,
      sigma1 = 1.2 * const_fwhm/2.355,
      a1     = 2   * a1
    )
  }

  res = try(
    nls(
      y ~ a1/(sqrt(2*pi)*sigma1)*exp(-1/2*(x-mu1)^2/sigma1^2),
      start = start,
      lower = lower,
      upper = upper,
      algorithm = 'port',
      weights = weights,
      control = list(tol=1e-5,warnOnly=TRUE)
    ),
    silent = TRUE
  )

  if(class(res) != 'try-error') {
    if(res$convergence != 0) {
      # Attempt second pass from pass1 optimum
      start = as.list(coef(summary(res))[,"Estimate"])
      res = try(
        nls(
          y ~ a1/(sqrt(2*pi)*sigma1)*exp(-1/2*(x-mu1)^2/sigma1^2),
          start = start,
          lower = lower,
          upper = upper,
          algorithm = 'port',
          weights = weights,
          control = list(tol=1e-5)
        ),
        silent = TRUE
      )
    }
  }

  return(res)

}
fcor_overlap = function(
  mu1, sigma1, a1,
  x, y, weights,
  const_fwhm = NA,
  search_lims = c(2,6),
  ocv_thresh = 0.05
) {

  s2p = sqrt(2*pi)

  # 0. Define search range  out of mu1 +/- 2*sig1
  #    and within mu1 +/- 6*sig1
  sel = abs(x-mu1)/sigma1 >= search_lims[1] &
        abs(x-mu1)/sigma1 <= search_lims[2]
  xs = x[sel]
  ys = y[sel]
  points(xs,ys,col=2)

  # 1. Locate potential peak in search range
  im = which.max(ys)
  mu2    = xs[im]
  sigma2 = sigma1
  ymax   = ys[im]
  a2     = s2p * (ymax-a1*dnorm(mu2,mu1,sigma1)) * sigma2

  # 2. Estimate overlap coefficient value (OCV)
  # OCV = a2 * (pnorm(mu1 + 2 * sigma1, mu2, sigma2) -
  #             pnorm(mu1 - 2 * sigma1, mu2, sigma2))
  # print(OCV)

  # 3. Test OCV and perform 2-peaks fit if necessary
  # if(OCV <= ocv_thresh * a1)
  #   return (NULL)

  ## 2-peaks fit

  # First pass
  lower = upper = NULL
  start = c(
    mu1    = mu1,
    sigma1 = sigma1,
    a1     = a1,
    mu2    = mu2,
    sigma2 = sigma2,
    a2     = a2
  )
  if(!is.na(const_fwhm)) {
    lower = c(
      mu1    = mu1 - const_fwhm/2,
      sigma1 = 0.8 * const_fwhm/2.355,
      a1     = 0.5 * a1,
      mu2    = mu2 - const_fwhm/2,
      sigma2 = 0.8 * const_fwhm/2.355,
      a2     = 0.5 * a2
    )
    start = c(
      mu1    = mu1,
      sigma1 = sigma1,
      a1     = a1,
      mu2    = mu2,
      sigma2 = sigma2,
      a2     = a2
    )
    upper = c(
      mu1    = mu1 + const_fwhm/2,
      sigma1 = 1.2 * const_fwhm/2.355,
      a1     = 2   * a1,
      mu2    = mu2 + const_fwhm/2,
      sigma2 = 1.2 * const_fwhm/2.355,
      a2     = 2   * a2
    )
  }

  res = try(
    nls(
      y ~ a1/(sqrt(2*pi)*sigma1)*exp(-1/2*(x-mu1)^2/sigma1^2) +
          a2/(sqrt(2*pi)*sigma2)*exp(-1/2*(x-mu2)^2/sigma2^2),
      start = start,
      lower = lower,
      upper = upper,
      algorithm = 'port',
      weights = weights,
      control = list(tol=1e-5,warnOnly=TRUE)
    ),
    silent = TRUE
  )

  if(class(res) != 'try-error') {
    if(res$convergence != 0) {
      # Attempt second pass from pass1 optimum
      start = as.list(coef(summary(res))[,"Estimate"])
      res = try(
        nls(
          y ~ a1/(sqrt(2*pi)*sigma1)*exp(-1/2*(x-mu1)^2/sigma1^2) +
              a2/(sqrt(2*pi)*sigma2)*exp(-1/2*(x-mu2)^2/sigma2^2),
          start = start,
          lower = lower,
          upper = upper,
          algorithm = 'port',
          weights = weights,
          control = list(tol=1e-5)
        ),
        silent = TRUE
      )
    }
  }

  return(res)

}

# x = seq(0,10, by=0.1)
#
# const_fwhm = 0.7
#
# mu1    = 3
# sigma1 = 0.7
# a1     = 1
# mu2    = 6
# sigma2 = 0.7
# a2     = 1
#
# y = a1 * dnorm(x,mu1,sigma1) +
#   a2 * dnorm(x,mu2,sigma2) +
#   2*a2 * dnorm(x,mu2+3,sigma2) +
#   rnorm(length(x),0,0.05)
# weights = rep(1,length(y))
# plot(x,y,type ='l')
#
# res0 = fit_1peak(
#   mu1, sigma1, a1,
#   x, y, weights,
#   const_fwhm = const_fwhm
# )
# print(res0)
# v = coefficients(res0)
# lines(x,v['a1']*dnorm(x,v['mu1'],v['sigma1']), col= 4, lty=2)

const_fwhm = 1.3
sel = CV > -5 & CV <5
x = CV[sel]
y = mMStot[sel]
weights = rep(1,length(y))
plot(x,y, type ='l')

v = coefficients(fitOut$res)
mu1 = 0
sigma1 = 0.5 #v['sigma']
a1 = 100

res0 = fit_1peak(
  mu1, sigma1, a1,
  x, y, weights,
  const_fwhm = const_fwhm
)
# print(res0)
v = coefficients(res0)
print(v)
lines(x,v['a1']*dnorm(x,v['mu1'],v['sigma1']), col= 4, lty=2)
cat('1-peak: fwhm=',2.355*v['sigma1'],'\n\n')

if (2.355*v['sigma1'] >= 1.2*const_fwhm ) {
  res = fcor_overlap(
    mu1, sigma1, a1,
    x, y, weights,
    const_fwhm = const_fwhm,
    search_lims = c(2,8),
    ocv_thresh = 0.01
  )
  # print(res)
  if(!is.null(res)) {
    v = coefficients(res)
    print(v)
    lines(x,v['a1']*dnorm(x,v['mu1'],v['sigma1']), col= 4, lty=1, lwd = 2)
    lines(x,v['a2']*dnorm(x,v['mu2'],v['sigma2']), col= 3, lty=1, lwd = 2)
    lines(x,v['a1']*dnorm(x,v['mu1'],v['sigma1']) +
            v['a2']*dnorm(x,v['mu2'],v['sigma2']), col= 2, lty=1, lwd = 2)
  }
}
cat('2-peak: fwhm=',2.355*v['sigma1'],'\n')
