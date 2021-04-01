# Code Setup ####
options( warn=0,
         stringsAsFactors = FALSE )

# Install packages if necessary
## CRAN packages
libs <- c('xtable','mixtools','inlmisc',
          'rlist','repmis','assertive',
          'data.table')
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE #,
      # repos = 'https://cran.irsn.fr'
    )
  }
}

## Load packages
repmis::LoadandCite(libs) # ,file='../article/packages.bib')

# Set graphical params
## For PNG figures
gPars = list(
  cols     = rev(inlmisc::GetColors(8))[1:7],
  cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.2))[1:7],
  cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.5))[1:7],
  pty      = 's',
  mar      = c(3,3,3,.5),
  mgp      = c(2,.75,0),
  tcl      = -0.5,
  lwd      = 4.0,
  cex      = 4.0,
  cex.leg  = 0.7,
  reso     = 1200  # (px) base resolution for png figs
)
## For local plots
gParsLoc = gPars
gParsLoc$cex = 1
gParsLoc$lwd = 1.5

# Expose gPars list
for (n in names(gPars))
  assign(n, rlist::list.extract(gPars, n))

# Define Data and Results repositories
dataRepo = '../data/'
figRepo  = '../results/figs/'
tabRepo  = '../results/tables/'

# Default user tag
userTag  = ''

sink(file ='./sessionInfo.txt')
print(sessionInfo(), locale=FALSE)
sink()

# Functions ####

getMS = function(file, ms_type = 'esquire') {

  cat('\n>>> Reading',ms_type,'MS in file:',file)
  MS0 = as.data.frame(
    data.table::fread( # Much faster than read.table
      file = file ,
      header = FALSE ,
      sep = ',' ,
      stringsAsFactors = FALSE
    )
  )
  time = MS0[, 1]

  if(ms_type == 'esquire') {
    range_mz = as.numeric(unlist(strsplit(MS0[1, 7], split = '-')))
    n_del_mz = as.numeric(unlist(strsplit(MS0[1, 8], split = '/')))
    nchan    = n_del_mz[1]
    del_mz   = n_del_mz[2]
    mz       = range_mz[1] + (0:(nchan - 1)) * del_mz
    MS       = as.matrix(MS0[, -(1:8)],
                         ncol = length(mz),
                         byrow = FALSE)
  } else {
    mySplit = function(x) as.numeric(unlist(strsplit(x, split = ' ')))
    mz = c()
    msLen = ncol(MS0)-8
    MS = matrix(NA, nrow = nrow(MS0), ncol = msLen)
    for (j in 1:nrow(MS0)) {
      dbl = vapply(
        MS0[j, 9:(msLen + 8)],
        FUN = mySplit,
        FUN.VALUE = numeric(2)
      )
      mz      = dbl[1, ]
      MS[j, ] = dbl[2, ]
    }
  }

  return(
    list(
      mz = mz,
      time = time,
      MS = MS
    )
  )
}

makeTag <- function(CVTable, msTable, userTag) {
  date =
    strsplit(
      strsplit(
        CVTable,
        split = ' '
      )[[1]][2],
      split = '-'
    )[[1]][1]

  tag = paste0(
    date,'_',
    strsplit(msTable, split='\\.')[[1]][1],
    ifelse(userTag=='','','_'),userTag
  )
}
readTasksFile <- function(file) {
 read.table(
    file = file,
    header = TRUE,
    sep = ',',
    check.names = FALSE
  )
}
gatherResults <- function(Tasks, tabRepo, userTag) {
  D = NULL
  for(task in 1:nrow(Tasks)) {

    msTable = Tasks[task,'MS_file']
    CVTable = Tasks[task,'DMS_file']
    # dilu    = Tasks[task,'dilu']

    # Build tag
    tag = makeTag(CVTable, msTable, userTag)
    file = paste0(tabRepo,tag,'_results.csv')
    if(!file.exists(file))
      stop(paste0('Missing file:',file))

    M = read.csv(file = file, check.names = FALSE)
    # M = cbind(M,dilu)
    D = rbind(D,M)
  }
  return(D)
}
readTargetsFile <- function(file) {
  M = read.table(
    file = file,
    header = TRUE,
    sep = ',',
    dec = '.',
    check.names = FALSE,
    fill = TRUE
  )
  # Backwards separator compatibility
  if(ncol(M) < 4)
    M = read.table(
      file = file,
      header = TRUE,
      sep = ';',
      dec = '.',
      check.names = FALSE,
      fill = TRUE
    )
  # Backwards colnames compatibility
  if('m/z_exact' %in% colnames(M))
    colnames(M)[colnames(M)=='m/z_exact'] = 'm/z_ref'

  return(M)
}
peak_shape = function(x, p) {
  if(length(p)==3)
    p['A'] / (sqrt(2*pi) * p['sigma'] ) *
    exp(-1/2*(x-p['mu'])^2/p['sigma']^2)
  else
    p['A'] / (sqrt(2*pi) * p['sy'] ) * # Marginalized on x
    exp(-1/2*(x-p['my'])^2/p['sy']^2)
}
plotPeak = function(
  mz, CV, MS,
  fitOut,
  mex = NA,
  leg = NA,
  val = NA,
  tag = NA,
  type = 'CV',
  CV0 = NA,
  gPars
) {

  nCV = length(CV)

  # Expose fitOut list
  for (n in names(fitOut))
    assign(n,rlist::list.extract(fitOut,n))

  mzlim  = c(mz1,mz2)
  CVlim  = range(CV)
  if(type == 'CV')
    CVlimf = range(CVf)
  else
    CVlimf = range(CV)

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 2),
    mar   = mar,
    mgp   = mgp,
    tcl   = tcl,
    lwd   = lwd,
    cex   = cex
  )

  ## 1. image
  sel1 = apply(cbind(CV-CVlim[1],CV-CVlim[2]),1,prod) <= 0
  sel2 = apply(cbind(mz-mzlim[1],mz-mzlim[2]),1,prod) <= 0

  image(
    CV[sel1], mz[sel2], MS[sel1,sel2],
    xlim = CVlim,
    xlab = 'CV',
    ylim = mzlim,
    ylab = 'm/z'
  )
  grid()

  if(class(res) != 'try-error') {
    v   = summary(res)$parameters[,"Estimate"]
    if(length(v) > 3) {
      # Fit 2D
      p = matrix(
        predict(res),
        ncol = ncol(MSloc),
        nrow = nrow(MSloc),
        byrow = TRUE
      )
      contour(CVf, mzloc, p, col = cols_tr2[5], add=TRUE)
    }
    if(type == 'CV')
      rect(CVlimf[1],mz1,CVlimf[2],mz2,
           col = cols_tr[4], border=NA)
  }


  if (!is.na(leg))
    title(main = leg, line = 0.5)

  if (!is.na(tag))
    mtext(
      text = tag,
      side = 3,
      line = 2,
      cex = cex
    )

  if (!is.na(CV0))
    abline(v = CV0, lty = 1, col = cols[2])

  if (!is.na(mex))
    abline(h = mex, lty = 1, col = cols[2])

  if (type != 'CV') {
    abline(h = v[1], lty = 2, col = cols[2])
  } else {
    if( length(v) == 3) {
      abline(v = v[1], lty = 2, col = cols[2])
      if(!is.na(mz0))
        abline(h = mz0, lty = 2, col = cols[2])
    } else {
      abline(h = v[1], lty = 2, col = cols[2])
      abline(v = v[3], lty = 2, col = cols[2])
    }
  }

  if(type == 'CV') {
    ## 2. CV profile
    xmod = seq(min(CV),max(CV),length.out = 1000)
    if(class(res) != 'try-error') {
      v   = summary(res)$parameters[,"Estimate"]
      vmod = peak_shape(xmod, v)
    } else {
      v = NA
      vmod = rep(NA,length(xmod))
    }
    ylim = range(c(mMS,vmod), na.rm = TRUE, finite = TRUE)
    plot(
      CV, mMStot,
      type = 'l', #pch = 16,
      col  = cols[4],
      xlim = CVlim,
      xlab = 'CV',
      ylim = ylim,
      ylab = 'a.u.'
    )
    title(main = 'Integ. CV profile', line = 0.5)
    rect(CVlimf[1], -100, CVlimf[2], ylim[2]*1.2,
         col = cols_tr[4], border=NA)

    ## 2.1 Gaussian fit
    if(!any(is.na(v))) {
      lines(xmod, vmod, col = cols[6])
      if(length(v) == 3)
        abline(v = v[1], lty = 2, col = cols[2])
      else
        abline(v = v[3], lty = 2, col = cols[2])
    }

    if (!is.na(CV0))
      abline(v = CV0, lty = 1, col = cols[2])

    ## Add fit results
    if (!is.na(val))
      legend(
        x=CVlim[1], y=ylim[2],
        yjust = 1.75,
        title  = val,
        legend = '',
        # inset  = 0.12,
        bty    = 'n',
        cex    = 0.8)

  } else {

    ## 2. MS profile
    xmod = seq(mz1,mz2,length.out = 1000)
    if(class(res) != 'try-error') {
      v   = summary(res)$parameters[,"Estimate"]
      vmod = peak_shape(xmod, v)
    } else {
      v = NA
      vmod = rep(NA,length(xmod))
    }
    ylim = c(
      0,
      1.2*max(c(MSl,vmod), na.rm = TRUE, finite = TRUE)
    )
    plot(
      mzl, MSl,
      type = 'l', #pch = 16,
      col  = cols[4],
      xlim = mzlim,
      xlab = 'm/z',
      ylim = ylim,
      ylab = 'a.u.',
      yaxs = 'i'
    )
    title(main = 'MS profile', line = 0.5)

    ## 2.1 Gaussian fit
    if(!any(is.na(v))) {
      lines(xmod, vmod, col = cols[6])
      abline(v = v[1], col = cols[2], lty = 2)
    }

    if (!is.na(mex))
      abline(v = mex, lty = 1, col = cols[2])

    ## Add fit results
    if (!is.na(val))
      legend('right',
        # x=mzlim[1], y=ylim[2],
        yjust = 1.75,
        title  = val,
        legend = '',
        # inset  = 0.05,
        bty    = 'n',
        cex    = 0.8)

  }

}
plotMaps = function(
  mz, CV, MS,
  mex = NA,
  leg = NA,
  tag = NA,
  mzlim = range(mz),
  CVlim = range(CV),
  zlim = c(0,max(MS)/2),
  logz = TRUE,
  gPars
) {

  nCV = length(CV)

  # Expose gPars list
  for (n in names(gPars))
    assign(n,rlist::list.extract(gPars,n))

  par(
    mfrow = c(1, 2),
    mar   = mar,
    mgp   = mgp,
    tcl   = tcl,
    lwd   = lwd,
    cex   = cex
  )

  M = MS
  if(logz) {
    M = log10(MS)
    zlim =range(M, finite = TRUE)
  }

  ## 1. image

  image(
    CV, mz, M,
    xlab = 'CV',
    ylab = 'm/z',
    zlim = zlim
  )
  grid()

  if(!is.na(leg))
    title(main = leg, line = 0.5)

  if(!is.na(tag))
    mtext(text = tag, side = 3, line = 2, cex = cex)

  if(!any(is.na(mex)))
    abline(h=mex,lty=1,col=cols_tr2[5])

  ## 1. image
  sel1 = apply(cbind(CV-CVlim[1],CV-CVlim[2]),1,prod) <= 0
  sel2 = apply(cbind(mz-mzlim[1],mz-mzlim[2]),1,prod) <= 0

  image(
    CV[sel1], mz[sel2], M[sel1,sel2],
    xlim = CVlim,
    xlab = 'CV',
    ylim = mzlim,
    ylab = 'm/z',
    zlim = zlim
  )
  grid()

  if(!is.na(leg))
    title(main = leg, line = 0.5)

  if(!any(is.na(mex)))
    abline(h=mex,lty=1,col=cols_tr2[5])

}
fit1D_MS <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  weighted = FALSE,
  const_fwhm = NA
) {

  # Select CV0
  iCV = which(CV >= CV0)[1]

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # mz profile
  MStot = MS[iCV,]
  mzl = mz[selMz]
  MSl = MS[iCV, selMz]

  # Normal fit
  weights = rep(1,length(MSl))
  if(weighted){ # Poisson
    weights = 1 / MSl
    vm = min(MSl[MSl>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  mu = mzl[which.max(MSl)]
  lower = NULL
  s2p = sqrt(2*pi)
  sigma0 = 0.1/2.355
  A0 = s2p * sigma0 * max(MSl)
  start = c(
    mu    = mu,
    sigma = sigma0,
    A     = A0
  )
  upper = NULL
  if(!is.na(const_fwhm)) {
    sigma0 = const_fwhm/2.355
    A0 = s2p * sigma0 * max(MSl)
    lower = c(
      mu    = mu - dmz,
      sigma = 0.8 * sigma0,
      A     = 0.5 * A0
    )
    start = c(
      mu    = mu,
      sigma = sigma0,
      A     = A0
    )
    upper = c(
      mu    = mu + dmz,
      sigma = 1.2 * sigma0,
      A     = 1.5 * A0
    )
  }

  res = try(
    nls(
      MSl ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(mzl-mu)^2/sigma^2),
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
          MSl ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(mzl-mu)^2/sigma^2),
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
    # if(class(res) != 'try-error')
    #   print(cbind(lower,coef(summary(res))[,"Estimate"],upper))
  }

  return(
    list(
      mzl = mzl,
      MSl = MSl,
      res = res,
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      CVf = CV,
      mMStot = MStot
    )
  )
}
trapz = function (x, y) {
  if(length(x) == 0)
    return(0)
  idx = 2:length(x)
  return(
    as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
  )
}
fit1D_old <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  del_mz,
  weighted = FALSE,
  const_fwhm = NA,
  refine_CV0 = FALSE,
  correct_overlap = FALSE
) {
print("enter...")
  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # Integrated CV profile
  mMStot = rowSums(MS[, selMz]) * del_mz

  # Select CV window
  if(is.na(CV0))
    CV0 = CV[which.max(mMStot)]
  CV1 = CV0 - dCV/2
  CV2 = CV0 + dCV/2
  selCV = CV >= CV1 & CV <= CV2

  if(refine_CV0){
    im = which.max(mMStot[selCV])
    CV0 = CV[selCV][im]
    CV1 = CV0 - dCV/2
    CV2 = CV0 + dCV/2
    selCV = CV >= CV1 & CV <= CV2
  }
  CVf = CV[selCV]

  # Refine mz window
  MSloc  = MS[selCV, selMz]
  mMS = rowSums(MSloc)*del_mz # Sum over selected mz
  mz_max = which.max(MSloc[which.max(mMS),])
  mz0 = mz[selMz][mz_max]
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz

  selMz  = mz >= mz1 & mz <= mz2 # Select mz area

  # Normal fit
  mMS = rowSums(MS[selCV, selMz])*del_mz # Sum over selected mz
  weights = rep(1,length(mMS))
  if(weighted){ # Poisson
    weights = 1 / mMS
    vm = min(mMS[mMS>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  lower = NULL
  s2p = sqrt(2*pi)
  sigma0 = 0.7/2.355
  A0 = s2p * sigma0 * max(mMS)
  start = c(
    mu    = CVf[which.max(mMS)],
    sigma = sigma0,
    A     = A0
  )
  upper = NULL
  if(!is.na(const_fwhm)) {
    sigma0 = const_fwhm/2.355
    A0 = s2p * sigma0 * max(mMS)
    lower = c(
      mu    = CV0 - dCV/10,
      sigma = 0.8 * sigma0,
      A     = 0.5 * A0
    )
    start = c(
      mu    = CV0,
      sigma = sigma0,
      A     = A0
    )
    upper = c(
      mu    = CV0 + dCV/10,
      sigma = 1.2 * sigma0,
      A     = 1.5 * A0
    )
  }

  res = try(
    nls(
      mMS ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(CVf-mu)^2/sigma^2),
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
          mMS ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(CVf-mu)^2/sigma^2),
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
    # if(class(res) != 'try-error')
    #   print(cbind(lower,coef(summary(res))[,"Estimate"],upper))
  }

  return(
    list(
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      res = res,
      mMS = mMS,
      mMStot = mMStot,
      CVf = CVf
    )
  )
}
fit1D <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  weighted = FALSE,
  const_fwhm = NA,
  refine_CV0 = FALSE,
  correct_overlap = FALSE
) {

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # Integrated CV profile
  mzloc = mz[selMz]
  mMStot = apply(MS[, selMz], 1, function(x) trapz(mzloc,x))

  # Select CV window
  if(is.na(CV0))
    CV0 = CV[which.max(mMStot)]
  CV1 = CV0 - dCV/2
  CV2 = CV0 + dCV/2
  selCV = CV >= CV1 & CV <= CV2

  if(refine_CV0){
    im = which.max(mMStot[selCV])
    CV0 = CV[selCV][im]
    CV1 = CV0 - dCV/2
    CV2 = CV0 + dCV/2
    selCV = CV >= CV1 & CV <= CV2
  }
  CVf = CV[selCV]

  # Refine mz window
  MSloc  = MS[selCV, selMz]
  mzloc  = mz[selMz]
  mMS    = apply(MSloc, 1, function(x) trapz(mzloc, x))

  mz_max = which.max(MSloc[which.max(mMS),])
  mz0 = mz[selMz][mz_max]
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area

  # Normal fit
  mzloc = mz[selMz]
  mMS = apply(MS[selCV, selMz], 1, function(x) trapz(mzloc,x))

  # Weighting scheme
  weights = rep(1,length(mMS))
  if(weighted){ # Poisson
    weights = 1 / mMS
    vm = min(mMS[mMS>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  lower = NULL
  s2p = sqrt(2*pi)
  sigma0 = 0.7/2.355
  A0 = s2p * sigma0 * max(mMS)
  start = c(
    mu    = CVf[which.max(mMS)],
    sigma = sigma0,
    A     = A0
  )
  upper = NULL
  if(!is.na(const_fwhm)) {
    sigma0 = const_fwhm/2.355
    A0 = s2p * sigma0 * max(mMS)
    lower = c(
      mu    = CV0 - dCV/10,
      sigma = 0.8 * sigma0,
      A     = 0.5 * A0
    )
    start = c(
      mu    = CV0,
      sigma = sigma0,
      A     = A0
    )
    upper = c(
      mu    = CV0 + dCV/10,
      sigma = 1.2 * sigma0,
      A     = 1.5 * A0
    )
  }

  res = try(
    nls(
      mMS ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(CVf-mu)^2/sigma^2),
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
          mMS ~ A/(sqrt(2*pi)*sigma)*exp(-1/2*(CVf-mu)^2/sigma^2),
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
    # if(class(res) != 'try-error')
    #   print(cbind(lower,coef(summary(res))[,"Estimate"],upper))
  }

  return(
    list(
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      res = res,
      mMS = mMS,
      mMStot = mMStot,
      CVf = CVf
    )
  )
}
fit2D <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  weighted = FALSE,
  const_fwhm = NA,
  refine_CV0 = FALSE,
  correct_overlap = FALSE
) {

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # Integrated CV profile
  mzloc = mz[selMz]
  mMStot = apply(MS[, selMz], 1, function(x) trapz(mzloc,x))

  # Select CV window
  if(is.na(CV0))
    CV0 = CV[which.max(mMStot)]
  CV1 = CV0 - dCV/2
  CV2 = CV0 + dCV/2
  selCV = CV >= CV1 & CV <= CV2

  if(refine_CV0){
    im = which.max(mMStot[selCV])
    CV0 = CV[selCV][im]
    CV1 = CV0 - dCV/2
    CV2 = CV0 + dCV/2
    selCV = CV >= CV1 & CV <= CV2
  }
  CVf = CV[selCV]

  # Refine mz window
  MSloc = MS[selCV, selMz]
  mzloc = mz[selMz]
  mMS = apply(MSloc, 1, function(x) trapz(mzloc,x))
  mz_max = which.max(MSloc[which.max(mMS),])
  mz0 = mz[selMz][mz_max]
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area

  MSloc = MS[selCV, selMz]
  mzloc = mz[selMz]
  mMS = apply(MSloc, 1, function(x) trapz(mzloc,x))

  # Normal fit
  grid = expand.grid(x = mzloc, y = CVf)
  x = grid$x
  y = grid$y
  z = as.vector(t(MSloc))

  weights = rep(1,length(z))
  if(weighted){ # Poisson
    weights = 1 / z
    vm = min(z[z>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  maxz = which.max(z)

  if(!is.na(const_fwhm)) {
    sx0 = 0.2
    sy0 = const_fwhm / 2.355
    A0  = 2 * pi * sx0 * sy0 * max(z)
    lower = c(
      mx  = x[maxz] - dmz / 2,
      sx  = sx0 / 4,
      my  = CV0 - dCV / 10,
      sy  = 0.8 * sy0,
      A   = 0.5 * A0
    )
    start = c(
      mx  = x[maxz],
      sx  = sx0,
      my  = CV0,
      sy  = sy0,
      A   = A0
    )
    upper = c(
      mx  = x[maxz] + dmz / 2,
      sx  = 4 * sx0,
      my  = CV0 + dCV / 10,
      sy  = 1.2 * sy0,
      A   = 1.5 * A0
    )

  } else {
    sx0 = 0.2
    sy0 = 0.7/2.355
    A0  = 2 * pi * sx0 * sy0 * max(z)

    lower = NULL
    start = c(
      mx  = x[maxz],
      sx  = sx0,
      my  = y[maxz],
      sy  = sy0,
      A   = A0
    )
    upper = NULL
  }

  res = try(
    nls(
      z ~ A/(2*pi*sx*sy)*exp(
        -0.5 * ( (x-mx)^2/sx^2 + (y-my)^2/sy^2 ) ),
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
          z ~ A/(2*pi*sx*sy)*exp(
            -0.5 * ( (x-mx)^2/sx^2 + (y-my)^2/sy^2 ) ),
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

  return(
    list(
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      res = res,
      mMS = mMS,
      mMStot = mMStot,
      CVf = CVf,
      mzloc = mzloc,
      MSloc = MSloc
    )
  )
}
getPars1D <- function(res, fit_dim) {
  # Get best params and uncertainty
  v   = summary(res)$parameters[,"Estimate"]   # Best params
  u_v = summary(res)$parameters[,"Std. Error"] # Uncertainty

  # Transform params to quantities of interest
  mu     = v['mu']
  u_mu   = u_v['mu']
  fwhm   = 2.355 * v['sigma']
  u_fwhm = 2.355 * u_v['sigma']
  area   = v['A']
  u_area = u_v['A']

  if (fit_dim == 1)
    return(
      list(
        v         = v,
        u_v       = u_v,
        mzopt     = NA,
        u_mz      = NA,
        cvopt     = mu,
        u_cv      = u_mu,
        fwhm_mz   = NA,
        u_fwhm_mz = NA,
        fwhm_cv   = fwhm,
        u_fwhm_cv = u_fwhm,
        area      = area,
        u_area    = u_area
      )
    )
  else
    return(
      list(
        v         = v,
        u_v       = u_v,
        mzopt     = mu,
        u_mz      = u_mu,
        cvopt     = NA,
        u_cv      = NA,
        fwhm_mz   = fwhm,
        u_fwhm_mz = u_fwhm,
        fwhm_cv   = NA,
        u_fwhm_cv = NA,
        area      = area,
        u_area    = u_area
      )
    )

}
getPars2D <- function(res) {
  # Get best params and uncertainty
  v   = summary(res)$parameters[,"Estimate"]   # Best params
  u_v = summary(res)$parameters[,"Std. Error"] # Uncertainty

  # Transform params to quantities of interest
  mu1     = v['mx']
  u_mu1   = u_v['mx']
  mu2     = v['my']
  u_mu2   = u_v['my']
  fwhm1   = 2.355 * v['sx']
  u_fwhm1 = 2.355 * u_v['sx']
  fwhm2   = 2.355 * v['sy']
  u_fwhm2 = 2.355 * u_v['sy']
  area    = v['A']
  u_area  = u_v['A']

  return(
    list(
      v         = v,
      u_v       = u_v,
      mzopt     = mu1,
      u_mz      = u_mu1,
      cvopt     = mu2,
      u_cv      = u_mu2,
      fwhm_mz   = fwhm1,
      u_fwhm_mz = u_fwhm1,
      fwhm_cv   = fwhm2,
      u_fwhm_cv = u_fwhm2,
      area      = area,
      u_area    = u_area
    )
  )

}
getPars = function(res, dimfit){
  if (dimfit == 2)
    getPars2D(res)
  else
    getPars1D(res, dimfit)
}
fwm = function(x, ux) {
  # Weigthed mean and its uncertainty using Cochran's ANOVA for weights
  # Ref: C. Rivier et al. (2014) Accredit. Qual. Assur. 19, 269–274
  # (doi:https://doi.org/10.1007/s00769-014-1066-3)

  # u2 is the part of the variance not explained by ux (parametric unc.)
  u2 = max(0, var(x) - mean(ux ^ 2)) # Cochran's ANOVA
  ut2 = u2 + ux^2                    # Combined variance for each point

  # Variance weights
  w = (1 / ut2) / sum(1 / ut2)

  # Standard weighted mean and uncertainty
  wm  = weighted.mean(x, w)
  uwm = sqrt(1 / sum(1 / ut2))

  return(list(wm = wm, uwm = uwm))
}
#' compressMS
#'
#' Down-sample a FT-ICR MS by choosing points closest to a preset
#' regular grid.
#'
#' @param file_in
#' @param file_out
#' @param ms_type
#' @param mzMin
#' @param mzMax
#' @param dmz
#'
#' @return Nothing. Used for its side effects.
#' @export
#'
#' @examples
compressMS = function(file_in,
                      file_out,
                      ms_type  = 'fticr',
                      mzMin    = 70,
                      mzMax    = 250,
                      dmz      = 0.002,
                      test     = FALSE,
                      compMode = 'grid+targets',
                      targets  = NA,
                      dmzTarget = 1
) {

  if(ms_type == 'esquire') {
    cat('\n>>> Nothing to do ! <<<\n')

  } else {

    cat('\n>>> Processing',ms_type,'MS in file:',file_in,'\n')

    #  Get first row
    MS0 = data.table::fread(
      file = file_in ,
      header = FALSE ,
      sep = ',' ,
      stringsAsFactors = FALSE,
      nrows = 1,
      data.table = FALSE
    )
    cat('>>> Size of initial MS in memory:',
        format(object.size(MS0),units='Mb'),'\n')

    mySplit = function(x)
      as.numeric(unlist(strsplit(x, split = ' ')))

    # Build filter on 1rst MS, assuming all MS are on same grid
    ## Get mz
    mz = c()
    msLen = ncol(MS0)-8
    dbl = vapply(
      MS0[1, 9:(msLen + 8)],
      FUN = mySplit,
      FUN.VALUE = numeric(2)
    )
    mz = dbl[1, ]
    if(max(mz) < mzMin | min(mz) > mzMax) {
      message('>>> Error: mzMin/mzMax incompatible with mz range of data')
      stop(call. = FALSE)
    }
    ## Filter mz
    mz0 = mz
    selClosest = 1:length(mz)

    if(grepl('grid',compMode)) {
      # Regular grid

      mz0 = seq(
        round(max(mzMin,min(mz))),
        round(min(mzMax,max(mz))),
        by = dmz)

      # Pointer to  mz values closest to mz0
      # (timing identical to a vapply implementation)
      selClosest = vector(length=length(mz0))
      for(i in seq_along(mz0))
        selClosest[i] = which.min(abs(mz-mz0[i]))

      # Get rid of doubles
      selClosest = unique(selClosest)
    }

    if(grepl('targets',compMode) &
       !any(is.na(targets))) {
      # Targets-based grid

      mz0 = mz[selClosest]
      sel = c()
      for(targ in targets) {
        sel = c(sel,
                which(mz0 >= targ - dmzTarget &
                        mz0 <= targ + dmzTarget)
        )
      }
      # Filter and remove doubles (possible targets overlap)
      selClosest = selClosest[unique(sel)]
    }

    cat('>>> Compress m/z from',length(mz),'to',length(selClosest),
        '(',round(100*(1-length(selClosest)/length(mz))),'%)\n' )

    newRange = range(mz[selClosest])
    cat('>>> New range:',paste0(newRange,collapse='-'),'\n')

    notFinished = TRUE
    j=0
    while(notFinished) {
      j = j+1
      MS1 = c(
        MS0[1,1:6],
        paste0(newRange,collapse ='-'),
        length(selClosest),
        MS0[1,selClosest + 8]
      )
      if(j==1)
        cat('>>> Size of filtered MS in memory:',
            format(object.size(MS1),units='Mb'),'\n')
      if(test) {
        stop(call. = FALSE)
      }
      ## Save to file
      data.table::fwrite(
        MS1,
        file = file_out,
        append = j > 1,
        col.names = FALSE
      )
      cat('+')
      if(j%%20 == 0)
        cat('<',j,'\n')
      MS0 = try(
        # Get new row
        data.table::fread(
          file = file_in ,
          header = FALSE ,
          sep = ',' ,
          stringsAsFactors = FALSE,
          nrows = 1,
          skip  = j,
          data.table = FALSE
        ),
        silent = TRUE
      )
      if(class(MS0) == 'try-error')
        notFinished = FALSE
    }
    cat('\n>>> Finished <<<\n')
  }
}
