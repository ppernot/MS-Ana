# Code Setup ####
options(
  warn = 0,
  stringsAsFactors = FALSE
)

# Install packages if necessary

## CRAN packages
libs <- c('xtable','mixtools','inlmisc',
          'rlist','repmis','assertive',
          'data.table','chemCal')
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


sink(file ='./sessionInfo.txt')
print(sessionInfo(), locale=FALSE)
sink()

# Functions ####
setgPars = function() {
  list(
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
    reso     = 1200
  )
}
getPeakSpecs = function(ms_type) {

  if(! ms_type %in% c('esquire','fticr'))
    stop('>>> Bad ms_type !')

  if(ms_type == 'fticr') {

    # Acceptable range for MS peak width
    fwhm_mz_min = 0.001
    fwhm_mz_max = 0.02

    # Nominal value for peak fit
    fwhm_mz_nom = 0.01

    # Width of mz window around exact mz for signal averaging
    dmz = 5 * fwhm_mz_nom # ~ 10-sigma interval

    # Constant baseline over all matrix; median estimator
    baseline_cor = 'global+median'

  } else { # Default: Esquire

    # Acceptable range for MS peak width
    fwhm_mz_min = 0.1
    fwhm_mz_max = 0.5

    # Nominal value for peak fit
    fwhm_mz_nom = 0.2

    # Width of mz window around exact mz for signal averaging
    dmz = 5 * fwhm_mz_nom

    # Constant baseline over all matrix; median estimator
    baseline_cor = NULL

  }

  # Acceptable range for CV peak width
  fwhm_cv_min = 0.5
  fwhm_cv_max = 1.5

  # Nominal value for CV peak's FWHM constraint
  fwhm_cv_nom = 1.0

  # Width of CV window around reference CV for peak fit
  dCV = 2.5

  # Min area for validation of  peak
  area_min = 10

  return(
    list(
      fwhm_mz_min  = fwhm_mz_min,
      fwhm_mz_max  = fwhm_mz_max,
      fwhm_mz_nom  = fwhm_mz_nom,
      dmz          = dmz,
      baseline_cor = baseline_cor,
      fwhm_cv_min  = fwhm_cv_min,
      fwhm_cv_max  = fwhm_cv_max,
      fwhm_cv_nom  = fwhm_cv_nom,
      dCV          = dCV,
      area_min     = area_min
    )
  )

}
getMS = function(file, ms_type = 'esquire') {

  cat('\n>>> Reading',ms_type,'MS in file:',file,'\n')
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
bslCorMS = function(MS, baseline_cor = 'median') {
  # Baseline correction based on mode of the data distribution
  # in the hypothesis of a constant shift
  # FTICR has gamma noise distribution (doi:10.1016/j.aca.2009.10.043)
  # The median seems to be a good (constant) baseline estimator

  if(!is.null(baseline_cor)) {

    if(grepl('median',baseline_cor)) {
      cor_fun = median
    } else {
      cor_fun = mean
    }

    if(grepl('perCV',baseline_cor)) {
      cat('>>> Baseline correction of individual MSs...\n')
      # Operate correction on individual MSs
      for(i in 1:nrow(MS))
        MS[i,] = MS[i,]- cor_fun(MS[i,])
    } else {
      cat('>>> Global baseline correction...\n')
      # Global correction
      MS = MS - cor_fun(MS)
    }

  }

  return(MS)
}
getDayTag = function(tab) {
  # Scrap date from dims file name
  strsplit(
    strsplit(
      tab,
      split = ' '
    )[[1]][2],
    split = '-'
  )[[1]][1]
}
makeTag <- function(CVTable, msTable, userTag) {
  date = getDayTag(CVTable)
  paste0(
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
  expandFactor = 5,
  gPars
) {

  nCV = length(CV)

  # Expose fitOut list
  for (n in names(fitOut))
    assign(n,rlist::list.extract(fitOut,n))

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

  if(type == 'CV') {
    CVp = CV
    MSp = MS
  }else{
    # Fill matrix CV-wise for better image
    CVp = seq(min(CV),max(CV),length.out=expandFactor*length(CV))
    MSp = matrix(0,nrow=length(CVp),ncol=ncol(MS))
    for(i in 1:length(CV)) {
      iCV = which.min(abs(CVp-CV[i])) # Closest point
      CVp[iCV] = CV[i]
      MSp[iCV,] = MS[i,]
    }
  }

  mzlim  = c(mz1,mz2)
  CVlim  = range(CV)
  if(type == 'CV')
    CVlimf = range(CVf)
  else
    CVlimf = range(CVp)

  ## 1. image
  sel1 = apply(cbind(CVp-CVlim[1],CVp-CVlim[2]),1,prod) <= 0
  sel2 = apply(cbind(mz-mzlim[1],mz-mzlim[2]),1,prod) <= 0

  image(
    CVp[sel1], mz[sel2], MSp[sel1,sel2],
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
      rect(CVlimf[1], mz1, CVlimf[2], mz2,
           col = cols_tr[4], border = NA)
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
    } else {
      abline(h = v[1], lty = 2, col = cols[2])
      abline(v = v[3], lty = 2, col = cols[2])
    }
    if(!is.na(mz0))
      abline(h = mz0, lty = 2, col = cols[2])
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
        x=CVlim[1],
        y=ylim[2]-0.75*strheight(val,units='user',cex=0.8),
        yjust = 0,
        # title  = val,
        legend = val,
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
      legend(
        x = mzlim[1],
        y = ylim[2] - strheight(val, units = 'user', cex = 0.8),
        yjust = 0,
        legend = val,
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
trapz = function (x, y) {
  if(length(x) == 0)
    return(0)
  idx = 2:length(x)
  return(
    as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
  )
}
fit1D_MS <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  fwhm_mz_nom,
  weighted = FALSE
) {

  # Select CV0
  iCV = which(CV >= CV0)[1]

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # mz profile
  MStot = MS[iCV,]
  mzl = as.vector(mz[selMz])
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
  sigma0 = fwhm_mz_nom/2.355
  A0 = sqrt(2*pi) * sigma0 * max(MSl)
  lower = c(
    mu    = mu - dmz,
    sigma = 0.1 * sigma0,
    A     = 0.5 * A0
  )
  start = c(
    mu    = mu,
    sigma = sigma0,
    A     = A0
  )
  upper = c(
    mu    = mu + dmz,
    sigma = 2.0 * sigma0,
    A     = 1.5 * A0
  )

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
fit1D <- function(
  mz0, CV0,
  dmz, dCV,
  mz, CV, MS,
  fwhm_cv_nom,
  weighted = FALSE,
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
  sigma0 = fwhm_cv_nom/2.355
  A0 = sqrt(2*pi) * sigma0 * max(mMS)
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
  fwhm_mz_nom,
  fwhm_cv_nom,
  weighted = FALSE,
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
  x = as.vector(grid$x)
  y = as.vector(grid$y)
  z = as.vector(t(MSloc))

  weights = rep(1,length(z))
  if(weighted){ # Poisson
    weights = 1 / z
    vm = min(z[z>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  # First pass
  maxz = which.max(z)
  mx = x[maxz]
  sx0 = fwhm_mz_nom / 2.355
  sy0 = fwhm_cv_nom / 2.355
  A0  = 2 * pi * sx0 * sy0 * max(z)
  lower = c(
    mx  = mx - dmz / 2,
    sx  = sx0 / 4,
    my  = CV0 - dCV / 10,
    sy  = 0.8 * sy0,
    A   = 0.5 * A0
  )
  start = c(
    mx  = mx,
    sx  = sx0,
    my  = CV0,
    sy  = sy0,
    A   = A0
  )
  upper = c(
    mx  = mx + dmz / 2,
    sx  = 4 * sx0,
    my  = CV0 + dCV / 10,
    sy  = 1.2 * sy0,
    A   = 1.5 * A0
  )

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
      cat('>>> Grid selection\n')
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
      cat('>>> Targets selection\n')

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

    # Reorder to compensate for targets +/- random order
    selClosest = sort(selClosest)

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
        stop('>>> Stopped by test=TRUE',call. = FALSE)
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
# From discontinued FTICRMS package
# install.packages("~/Téléchargements/FTICRMS_0.8.tar.gz",
#                 repos = NULL, type = "source")
# Does not run...
baselineMS = function (spect, init.bd,
                       sm.par = 1e-11,
                       sm.ord = 2,
                       max.iter = 20,
                       tol = 5e-08,
                       sm.div = NA,
                       sm.norm.by = c(
                         "baseline",
                         "overestimate",
                         "constant"),
                       neg.div = NA,
                       neg.norm.by = c(
                         "baseline",
                         "overestimate",
                         "constant"),
                       rel.conv.crit = TRUE,
                       zero.rm = TRUE,
                       halve.search = FALSE)
{
    require(Matrix)
    L <- length(spect)
    if (zero.rm && any(spect == 0)) {
      wh <- which(spect == 0)
      warning(paste(length(wh), "entries are equal to zero; replacing them with average value of surrounding points."))
      wh <- by(wh, cumsum(c(1, diff(wh) > 1)), c)
      for (i in wh) {
        if (1 %in% i) {
          spect[i] <- spect[max(i) + 1]
        }
        else if (L %in% i) {
          spect[i] <- spect[min(i) - 1]
        }
        else {
          spect[i] <- mean(spect[range(i) + c(-1, 1)])
        }
      }
    }
    neg.norm.by <- match.arg(neg.norm.by)
    sm.norm.by <- match.arg(sm.norm.by)
    if (sm.norm.by == "constant" || neg.norm.by == "constant") {
      bin.ends <- round((0:1024)/1024 * L)
      ss <- sapply(1:1024, function(x) {
        .biweight.FTICRMS(spect[(bin.ends[x] + 1):bin.ends[x + 1]],
                          K = 9)$scale/qnorm(0.75)
      })
      ss <- .biweight.FTICRMS(ss, K = 9)$center
    }
    if (is.na(sm.div)) {
      sm.div <- switch(sm.norm.by, baseline = 0.5223145,
                       overestimate = 1,
                       constant = ss)
    }
    sm.fact <- L^(2 * sm.ord) * sm.par/sm.div
    if (is.na(neg.div)) {
      neg.div <- switch(neg.norm.by, baseline = 0.4210109,
                        overestimate = 1, constant = ss/sqrt(pi/2))
    }
    if (missing(init.bd)) {
      bd <- rep(median(spect, na.rm = TRUE), L)
    }
    else {
      bd <- init.bd
    }
    sm.ord.even <- sm.ord - (sm.ord%%2)
    L1 <- L - sm.ord.even
    M <- t(
      new(
        "dgCMatrix",
        Dim = as.integer(c(L, L1)),
        i = as.integer(outer(0:sm.ord.even, 0:(L1 - 1), "+")),
        p = as.integer(cumsum(c(
          0, rep(sm.ord.even + 1, L1)
        ))),
        x = as.numeric(rep(choose(sm.ord.even, 0:sm.ord.even) *
                             (-1) ^
                             (0:sm.ord.even), L1))))
    if (sm.ord %% 2 == 1) {
      M <-
        new(
          "dgCMatrix",
          Dim = as.integer(c(L1, L1)),
          i = as.integer(c(0,
                           1, outer(
                             c(0, 2), 0:(L1 - 3), "+"
                           ), L1 - 2, L1 -
                             1)),
          p = as.integer(cumsum(c(0, rep(
            2, L1
          )))),
          x = as.numeric(c(-1, -0.5, 1, rep(c(
            -0.5, 0.5
          ), L1 -
            3), -1, 0.5, 1))
        ) %*% M
    }
    if (sm.norm.by == "constant") {
      Mtmp <- sm.fact * crossprod(M)
      rm(M)
      tmpgc <- ""
      while (!identical(tmpgc, tmpgc <- gc())) {
      }
    }
    indicator <- (bd > spect)
    indicator[is.na(indicator)] <- FALSE
    indicator0 <- rep(TRUE, L)
    bd0 <- Inf
    changed <- c()
    iter <- 0
    if (halve.search) {
      hs <- c()
    }
    else {
      hs <- NA
    }
    if (rel.conv.crit) {
      step.div <- bd
    }
    else {
      step.div <- 1
    }
    while (mean(((bd - bd0) / step.div) ^ 2) > tol &&
           iter < max.iter) {
      indicator0 <- indicator
      bd0 <- bd
      if (sm.norm.by == "baseline") {
        Mtmp <- crossprod(
          new(
            "dsCMatrix",
            Dim = as.integer(c(L1,
                               L1)),
            x = as.numeric(sqrt(sm.fact / bd[floor(1 +
                                                     (L - L1) /
                                                     2):floor((L + L1) / 2)])),
            p = as.integer(0:L1),
            i = as.integer(0:(L1 - 1))
          ) %*% M
        )
      }
      else if (sm.norm.by == "overestimate") {
        Mtmp <- crossprod(
          new(
            "dsCMatrix",
            Dim = as.integer(c(L1,
                               L1)),
            x = as.numeric(sqrt(
              sm.fact / ifelse(bd >
                                 spect, as.numeric(bd - spect), 1)[floor(1 + (L -
                                                                                L1) /
                                                                           2):floor((L + L1) / 2)]
            )),
            p = as.integer(0:L1),
            i = as.integer(0:(L1 - 1))
          ) %*% M
        )
      }
      if (neg.norm.by == "baseline") {
        B <- Mtmp %*% bd - 1/2 +
          ifelse(is.na(spect), 0,
                 as.numeric((bd - spect) * indicator /
                              bd / neg.div))
        M1 <- Mtmp + new(
          "dsCMatrix",
          Dim = as.integer(c(L,
                             L)),
          x = as.numeric(indicator / bd / neg.div),
          p = as.integer(0:L),
          i = as.integer(0:(L - 1))
        )
      }
      else if (neg.norm.by == "overestimate") {
        B <- Mtmp %*% bd - 1 / 2 + ifelse(is.na(spect), 0,
                                          as.numeric(indicator / neg.div))
        M1 <- Mtmp + new(
          "dsCMatrix",
          Dim = as.integer(c(L,
                             L)),
          x = as.numeric(
            indicator / neg.div / ifelse(as.numeric(bd) >
                                           spect, as.numeric(bd - spect), 1)
          ),
          p = as.integer(0:L),
          i = as.integer(0:(L - 1))
        )
      }
      else if (neg.norm.by == "constant") {
        M1 <- Mtmp + new(
          "dsCMatrix",
          Dim = as.integer(c(L,
                             L)),
          x = as.numeric(indicator / neg.div),
          p = as.integer(0:L),
          i = as.integer(0:(L - 1))
        )
        B <- Mtmp %*% bd - 1 / 2 + ifelse(is.na(spect), 0,
                                          as.numeric((bd - spect) * indicator /
                                                       neg.div))
      }
      if (sm.ord > 2) {
        M1@factors <- list(spdCholesky = Cholesky(M1, Imult = 1e-17))
      }
      bd.step <- as.numeric(solve(M1, B))
      curr.val <-
        sum(bd) - sum(bd * (Mtmp %*% bd)) - sum(indicator *
                                                  (bd - spect))
      if (rel.conv.crit) {
        step.div <- bd
      }
      else {
        step.div <- 1
      }
      if (halve.search) {
        hs <- c(hs, 0)
        while (sum((bd.step / step.div) ^ 2) > tol && sum((bd -
                                                           bd.step)) - sum((bd - bd.step) * (Mtmp %*% (bd -
                                                                                                       bd.step))) - sum(indicator * ((bd - bd.step) -
                                                                                                                                     spect)) <= curr.val) {
          bd.step <- bd.step / 2
          hs[length(hs)] <- hs[length(hs)] + 1
        }
      }
      bd <- bd - bd.step
      indicator <- ifelse(is.na(spect), 0, as.numeric(bd >
                                                        spect))
      changed <- c(changed, sum(xor(indicator, indicator0)))
      iter <- iter + 1
    }
    if (mean(((bd - bd0) / step.div) ^ 2) > tol) {
      warning(
        paste(
          "Iteration limit of",
          max.iter,
          "reached without convergence to specified tolerance."
        )
      )
    }
    list(
      baseline = as.numeric(bd),
      iter = iter,
      changed = changed,
      hs = hs
    )
}

dmsAnalysis = function(
  ms_type         = c('esquire','fticr')[1],
  taskTable       = NA,
  tgTable         = NA,
  dataRepo        = '../data/',
  figRepo         = '../results/figs/',
  tabRepo         = '../results/tables/',
  fit_dim         = 1,
  filter_results  = TRUE,
  userTag         = paste0('fit_dim_',fit_dim),
  save_figures    = TRUE,
  plot_maps       = FALSE,
  fallback        = TRUE,
  correct_overlap = FALSE,
  weighted_fit    = FALSE,
  refine_CV0      = TRUE,
  debug           = FALSE
) {
  # Changes ----
  #
  # 2020_07_16 [PP]
  # - Added 2 new colums to results table
  #   to store m/z and u_m/z from 2D fit
  # - Integrate fast method (fit_dim = 0)
  # - Added data path management
  # - Added optional tag ('userTag'):
  #   presently defined by 'fit_dim' value,
  #   to avoid overwriting of results files
  #   when trying different 'fit_dim' options
  # 2020_07_20 [PP]
  # - Replaced '=' by '_' in userTag (Windows pb.)
  # - Reparameterized gaussians with area replacing height
  #   (avoids covariances in estimation of u_area )
  # - Suppressed 'rho' param in 2D gaussians
  # 2020_07_21 [PP]
  # - Corrected typo in calculation of u_area in 'getPars1D()'
  # 2020_07_24 [PP]
  # - Save ctrl params as metadata
  # 2020_09_24 [PP]
  # - All input files should now be comma-delimited
  # 2021_03_15 [PP]
  # - Adapted code to FT-ICR MS files:
  #   * new 'ms_type' flag
  #   * new 'getMS()' function to handle ms_type cases
  #   * new 'trapz()' function to handle MS integration
  #     with irregular mz grids
  #   * adapted 'fit1D()' and 'fit2D()' to use 'trapz()'
  # 2021_04_06 [PP]
  # - Added baseline correction function
  #   * new 'baseline_cor' flag
  #   * new 'bslCorMS()' function
  # 2021_04_07 [PP]
  # - Changed management of peak specifications to facilitate
  #   modifications of apparatus characteristics
  #   * new 'getPeakSpecs()' functions
  #   * changed 'const_fwhm' to 'fwhm_cv_nom'
  #   * new 'fwhm_mz_nom' variable
  #   * changed logic on fwhm fit constraints
  #     Fit is now always constrained in fit1D_MS(),
  #     fit_1D() and fit_2D()
  #
  #---

  # Graphical params for external and local plots
  gPars = setgPars()
  gParsLoc = gPars
  gParsLoc$cex = 1
  gParsLoc$lwd = 1.5

  # Get apparatus-dependent specifications ----
  peakSpecs   = getPeakSpecs(ms_type)

  fwhm_mz_min = peakSpecs$fwhm_mz_min
  fwhm_mz_max = peakSpecs$fwhm_mz_max
  fwhm_mz_nom = peakSpecs$fwhm_mz_nom
  dmz         = peakSpecs$dmz

  fwhm_cv_min = peakSpecs$fwhm_cv_min
  fwhm_cv_max = peakSpecs$fwhm_cv_max
  fwhm_cv_nom = peakSpecs$fwhm_cv_nom
  dCV         = peakSpecs$dCV

  baseline_cor = peakSpecs$baseline_cor
  area_min     = peakSpecs$area_min

  # Gather run params for reproducibility ----
  ctrlParams = list(
    userTag        = userTag,
    ms_type        = ms_type,
    taskTable      = taskTable,
    tgTable        = tgTable,
    dataRepo       = dataRepo,
    filter_results = filter_results,
    fwhm_mz_min    = fwhm_mz_min,
    fwhm_mz_max    = fwhm_mz_max,
    fwhm_mz_nom    = fwhm_mz_nom,
    dmz            = dmz,
    fwhm_cv_min    = fwhm_cv_min,
    fwhm_cv_max    = fwhm_cv_max,
    fwhm_cv_nom    = fwhm_cv_nom,
    dCV            = dCV,
    area_min       = area_min,
    fit_dim        = fit_dim,
    fallback       = fallback,
    weighted_fit   = weighted_fit,
    refine_CV0     = refine_CV0,
    baseline_cor   = baseline_cor
  )

  # Check sanity of parameters ----
  assertive::assert_all_are_existing_files(dataRepo)
  assertive::assert_all_are_existing_files(figRepo)
  assertive::assert_all_are_existing_files(tabRepo)

  file = paste0(dataRepo, tgTable)
  assertive::assert_all_are_existing_files(file)

  file = paste0(dataRepo, taskTable)
  assertive::assert_all_are_existing_files(file)

  assertive::assert_is_numeric(fwhm_mz_min)
  if(!assertive::is_positive(fwhm_mz_min))
    stop(paste0('Erreur: fwhm_mz_min =',
                fwhm_mz_min,' should be positive'))

  assertive::assert_is_numeric(fwhm_mz_max)
  if(!assertive::is_positive(fwhm_mz_max))
    stop(paste0('Erreur: fwhm_mz_max =',
                fwhm_mz_max,' should be positive'))

  assertive::assert_is_numeric(fwhm_cv_min)
  if(!assertive::is_positive(fwhm_cv_min))
    stop(paste0('Erreur: fwhm_cv_min =',
                fwhm_cv_min,' should be positive'))

  assertive::assert_is_numeric(fwhm_cv_max)
  if(!assertive::is_positive(fwhm_cv_max))
    stop(paste0('Erreur: fwhm_cv_max =',
                fwhm_cv_max,' should be positive'))

  assertive::assert_is_numeric(area_min)
  if(!assertive::is_positive(area_min))
    stop(paste0('Erreur: area_min =',area_min
                ,' should be positive'))

  assertive::assert_is_numeric(dmz)
  if(!assertive::is_positive(dmz))
    stop(paste0('Erreur: dmz =',dmz,' should be positive'))

  assertive::assert_is_numeric(dCV)
  if(!assertive::is_positive(dCV))
    stop(paste0('Erreur: dCV =',dCV,' should be positive'))

  # Get targets ----
  targets = readTargetsFile(paste0(dataRepo, tgTable))
  empty = rep(NA,nrow(targets))
  if(!'CV_ref' %in% colnames(targets))
    targets = cbind(targets,CV_ref=empty)

  # Get tasks list ----
  Tasks = readTasksFile(paste0(dataRepo, taskTable))

  # Check that files exist before proceeding
  if('path' %in% colnames(Tasks)) {
    files = paste0(dataRepo,Tasks[,'path'],Tasks[,'MS_file'])
  } else {
    files = paste0(dataRepo,Tasks[,'MS_file'])
  }
  assertive::assert_all_are_existing_files(files)

  files = paste0(dataRepo,Tasks[,'DMS_file'])
  assertive::assert_all_are_existing_files(files)

  # Loop over tasks ----
  dilu = NA
  for(task in 1:nrow(Tasks)) {

    msTable = Tasks[task,'MS_file']
    CVTable = Tasks[task,'DMS_file']
    if('dilu' %in% colnames(Tasks))
      dilu    = Tasks[task,'dilu']
    dataPath = ''
    if('path' %in% colnames(Tasks))
      if(!is.na(Tasks[task,'path']))
        dataPath = Tasks[task,'path']

    # Build tag
    tag  = makeTag(CVTable, msTable, userTag)

    ## Get MS ----
    file = paste0(dataRepo, dataPath, msTable)
    lMS  = getMS(file, ms_type)
    time = lMS$time
    mz   = lMS$mz
    MS   = lMS$MS
    rm(lMS) # Clean-up memory

    # Get CV ----
    file = paste0(dataRepo, CVTable)
    CV0 = read.table(
      file = file,
      header = FALSE,
      sep = '\t',
      stringsAsFactors = FALSE
    )
    CV = rev(CV0[, 4]) # We want increasing CVs

    ## Ensure CV & MS tables conformity
    t0   = Tasks[task,'t0']
    it   = which.min(abs(time - t0))
    selt = which(time >= time[it])
    nt   = length(selt)

    CV0   = Tasks[task,'CV0']
    iCV   = which.min(abs(CV - CV0))
    selCV = which(CV <= CV[iCV])
    nCV   = length(selCV)

    ncut  = min(nt,nCV)
    selt  = selt[1:ncut]
    selCV = rev(rev(selCV)[1:ncut])

    time = time[selt]
    MS   = MS[selt,]
    MS   = apply(MS, 2, rev) # reverse column to conform with CV
    CV   = CV[selCV]
    nCV  = length(CV)

    # Baseline correction ----
    MS = bslCorMS(MS, baseline_cor)

    ## Initialize results table
    resu = cbind(
      targets,
      empty,empty,
      empty,empty,
      empty,empty,
      empty,empty,
      empty,empty,
      empty,empty,
      empty)
    colnames(resu) = c(
      colnames(targets),
      'm/z',     'u_m/z',
      'CV',      'u_CV',
      'FWHM_m/z','u_FWHM_m/z',
      'FWHM_CV', 'u_FWHM_CV',
      'Area',    'u_Area',
      'fit_dim',  'dilu',
      'tag'
    )

    if( fit_dim == 0) {
      xic = matrix(mz,ncol=1)
      colnames(xic) = 'm/z'
      xfi = matrix(mz,ncol=1)
      colnames(xfi) = 'm/z'
    } else {
      xic = cbind(time,rev(CV))
      colnames(xic) = c('time','CV')
      xfi = cbind(time,rev(CV))
      colnames(xfi) = c('time','CV')
    }

    # Loop over targets ----
    for( it in 1:nrow(targets) ) {

      mz0 = targets[it,'m/z_ref']
      CV0 = targets[it,'CV_ref']

      if(fit_dim == 2) {
        # 2D fit of peaks
        fitOut = fit2D(
          mz0, CV0,
          dmz, dCV,
          mz, CV, MS,
          fwhm_mz_nom = fwhm_mz_nom,
          fwhm_cv_nom = fwhm_cv_nom,
          weighted = weighted_fit,
          refine_CV0 = refine_CV0,
          correct_overlap = correct_overlap
        )
        dimfit = 2
        if(class(fitOut$res) == 'try-error' & fallback) {
          # 1D fit of peaks
          fitOut = fit1D(
            mz0, CV0,
            dmz, dCV,
            mz, CV, MS,
            fwhm_cv_nom = fwhm_cv_nom,
            weighted = weighted_fit,
            refine_CV0 = refine_CV0,
            correct_overlap = correct_overlap
          )
          dimfit = 1
        }

      } else if (fit_dim == 1) {
        # 1D fit of peaks; fixed m/z
        fitOut = fit1D(
          mz0, CV0,
          dmz, dCV,
          mz, CV, MS,
          fwhm_cv_nom = fwhm_cv_nom,
          weighted = weighted_fit,
          refine_CV0 = refine_CV0,
          correct_overlap = correct_overlap
        )
        dimfit = 1

      } else {
        # 1D m/z fit of peak; fixed CV
        fitOut = fit1D_MS(
          mz0, CV0,
          dmz, dCV,
          mz, CV, MS,
          weighted = weighted_fit,
          fwhm_mz_nom = fwhm_mz_nom
        )
        dimfit = 0
      }

      for (n in names(fitOut))
        assign(n,rlist::list.extract(fitOut,n))

      targets[it,'m/z_ref'] = mz0

      if(class(res)=="try-error") {
        # Fit failed => no fit params
        v       = NA
        mzopt   = NA
        cvopt   = NA
        fwhm_mz = NA
        fwhm_cv = NA
        area    = NA
        warning = TRUE

      } else {
        v   = summary(res)$parameters[,"Estimate"]
        peakPars = getPars(res,dimfit)
        for (n in names(peakPars))
          assign(n,rlist::list.extract(peakPars,n))

        # Quality control
        if(filter_results &
           (
             ifelse(
               !is.na(fwhm_cv),
               fwhm_cv <= fwhm_cv_min | fwhm_cv >= fwhm_cv_max,
               FALSE
             ) |
             ifelse(
               !is.na(fwhm_mz),
               fwhm_mz <= fwhm_mz_min | fwhm_mz >= fwhm_mz_max,
               FALSE
             ) |
             area <= area_min
           )
        ) {
          warning = TRUE
          # Do not store results

        } else {
          warning = FALSE
          # Store in results table
          resu[it,'m/z']        = signif(mzopt,6)
          resu[it,'u_m/z']      = signif(u_mz,2)
          resu[it,'CV']         = signif(cvopt,4)
          resu[it,'u_CV']       = signif(u_cv,2)
          resu[it,'FWHM_m/z']   = signif(fwhm_mz,3)
          resu[it,'u_FWHM_m/z'] = signif(u_fwhm_mz,2)
          resu[it,'FWHM_CV']    = signif(fwhm_cv,3)
          resu[it,'u_FWHM_CV']  = signif(u_fwhm_cv,2)
          resu[it,'Area']       = signif(area,3)
          resu[it,'u_Area']     = signif(u_area,2)
        }
      }
      resu[it,'fit_dim'] = dimfit
      resu[it,'dilu'] = dilu
      resu[it,'tag'] = tag

      # Plot data and fit results
      pars = paste0(
        ifelse (warning, '** WARNING **\n','') ,
        ifelse(dimfit == 1,
               '',
               paste0('m/z = ', signif(mzopt,6),'\n')),
        ifelse(dimfit == 0,
               '',
               paste0('CV = ', signif(cvopt,4),'\n')),
        ifelse(dimfit ==1,
               '',
               paste0('FWHM_m/z = ', signif(fwhm_mz,3),'\n')),
        ifelse(dimfit ==0,
               '',
               paste0('FWHM_CV = ', signif(fwhm_cv,3),'\n')),
        'Area = ', signif(area,3)
      )

      if (class(fitOut$res) != "try-error")
        print(coefficients(fitOut$res))

      plotPeak(
        mz, CV, MS,
        fitOut,
        mex = targets[it,'m/z_ref'],
        leg = targets[it,'Name'],
        tag = tag,
        val = pars,
        type = ifelse(dimfit==0,'m/z','CV'),
        CV0 = CV0,
        gPars = gParsLoc
      )

      if(save_figures) {
        png(filename = paste0(figRepo, tag, '_', targets[it, 1], '.png'),
            width    = 2*gPars$reso,
            height   =   gPars$reso )
        plotPeak(
          mz, CV, MS,
          fitOut,
          mex = targets[it,'m/z_ref'],
          leg = targets[it,'Name'],
          tag = tag,
          val = pars,
          type = ifelse(dimfit==0,'m/z','CV'),
          CV0 = CV0,
          gPars = gPars
        )
        dev.off()
      }

      # Save XIC and fit file
      nam0 = colnames(xic)
      if(fit_dim == 0) {
        xic = cbind(xic,mMStot)
        fit = peak_shape(mz, v)
        xfi = cbind(xfi,fit)
      } else {
        xic = cbind(xic,rev(mMStot))
        fit = peak_shape(CV, v)
        xfi = cbind(xfi,rev(fit))
      }
      colnames(xic) = c(nam0,targets[it,1])
      colnames(xfi) = c(nam0,targets[it,1])

      # if(debug) break
    }

    # Global Heat maps
    if(plot_maps) {
      mex = targets[,'m/z_ref']
      plotMaps(
        mz, CV, MS,
        mex = mex,
        leg = 'log10(MS)',
        tag = tag,
        mzlim = c(min(mex)-5*dmz,max(mex)+5*dmz),
        CVlim = range(CV),
        logz = TRUE,
        gPars = gParsLoc
      )
      if(save_figures) {
        png(
          filename = paste0(figRepo, tag,
                            '_heatmaps.png'),
          width    = 2*gPars$reso,
          height   =   gPars$reso )
        plotMaps(
          mz, CV, MS,
          mex = mex,
          leg = 'log10(MS)',
          tag = tag,
          mzlim = c(min(mex)-5*dmz,max(mex)+5*dmz),
          CVlim = range(CV),
          logz = TRUE,
          gPars = gPars
        )
        dev.off()
      }
    }

    # Save results
    write.csv(resu,
              file = paste0(tabRepo, tag, '_results.csv'),
              row.names = FALSE)
    write.csv(xic,
              file  = paste0(tabRepo, tag, '_XIC.csv'),
              row.names = FALSE)
    write.csv(xfi,
              file  = paste0(tabRepo, tag, '_fit.csv'),
              row.names = FALSE)

    # Metadata
    rlist::list.save(
      ctrlParams,
      file = file.path(
        tabRepo,
        paste0(tag,'_ctrlParams.yaml')
      )
    )

    if(debug){
      message('Ended prematurely (debug)...')
      stop(call. = FALSE)
    }

    # End of targets loop ----
  }

  # End of tasks loop ----
}
#' Truncate value and uncertainty to consistent number of digits.
#'
#' @param y (numeric) value
#' @param uy (numeric) uncertainty on `y`
#' @param numDig (numeric) number of digits to keep on `uy`
#'
#' @return A list with strings of truncated values of `y` and `uy`.
#' @export
#'
#' @examples
formatUnc = function(y, uy, numDig = 2) {

  if (!is.finite(y) | !is.finite(uy) | uy <= 0)
    return(
      list(y  = y, uy = uy)
    )

  # Get scales
  n0 = 1 + floor(log10(abs(y)))
  n1 = floor(log10(uy))

  # Format uncertainty
  fmt = switch(
    sign(n1) + 2, # Map (-1,0,1) to (1,2,3)
    paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f"),
    paste0("%", n0 - n1 + numDig - 1, ".", -n1 + numDig - 1, "f"),
    paste0("%", n0, ".0f")
  )
  short_y  = sprintf(fmt, y)
  short_uy = sprintf(fmt,uy)

  return(
    list(
      y  = short_y,
      uy = short_uy
    )
  )
}
formatUncVec = function(y, uy, numDig = 2) {
  ftab = matrix(NA, nrow = length(y), ncol = 2)
  colnames(ftab) = c('y', 'uy')
  for (i in seq_along(y)) {
    f = formatUnc(y[i], uy[i], numDig)
    ftab[i, 1] = as.numeric(f$y)
    ftab[i, 2] = as.numeric(f$uy)
  }
  return(ftab)
}
cv = function(x, ux) {
  if(x==0)
    return(NA)
  else
    return(100 * abs(ux / x))
}
lodCal = function(
  D,
  serialize = TRUE,
  weighted_mean = TRUE,
  checkDin = FALSE
) {

  xall    = D$x
  yall    = D$y

  fitList = list()

  if(serialize) {
    dayTags = D$t
    series  = unique(dayTags)

    # Treat data per day
    tabLod = tabLodu = tabLoq  = tabLoqu  = c()
    tabInt = tabSlop = tabIntu = tabSlopu = tabR2 = c()

    for(iS in seq_along(series)) {

      # Select data of the day
      sel = dayTags == series[iS]
      xd  = xall[sel]
      yd  = yall[sel]

      # Average over daily repeats
      yx = split(yd,xd)
      x  = unique(xd)
      y  = sapply(yx, mean)

      # Linear fit
      reg = lm(y ~ x)
      fitList[[series[iS]]] = reg
      linePars     = summary(reg)$coefficients[, 1]
      ulinePars    = summary(reg)$coefficients[, 2]
      tabR2[iS]    = summary(reg)$r.squared
      tabInt[iS]   = linePars[1]
      tabIntu[iS]  = ulinePars[1]
      tabSlop[iS]  = linePars[2]
      tabSlopu[iS] = ulinePars[2]

      p = predict(
        reg,
        newdata = data.frame(x=0), # in case x=0 is absent from data
        interval = 'prediction',
        level = 0.90)
      S90 = p[1, 3] - p[1, 1]

      slope  = linePars[2]
      uSlope = ulinePars[2]
      tabLod[iS]  = 2 * S90 / slope
      tabLodu[iS] = 2 * S90 / slope ^ 2 * uSlope
      tabLoq[iS]  = 3.04 * S90 / slope
      tabLoqu[iS] = 3.04 * S90 / slope ^ 2 * uSlope

      if(checkDin) {
        din = chemCal::lod(reg,method ='din')$x
        if(abs((tabLod[iS] - din)/din) > 1e-3)
          warning('LOD is different from reference value...')
      }
    }

    # Average daily LODs/LOQs
    if (weighted_mean) {
      # Weighted mean
      # tabw  = 1 / tabLodu ^ 2
      # tabw  = tabw / sum(tabw)
      # LOD   = sum(tabLod * tabw)
      # LOD.u = sqrt(sum(tabw * (tabLod ^ 2 + tabLodu ^ 2)) - LOD ^ 2)
      tmp    = fwm(tabLod, tabLodu)
      LOD    = tmp$wm
      LOD.u  = tmp$uwm
      LOD.cv = cv(LOD,LOD.u)

      # tabw  = 1 / tabLoqu ^ 2
      # tabw  = tabw / sum(tabw)
      # LOQ   = sum(tabLoq * tabw)
      # LOQ.u = sqrt(sum(tabw * (tabLoq ^ 2 + tabLoqu ^ 2)) - LOQ ^ 2)
      tmp    = fwm(tabLoq, tabLoqu)
      LOQ    = tmp$wm
      LOQ.u  = tmp$uwm
      LOQ.cv = cv(LOQ,LOQ.u)

    } else {
      # Arithmetic mean
      LOD   = mean(tabLod)
      LOD.u = sd(tabLod) / sqrt(length(tabLod))
      LOD.cv = cv(LOD,LOD.u)

      LOQ   = mean(tabLoq)
      LOQ.u = sd(tabLoq) / sqrt(length(tabLoq))
      LOQ.cv = cv(LOQ,LOQ.u)

    }

    tmp     = fwm(tabInt,tabIntu)
    Intm    = tmp$wm
    Intm.u  = tmp$uwm
    Intm.cv = cv(Intm,Intm.u)
    tmp     = fwm(tabSlop,tabSlopu)
    Slom    = tmp$wm
    Slom.u  = tmp$uwm
    Slom.cv = cv(Slom,Slom.u)

    tabres = data.frame(
      Int   = c(tabInt,Intm,Intm.u,Intm.cv),
      Int.u = c(tabIntu,NA,NA,NA),
      Slo   = c(tabSlop,Slom,Slom.u,Slom.cv),
      Slo.u = c(tabSlopu,NA,NA,NA),
      R2    = c(tabR2,NA,NA,NA),
      LOD   = c(tabLod,LOD,LOD.u,LOD.cv),
      LOD.u = c(tabLodu,NA,NA,NA),
      LOQ   = c(tabLoq,LOQ,LOQ.u,LOQ.cv),
      LOQ.u = c(tabLoqu,NA,NA,NA)
    )
    tabres = cbind(
      Series = c(series,'wMean','wSE','CV'),
      tabres)

  } else {

    series ='fitMean'

    # Unweighted means
    yx = split(yall,xall)
    x  = unique(xall)
    y  = sapply(yx, mean)

    # Linear fit
    reg = lm(y ~ x)
    fitList[[series]] = reg
    linePars     = summary(reg)$coefficients[, 1]
    ulinePars    = summary(reg)$coefficients[, 2]
    R2     =  summary(reg)$r.squared
    Int    = linePars[1]
    Intu   = ulinePars[1]
    slope  = linePars[2]
    uSlope = ulinePars[2]

    p = predict(
      reg,
      newdata = data.frame(x=0), # in case x=0 is absent from data
      interval = 'prediction',
      level = 0.90)
    S90 = p[1, 3] - p[1, 1]

    LOD    = 2 * S90 / slope
    LOD.u  = 2 * S90 / slope ^ 2 * uSlope
    LOQ    = 3.04 * S90 / slope
    LOQ.u  = 3.04 * S90 / slope ^ 2 * uSlope

    if(checkDin) {
      din = chemCal::lod(reg,method ='din')$x
      if(abs((LOD - din)/din) > 1e-3)
        warning('LOD is different from reference value...')
    }
    tabres = data.frame(
      Series= series,
      Int   = Int,
      Int.u = Intu,
      Slo   = slope,
      Slo.u = uSlope,
      R2    = R2,
      LOD   = LOD,
      LOD.u = LOD.u,
      LOQ   = LOQ,
      LOQ.u = LOQ.u
    )
  }

  return(
    list(
      LOD    = LOD,
      LOD.u  = LOD.u,
      LOQ    = LOQ,
      LOQ.u  = LOQ.u,
      tab    = tabres,
      series = series,
      fit    = fitList
    )
  )
}



