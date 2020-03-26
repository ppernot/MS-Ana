# User configuration params ####
taskTable = 'list_of_files_Viet.csv'
tgTable = 'list_of_targets_Plasma_Viet.csv'
filter_results = TRUE
fwhm_min = 0.5
fwhm_max = 1.5
area_min = 10
save_figures = TRUE
weighted_fit = FALSE
fit_dim = 2     # Fit 2D peaks
fallback = TRUE # Fallback on 1D fit if 2D fails
dmz = 1.0       # Width of mz window around
                # exact mz for signal averaging
dCV = 2.0       # Width of CV window around
                # reference CV for peak fit

# Code Setup ####
# Install packages if necessary
## CRAN packages
libs <- c('xtable','mixtools','inlmisc',
          'rlist','repmis','assertive')
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE,
      repos = 'https://cran.univ-paris1.fr'
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
gParsLoc$lwd = 2

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

## Functions
peak_shape = function(x,p) {
  if(length(p)==3)
    p[3]*exp(-1/2*(x-p[1])^2/p[2]^2)
  else
    sqrt(2*pi)*p['k']*abs(p['sx']) *
    exp(-1/2*(x-p['my'])^2/p['sy']^2)
}
plotPeak = function(
  mz, CV, MS, #CVf, mMS, vg,
  fitOut,
  mex = NA,
  leg = NA,
  val = NA,
  tag = NA,
  gPars
) {

  nCV = length(CV)

  # Expose fitOut list
  for (n in names(fitOut))
    assign(n,rlist::list.extract(fitOut,n))

  mzlim = c(mz1,mz2)
  CVlim = range(CV)

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
      # image(CVf, mzLoc, MSloc)
      p = matrix(
        predict(res),
        ncol=ncol(MSloc),
        nrow=nrow(MSloc),
        byrow = TRUE
      )
      contour(CVf, mzLoc, p, col = cols_tr2[5], add=TRUE)
    }
  }

  if(!is.na(leg))
    title(main = leg, line = 0.5)

  if(!is.na(tag))
    mtext(text = tag, side = 3, line = 2, cex = cex)

  if(!is.na(mex))
    abline(h=mex, lty=1, col=cols[2])

  ## 2. CV profile
  xmod = seq(min(CV),max(CV),length.out = 1000)
  if(class(res) != 'try-error') {
    v   = summary(res)$parameters[,"Estimate"]
    vmod = peak_shape(xmod,v)
  } else {
    v = NA
    vmod = rep(NA,length(xmod))
  }
  ylim = range(c(mMS,vmod), na.rm = TRUE, finite = TRUE)
  plot(
    CV, mMStot,
    type = 'p', pch = 16,
    col  = cols[4],
    xlim = CVlim,
    xlab = 'CV',
    ylim = ylim,
    ylab = 'a.u.'
  )
  title(main = 'Mean CV profile', line = 0.5)

  ## 2.1 Gaussian fit
  if(!any(is.na(v)))
    lines(xmod,vmod,col = cols[2])

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
fit1D <- function(
  mz0, CV0,
  dmz, dCV,
  mz, cv, MS,
  del_mz,
  weighted = FALSE
) {

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # Integrated CV profile
  mMStot = rowSums(MS[, selMz])*del_mz

  # Select CV window
  if(is.na(CV0))
    CV0 = CV[which.max(mMStot)]
  CV1 = CV0 - dCV/2
  CV2 = CV0 + dCV/2
  selCV = CV >= CV1 & CV <= CV2
  CVf = CV[selCV]

  # Refine mz window
  MSloc = MS[selCV, selMz]
  mMS = rowSums(MSloc)*del_mz # Sum over selected mz
  mz_max = which.max(MSloc[which.max(mMS),])
  mz0 = mz[selMz][mz_max]
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz

  selMz  = mz >= mz1 & mz <= mz2 # Select mz area

  # plot(mz[selMz],MS[which.max(mMStot),selMz])

  # Normal fit
  mMS = rowSums(MS[selCV, selMz])*del_mz # Sum over selected mz
  res = try(
    nls(
      mMS ~ k*exp(-1/2*(CVf-mu)^2/sigma^2),
      start = c(
        mu    = CVf[which.max(mMS)],
        sigma = 0.7,
        k     = max(mMS)
      )
    ),
    silent = TRUE
  )

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
  mz, cv, MS,
  del_mz,
  weighted = FALSE
) {

  # Select mz window
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area
  # Integrated CV profile
  mMStot = rowSums(MS[, selMz])*del_mz

  # Select CV window
  if(is.na(CV0))
    CV0 = CV[which.max(mMStot)]
  CV1 = CV0 - dCV/2
  CV2 = CV0 + dCV/2
  selCV = CV >= CV1 & CV <= CV2
  CVf = CV[selCV]

  # Refine mz window
  MSloc = MS[selCV, selMz]
  mMS = rowSums(MSloc)*del_mz # Sum over selected mz
  mz_max = which.max(MSloc[which.max(mMS),])
  mz0 = mz[selMz][mz_max]
  mz1 = mz0 - dmz/2 # min mz for averaging
  mz2 = mz0 + dmz/2 # max mz
  selMz  = mz >= mz1 & mz <= mz2 # Select mz area


  MSloc = MS[selCV, selMz]
  mzLoc = mz[selMz]
  mMS = rowSums(MSloc)*del_mz # Sum over selected mz

  # Normal fit
  grid = expand.grid(x = mzLoc, y = CVf)
  x = grid$x
  y = grid$y
  z = as.vector(t(MSloc))
  weights = rep(1,length(z))
  if(weighted){ # Poisson
    weights = 1 / z
    vm = min(z[z>0])
    weights[!is.finite(weights)] = 1 / vm
  }

  res = try(
    nls(
      z ~ k*exp(-0.5/(1-rho^2) * (
        (x-mx)^2/sx^2 + (y-my)^2/sy^2 -
        2*rho*(x-mx)*(y-my)/(sx*sy))),
      start = c(
        mx  = x[which.max(z)],
        sx  = 0.7,
        my  = y[which.max(z)],
        sy  = 0.4,
        rho = 0,
        k  = max(z)
      ),
      weights = weights
    ),
    silent = TRUE
  )

  return(
    list(
      mz0 = mz0,
      mz1 = mz1,
      mz2 = mz2,
      res = res,
      mMS = mMS,
      mMStot = mMStot,
      CVf = CVf,
      mzLoc = mzLoc,
      MSloc = MSloc
    )
  )
}
getPars1D <- function(res) {
  # Get best params and uncertainty
  v   = summary(res)$parameters[,"Estimate"]   # Best params
  u_v = summary(res)$parameters[,"Std. Error"] # Uncertainty

  # Transform params to quantities of interest
  mu     = v[1]
  u_mu   = u_v[1]
  fwhm   = 2.355 * abs(v[2])
  u_fwhm = 2.355 * u_v[2]
  area   = sqrt(2*pi) * abs(v[2]) * v[3]
  u_area = area * sqrt((u_v[2]/v[2])^2 + (u_v[3]/v[3])^2)
  return(
    list(
      v      = v,
      u_v    = u_v,
      mu     = mu,
      u_mu   = u_mu,
      fwhm   = fwhm,
      u_fwhm = u_fwhm,
      area   = area,
      u_area = u_area
    )
  )
}
getPars2D <- function(res) {
  # Get best params and uncertainty
  v   = summary(res)$parameters[,"Estimate"]   # Best params
  u_v = summary(res)$parameters[,"Std. Error"] # Uncertainty

  # Transform params to quantities of interest
  mu     = v['my']
  u_mu   = u_v['my']
  fwhm   = 2.355 * abs(v['sy'])
  u_fwhm = 2.355 * u_v['sy']
  area   = 2*pi * abs(v['sx'] * v['sy'] * v['k'])
  u_area = area * sqrt((u_v['sx']/v['sx'])^2 +
                       (u_v['sy']/v['sy'])^2 +
                       (u_v['k']/v['k'])^2
                       )

  return(
    list(
      v      = v,
      u_v    = u_v,
      mu     = mu,
      u_mu   = u_mu,
      fwhm   = fwhm,
      u_fwhm = u_fwhm,
      area   = area,
      u_area = u_area
    )
  )

}
getPars = function(res){

  v   = summary(res)$parameters[,"Estimate"]

  if(length(v)==3)
    out = getPars1D(res)
  else
    out = getPars2D(res)

  return(out)
}

# Check sanity of user params ####
file = paste0(dataRepo, tgTable)
assertive::assert_all_are_existing_files(file)

file = paste0(dataRepo, taskTable)
assertive::assert_all_are_existing_files(file)

assertive::assert_is_numeric(fwhm_min)
if(!assertive::is_positive(fwhm_min))
  stop(paste0('Erreur: fwhm_min =',fwhm_min,' should be positive'))

assertive::assert_is_numeric(fwhm_max)
if(!assertive::is_positive(fwhm_max))
  stop(paste0('Erreur: fwhm_max =',fwhm_max,' should be positive'))

assertive::assert_is_numeric(area_min)
if(!assertive::is_positive(area_min))
  stop(paste0('Erreur: area_min =',area_min,' should be positive'))

assertive::assert_is_numeric(dmz)
if(!assertive::is_positive(dmz))
  stop(paste0('Erreur: dmz =',dmz,' should be positive'))

assertive::assert_is_numeric(dCV)
if(!assertive::is_positive(dCV))
  stop(paste0('Erreur: dCV =',dCV,' should be positive'))


# Get targets ####
file = paste0(dataRepo, tgTable)
targets = read.table(
  file = file,
  header = TRUE,
  sep = ';',
  dec = '.',
  check.names = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

empty = rep(NA,nrow(targets))

if(!'CV_ref' %in% colnames(targets))
  targets = cbind(targets,CV_ref=empty)

# Get list of tasks ####
file = paste0(dataRepo, taskTable)
Tasks = read.table(
  file = file,
  header = TRUE,
  sep = ',',
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check that files exist before proceeding
files = paste0(dataRepo,Tasks[,'MS_file'])
assertive::assert_all_are_existing_files(files)

files = paste0(dataRepo,Tasks[,'DMS_file'])
assertive::assert_all_are_existing_files(files)

# Loop over tasks ####
for(task in 1:nrow(Tasks)) {

  msTable = Tasks[task,'MS_file']
  CVTable = Tasks[task,'DMS_file']

  # Build tag
  ## Extract date from CVTable
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
    strsplit(msTable, split='\\.')[[1]][1]
  )

  # Get MS ####
  file = paste0(dataRepo, msTable)
  MS0 = read.table(
    file = file,
    header = FALSE,
    sep = ',',
    stringsAsFactors = FALSE
  )

  time = MS0[, 1]
  range_mz = as.numeric(unlist(strsplit(MS0[1, 7], split = '-')))
  n_del_mz = as.numeric(unlist(strsplit(MS0[1, 8], split = '/')))
  nchan    = n_del_mz[1]
  del_mz   = n_del_mz[2]
  mz       = range_mz[1] + (0:(nchan - 1)) * del_mz
  MS       = as.matrix(MS0[, -(1:8)], ncol = length(mz), byrow = FALSE)

  # Get CV ####
  file = paste0(dataRepo, CVTable)
  CV0 = read.table(
    file = file,
    header = FALSE,
    sep = '\t',
    stringsAsFactors = FALSE
  )
  CV  = rev(CV0[, 4]) # We want increasing CVs

  # tCV = CV0[, 15] / 1000 / 60 # Convertit msec to minutes
  # tCV = rev(tCV - tCV[1])

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
  MS       = apply(MS, 2, rev) # reverse column to conform with CV
  CV   = CV[selCV]
  nCV  = length(CV)

  ## Initialize results table
  resu = cbind(targets,empty,empty,empty,empty,empty,empty,empty)
  colnames(resu) = c(
    colnames(targets),
    'CV','u_CV',
    'FWHM','u_FWHM',
    'Area','u_Area',
    'FitDim'
  )
  xic = cbind(time,rev(CV))
  colnames(xic) = c('time','CV')
  xfi = cbind(time,rev(CV))
  colnames(xfi) = c('time','CV')

  # Loop over targets ####
  for(it in 1:nrow(targets)) {

    mz0 = targets[it,'m/z_exact']
    CV0 = targets[it,'CV_ref']

    if(fit_dim == 2) {
      # 2D fit of peaks
      fitOut = fit2D(
        mz0, CV0,
        dmz, dCV,
        mz, cv, MS,
        del_mz,
        weighted = weighted_fit
      )
      if(class(fitOut$res) == 'try-error' & fallback)
        # 1D fit of peaks
        fitOut = fit1D(
          mz0, CV0,
          dmz, dCV,
          mz, cv, MS,
          del_mz,
          weighted = weighted_fit
        )
    } else {
      # 1D fit of peaks
      fitOut = fit1D(
        mz0, CV0,
        dmz, dCV,
        mz, cv, MS,
        del_mz,
        weighted = weighted_fit
      )
    }

    for (n in names(fitOut))
      assign(n,rlist::list.extract(fitOut,n))

    targets[it,'m/z_exact'] = mz0

    if(class(res)=="try-error") {
      # Fit failed => no fit params
      v    = NA
      mu   = NA
      fwhm = NA
      area = NA
      warning = TRUE
    } else {
      v   = summary(res)$parameters[,"Estimate"]
      peakPars = getPars(res)
      for (n in names(peakPars))
        assign(n,rlist::list.extract(peakPars,n))

      # Quality control
      if(
        filter_results &
        (fwhm <= fwhm_min | fwhm >= fwhm_max | area <= area_min)
      ) {
        warning = TRUE
        # Do not store results

      } else {
        warning = FALSE
        # Store in results table
        resu[it,5]  = signif(mu,4)
        resu[it,6]  = signif(u_mu,2)
        resu[it,7]  = signif(fwhm,3)
        resu[it,8]  = signif(u_fwhm,2)
        resu[it,9]  = signif(area,3)
        resu[it,10] = signif(u_area,2)
        resu[it,11] = ifelse(length(v)==3, 1, 2)
      }
    }

    # Plot data and fit results
    pars = paste0(
      ifelse (warning, '** WARNING **\n','') ,
      'CV = ',      signif(mu,4),
      '\n FWHM = ', signif(fwhm,3),
      '\n Area = ', signif(area,3)
    )
    plotPeak(
      mz, CV, MS,
      fitOut,
      mex = targets[it,'m/z_exact'],
      leg = targets[it,'Name'],
      tag = tag,
      val = pars,
      gPars = gParsLoc
    )
    if(save_figures) {
      png(
        filename = paste0(figRepo, tag, '_', targets[it,1], '.png'),
        width    = 2*gPars$reso,
        height   =   gPars$reso )
      plotPeak(
        mz, CV, MS,
        fitOut,
        mex = targets[it,'m/z_exact'],
        leg = targets[it,'Name'],
        tag = tag,
        val = pars,
        gPars = gPars
      )
      dev.off()
    }

    # Save XIC file
    nam0 = colnames(xic)
    xic = cbind(xic,rev(mMStot))
    colnames(xic) = c(nam0,targets[it,1])

    # Save Fit file
    fit = peak_shape(CV,v)
    nam0 = colnames(xfi)
    xfi = cbind(xfi,rev(fit))
    colnames(xfi) = c(nam0,targets[it,1])

  }

  # Global Heat maps
  mex = targets[,'m/z_exact']
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

  # res = resu[!is.na(resu[,'CV']),]
  print(resu)

  # Save results
  write.csv(resu,file = paste0(tabRepo,tag,'_results.csv'))

  write.csv(xic,file = paste0(tabRepo, tag, '_XIC.csv'))

  write.csv(xfi,file = paste0(tabRepo, tag, '_fit.csv'))

}

# END ####
