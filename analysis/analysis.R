# User configuration params ####
taskTable = 'pm.csv'
tgTable = 'targets.csv'
filter_results = FALSE
fwhm_min = 0.5
fwhm_max = 1.5
area_min = 10
save_figures = TRUE
dmz = 0.5 # width of mz window around exact mz for signal averaging
dCV = 2 # width of CV window around reference CV for peak fit

# Code Setup ####
# Install packages if necessary
## CRAN packages
libs <- c("xtable","mixtools","inlmisc","rlist","repmis")
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    install.packages(
      lib,
      dependencies = TRUE,
      repos = "https://cran.univ-paris1.fr"
    )
  }
}
## Load packages and generate biblio
repmis::LoadandCite(libs,file='../article/packages.bib')

# Set graphical params
## For PNG figures
gPars = list(
  cols     = rev(inlmisc::GetColors(8))[1:7],
  cols_tr  = rev(inlmisc::GetColors(8, alpha = 0.2))[1:7],
  cols_tr2 = rev(inlmisc::GetColors(8, alpha = 0.5))[1:7],
  pty      = 's',
  mar      = c(3,3,1.6,.5),
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
  p[3]*exp(-1/2*(x-p[1])^2/p[2]^2)
}
plotPeak = function(
  mz, CV, MS, CVf, mMS, vg,
  mex = NA,
  leg = '',
  res = NA,
  model = peak_shape,
  mzlim = range(mz),
  CVlim = range(CV),
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

  ## 1. image
  sel1 = apply(cbind(CV-CVlim[1],CV-CVlim[2]),1,prod) <= 0
  sel2 = apply(cbind(mz-mzlim[1],mz-mzlim[2]),1,prod) <= 0

  image(
    CV[sel1], mz[sel2], MS[sel1,sel2],
    xlim = CVlim,
    xlab = 'CV',
    ylim = mzlim,
    ylab = 'm/z',
    main = leg
  )
  grid()
  if(!is.na(mex))
    abline(h=mex,lty=1,col=cols[2])

  ## 2. CV profile
  xmod = seq(min(CVf),max(CVf),length.out = 1000)
  vmod = model(xmod,vg)
  plot(
    CVf,
    mMS,
    type = 'p', pch = 16,
    col  = cols[4],
    xlim = CVlim,
    xlab = 'CV',
    ylim = range(c(mMS,vmod)),
    ylab = 'a.u.',
    main = 'Mean CV profile'
  )
  ## 2.1 Gaussian fit
  lines(xmod,vmod,col = cols[2])

  ## Add fit results
  if (!is.na(res))
    legend(
      'topleft',
      title  = res,
      legend = '',
      inset  = 0.12,
      bty    = 'n',
      cex    = 0.8)

}

# Get targets ####
targets = read.table(
  file = paste0(dataRepo, tgTable),
  header = TRUE,
  sep = '\t',
  dec = '.',
  check.names = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

empty = rep(NA,nrow(targets))

if(!'CV_ref' %in% colnames(targets))
  targets = cbind(targets,CV_ref=empty)

# Get list of tasks ####
Tasks = read.table(
  file = paste0(dataRepo, taskTable),
  header = TRUE,
  sep = ',',
  check.names = FALSE,
  stringsAsFactors = FALSE
)

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
  MS0 = read.table(
    file = paste0(dataRepo,msTable),
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
  MS       = as.matrix(MS0[, -(1:8)],ncol=length(mz),byrow=FALSE)
  MS       = apply(MS,2,rev) # reverse column to conform with CV

  # Get CV ####
  CV0 = read.table(
    file = paste0(dataRepo, CVTable),
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
  CV   = CV[selCV]
  nCV  = length(CV)

  ## Initialize results table
  resu = cbind(targets,empty,empty,empty,empty,empty,empty)
  colnames(resu) = c(
    colnames(targets),
    'CV','u_CV',
    'FWHM','u_FWHM',
    'Area','u_Area'
  )
  xic = cbind(time,rev(CV))
  colnames(xic) = c('time','CV')
  xfi = cbind(time,rev(CV))
  colnames(xfi) = c('time','CV')

  # Loop over targets ####
  for(it in 1:nrow(targets)) {

    # Select mz window
    mz0 = targets[it,'m/z_exact']
    mz1 = mz0 - dmz/2 # min mz for averaging
    mz2 = mz0 + dmz/2 # max mz

    # Select CV window
    CV0 = targets[it,'CV_ref']
    if(is.na(CV0)) {
      # Use full range
      selCV = 1:nCV
      CVf = CV
    } else {
      CV1 = CV0 - dCV/2
      CV2 = CV0 + dCV/2
      selCV = CV >= CV1 & CV <= CV2
      CVf = CV[selCV]
    }

    # Estimate Profile
    selMz  = mz >= mz1 & mz <= mz2 # Select mz area
    mMS    = rowSums(MS[selCV, selMz]) # sum over selected mz
    mMStot = rowSums(MS[, selMz]) # full CV range for output in XIC file

    # Normal fit
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
    if(class(res)=="try-error")
      next # Skip results and figs

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

    # Quality control
    if(
      filter_results &
      (fwhm <= fwhm_min | fwhm >= fwhm_max | area <= area_min)
    )
      next # Leave 'it' line of table results with NAs

    # Store in results table
    resu[it,5] = signif(mu,4)
    resu[it,6] = signif(u_mu,2)
    resu[it,7] = signif(fwhm,3)
    resu[it,8] = signif(u_fwhm,2)
    resu[it,9] = signif(area,3)
    resu[it,10] = signif(u_area,2)

    # Plot data and fit results
    plotPeak(
      mz, CV, MS, CVf, mMS,
      vg = v,
      mex = targets[it,'m/z_exact'],
      leg = targets[it,'Name'],
      res = paste0(
        'CV = ',      signif(mu,4),
        '\n FWHM = ', signif(fwhm,3),
        '\n Area = ', signif(area,3)
      ),
      mzlim = c(mz1,mz2),
      CVlim = range(CVf),
      gPars = gParsLoc
    )

    if(save_figures) {
      png(
        filename = paste0(figRepo, tag, '_', targets[it,1], '.png'),
        width    = 2*gPars$reso,
        height   =   gPars$reso )
      plotPeak(
        mz, CV, MS, CVf, mMS,
        vg = v,
        mex = targets[it,'m/z_exact'],
        leg = targets[it,'Name'],
        res = paste0(
          'CV = ',      signif(mu,4),
          '\n FWHM = ', signif(fwhm,3),
          '\n Area = ', signif(area,3)
        ),
        mzlim = c(mz1,mz2),
        CVlim = range(CV),
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

  # res = resu[!is.na(resu[,'CV']),]
  print(resu)
  # Save results
  write.csv(resu,file = paste0(tabRepo,tag,'_results.csv'))

  write.csv(xic,file = paste0(tabRepo, tag, '_XIC.csv'))

  write.csv(xfi,file = paste0(tabRepo, tag, '_fit.csv'))

}

# END ####
