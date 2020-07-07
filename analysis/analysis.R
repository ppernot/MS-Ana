# User configuration params ####
taskTable = 'list_of_files_Viet_total.csv'
tgTable = 'list_of_targets_Plasma_Viet.csv'
filter_results = TRUE
fwhm_min = 0.5
fwhm_max = 1.5
area_min = 10
save_figures = TRUE
plot_maps = FALSE
weighted_fit = FALSE
refine_CV0 = TRUE
const_fwhm = 0.7
fit_dim = 2     # Fit 2D peaks
fallback = TRUE # Fallback on 1D fit if 2D fails
dmz = 1.0       # Width of mz window around
                # exact mz for signal averaging
dCV = 1.2       # Width of CV window around
                # reference CV for peak fit

# Code Setup ####
options(warn=0)

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

## Load functions
source('functions.R')

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

# Check sanity of parameters ####
assertive::assert_all_are_existing_files(dataRepo)
assertive::assert_all_are_existing_files(figRepo)
assertive::assert_all_are_existing_files(tabRepo)

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

# Get tasks list ####
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
        mz, CV, MS,
        del_mz,
        weighted = weighted_fit,
        refine_CV0 = refine_CV0,
        const_fwhm = const_fwhm
      )
      if(class(fitOut$res) == 'try-error' & fallback)
        # 1D fit of peaks
        fitOut = fit1D(
          mz0, CV0,
          dmz, dCV,
          mz, CV, MS,
          del_mz,
          weighted = weighted_fit,
          refine_CV0 = refine_CV0,
          const_fwhm = const_fwhm
        )
    } else {
      # 1D fit of peaks
      fitOut = fit1D(
        mz0, CV0,
        dmz, dCV,
        mz, CV, MS,
        del_mz,
        weighted = weighted_fit,
        refine_CV0 = refine_CV0,
        const_fwhm = const_fwhm
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
  if(plot_maps) {
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
  }

  res = resu[!is.na(resu[,'CV']),]
  print(res)

  # Save results
  write.csv(resu,file = paste0(tabRepo,tag,'_results.csv'))

  write.csv(xic,file = paste0(tabRepo, tag, '_XIC.csv'))

  write.csv(xfi,file = paste0(tabRepo, tag, '_fit.csv'))

}

# END ####
