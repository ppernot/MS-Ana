#===============================================
# 2020_07_16 [PP]
# - Added 2 new colums to results table
#   to store m/z and u_m/z from 2D fit
# - Integrate fast method (fit_dim = 0)
# - Added data path management
# - Added optional tag (userTag):
#   presently defined by fit_dim value,
#   to avoid overwriting of results files
#   when trying different fit_dim options
# 2020_07_20 [PP]
# - replaced '=' by '_' in userTag (Windows pb.)
# - reparameterized gaussians with area replacing height
#   (avoids covariances in estimation of u_area )
# - suppressed 'rho' param in 2D gaussians
# 2020_07_21 [PP]
# - corrected typo in calculation of u_area in getPars1D
# 2020_07_24 [PP]
# - save ctrl params as metadata
#===============================================
#
## Load packages and functions ####
source('functions.R')

# User configuration params ####

taskTable = 'files_quantification_2019July10.csv'
tgTable   = 'targets_paper.csv'

filter_results = TRUE
fwhm_mz_min = 0.1
fwhm_mz_max = 0.5
fwhm_cv_min = 0.5
fwhm_cv_max = 1.5
area_min    = 10

save_figures = TRUE
plot_maps    = FALSE

fit_dim  = 2    # 2: fit 2D peaks; 1: fit 1D CV line; 0: fit 1D m/z line
fallback = TRUE # Fallback on fit_dim=1 fit if 2D fit fails

weighted_fit = FALSE
refine_CV0   = TRUE
const_fwhm   = ifelse(fit_dim == 0,NA,0.7)

dmz = 1.0       # Width of mz window around
                # exact mz for signal averaging
dCV = 1.2       # Width of CV window around
                # reference CV for peak fit

debug = FALSE    # Stops after first task

userTag = paste0('fit_dim_',fit_dim)

ctrlParams = list(
  userTag        = userTag,
  taskTable      = taskTable,
  tgTable        = tgTable,
  filter_results = filter_results,
  fwhm_mz_min    = fwhm_mz_min,
  fwhm_mz_max    = fwhm_mz_max,
  fwhm_cv_min    = fwhm_cv_min,
  fwhm_cv_max    = fwhm_cv_max,
  area_min       = area_min,
  fit_dim        = fit_dim,
  fallback       = fallback,
  weighted_fit   = weighted_fit,
  refine_CV0     = refine_CV0,
  const_fwhm     = const_fwhm,
  dmz            = dmz,
  dCV            = dCV
)

# Check sanity of parameters ####
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
  stop(paste0('Erreur: area_min =',area_min,' should be positive'))

assertive::assert_is_numeric(dmz)
if(!assertive::is_positive(dmz))
  stop(paste0('Erreur: dmz =',dmz,' should be positive'))

assertive::assert_is_numeric(dCV)
if(!assertive::is_positive(dCV))
  stop(paste0('Erreur: dCV =',dCV,' should be positive'))


# Get targets ####
targets = readTargetsFile(paste0(dataRepo, tgTable))
empty = rep(NA,nrow(targets))
if(!'CV_ref' %in% colnames(targets))
  targets = cbind(targets,CV_ref=empty)

# Get tasks list ####
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

# Loop over tasks ####
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
  tag = makeTag(CVTable, msTable, userTag)

  # Get MS ####
  file = paste0(dataRepo, dataPath, msTable)
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
  MS       = as.matrix(MS0[, -(1:8)],
                       ncol = length(mz),
                       byrow = FALSE)

  # Get CV ####
  file = paste0(dataRepo, CVTable)
  CV0 = read.table(
    file = file,
    header = FALSE,
    sep = '\t',
    stringsAsFactors = FALSE
  )
  CV  = rev(CV0[, 4]) # We want increasing CVs

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
  xic = cbind(time,rev(CV))
  colnames(xic) = c('time','CV')
  if( fit_dim == 0) {
    xfi = cbind(1:length(mz),mz)
    colnames(xfi) = c('index','m/z')
  } else {
    xfi = cbind(time,rev(CV))
    colnames(xfi) = c('time','CV')
  }

  # Loop over targets ####
  for( it in 1:nrow(targets) ) {

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
      dimfit = 2
      if(class(fitOut$res) == 'try-error' & fallback) {
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
        dimfit = 1
      }

    } else if (fit_dim == 1) {
      # 1D fit of peaks; fixed m/z
      fitOut = fit1D(
        mz0, CV0,
        dmz, dCV,
        mz, CV, MS,
        del_mz,
        weighted = weighted_fit,
        refine_CV0 = refine_CV0,
        const_fwhm = const_fwhm
      )
      dimfit = 1

    } else {
      # 1D fit of peak; fixed CV
      fitOut = fit1D_MS(
        mz0, CV0,
        dmz, dCV,
        mz, CV, MS,
        weighted = weighted_fit,
        const_fwhm = const_fwhm
      )
      dimfit = 0
    }

    for (n in names(fitOut))
      assign(n,rlist::list.extract(fitOut,n))

    targets[it,'m/z_exact'] = mz0

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
      ifelse(dimfit > 0,
             '',
             paste0('m/z = ', signif(mzopt,6),'\n')),
      ifelse(dimfit == 0,
             '',
             paste0('CV = ', signif(cvopt,4),'\n')),
      ifelse(dimfit > 0,
             paste0('FWHM = ', signif(fwhm_cv,3),'\n'),
             paste0('FWHM = ', signif(fwhm_mz,3),'\n')),
      'Area = ', signif(area,3)
    )
    plotPeak(
      mz, CV, MS,
      fitOut,
      mex = targets[it,'m/z_exact'],
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
        mex = targets[it,'m/z_exact'],
        leg = targets[it,'Name'],
        tag = tag,
        val = pars,
        type = ifelse(dimfit==0,'m/z','CV'),
        CV0 = CV0,
        gPars = gPars
      )
      dev.off()
    }

    # Save XIC file
    nam0 = colnames(xic)
    xic = cbind(xic,rev(mMStot))
    colnames(xic) = c(nam0,targets[it,1])

    # Save Fit file
    if(fit_dim == 0)
      fit = peak_shape(mz, v)
    else
      fit = peak_shape(CV, v)
    nam0 = colnames(xfi)
    xfi = cbind(xfi,rev(fit))
    colnames(xfi) = c(nam0,targets[it,1])

    # if(debug) break
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

  # res = resu[!is.na(resu[,'CV']),]
  # print(res)

  # Save results
  write.csv(resu,file = paste0(tabRepo, tag, '_results.csv'))
  write.csv(xic,file  = paste0(tabRepo, tag, '_XIC.csv'))
  write.csv(xfi,file  = paste0(tabRepo, tag, '_fit.csv'))

  # Metadata
  rlist::list.save(
    ctrlParams,
    file = file.path(
      tabRepo,
      paste0(tag,'_ctrlParams.yaml')
    )
  )

  if(debug) stop()
}

# END ####
