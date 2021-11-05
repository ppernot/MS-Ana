#
# Check repeatability of analysis
#
#===============================================
# 2020_07_17 [PP]
# - Adapt to new naming conventions in analysis.R
# - Estimate (weighted) means for properties
# 2020_07_20 [PP]
# - replaced '=' by '_' in userTag (Windows pb.)
# - replaced ':' by '_' in date tag...
# - added ratio vs.dilu fig.
# - added correl. graph for areas AA vs. IS
# 2020_09_24 [PP]
# - changes weighted mean method to Cochran's ANOVA
# - added bargraph of mean ratios
# 2020_09_30 [PP]
# - added the calculation of means for IS species
# 2021_11_03 [PP]
# - added plots
# - management of exported digits
#===============================================
#
## Load packages and functions ####
source('functions.R')

# User configuration params ####
# taskTable = 'list_of_files_AA_2.csv'
# quantTable = 'listaacom_quantification.csv'

taskTable = 'files_quantification_2019July10.csv'
quantTable = 'targets_paper_quantification.csv'

fit_dim = 2
userTag = paste0('fit_dim_',fit_dim)

dataTag = format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

const_fwhm = 0.7
area_min   = 10

makePlots = TRUE

# Check sanity of parameters ####
assertive::assert_all_are_existing_files(dataRepo)
assertive::assert_all_are_existing_files(figRepo)
assertive::assert_all_are_existing_files(tabRepo)

file = paste0(dataRepo, taskTable)
assertive::assert_all_are_existing_files(file)

file = paste0(dataRepo, quantTable)
assertive::assert_all_are_existing_files(file)

# Get tasks list & gather corresponding results ####
Tasks = readTasksFile(paste0(dataRepo, taskTable))
D     = gatherResults(Tasks, tabRepo, userTag)

# Add results columns
D = cbind(
  D,
  ratio = NA, CV_ratio = NA, u_ratio = NA,
  LOD = NA, CV_LOD = NA, u_LOD = NA,
  LOQ = NA, CV_LOQ = NA, u_LOQ = NA,
  slope0 = NA, CV_slope0 = NA, u_slope0 = NA,
  slope = NA, CV_slope = NA, u_slope = NA,
  intercept = NA, CV_intercept = NA, u_intercept = NA,
  R2 = NA, R20 = NA
)

# Get targets ####
quant = readTargetsFile(paste0(dataRepo, quantTable))
targets = quant$Name

# Manage graph labels according to fit_dim
ylab = ifelse(fit_dim==0,'m/z','CV')
rlab = ifelse(fit_dim==0,'...m/z_ref','CV_ref')

if(makePlots) {
  pdf(
    file=paste0(
      figRepo,
      dataTag,
      ifelse(userTag == '', '', '_'),
      userTag,
      '_repAndQuant.pdf'),
    width = 7,
    height = 10)
  par(mfrow=c(3,2),
      mar = mar,
      pty = pty,
      mgp = mgp,
      tcl = tcl)
}

for(it in 1:length(targets)) {
  AA = targets[it]
  IS = quant[it, 'IS']

  selAA = D$Name == AA
  if(sum(selAA) == 0)
    next

  diluAA  = D[selAA, 'dilu']
  aireAA  = D[selAA, 'Area']
  daireAA = D[selAA, 'u_Area']

  selIS = D$Name == IS

  diluIS  = D[selIS, 'dilu']
  aireIS  = D[selIS, 'Area']
  daireIS = D[selIS, 'u_Area']

  if (any(diluAA != diluIS))
    stop('Pb in dilutions consistency')

  cAA = quant[it, 'CAA_Plasma'] +
    quant[it, 'CAA_ref'] * diluAA
  cIS = quant[it, 'CIS_ref']

  ratio = aireAA / aireIS
  dratio = ratio * sqrt(daireAA^2 / aireAA^2 + daireIS^2 / aireIS^2)
  f = formatUncVec(ratio,dratio)
  D[selAA,'ratio']    = f[,'y']
  D[selAA,'u_ratio']  = f[,'uy']
  D[selAA,'CV_ratio'] = signif(100 * abs(dratio/ratio),2)

  # Linear regressions
  io   = order(cAA)
  xo   = cAA[io]
  yo   = ratio[io]
  uyo  = dratio[io]
  wo   = 1/uyo^2
  reg  = lm(yo ~ xo, weights = wo)
  reg0 = lm(yo ~ 0 + xo, weights = wo)

  slope0 = coefficients(reg0)[1]
  uSlope0 = sqrt(diag(vcov(reg0)))
  f = formatUncVec(slope0,uSlope0)
  D[selAA,'slope0']    = f[,'y']
  D[selAA,'u_slope0']  = f[,'uy']
  D[selAA,'CV_slope0'] = signif(100 * abs(uSlope0 / slope0), 2)

  intercept  = coefficients(reg)[1]
  slope      = coefficients(reg)[2]
  uCoefs     = sqrt(diag(vcov(reg)))
  uIntercept = uCoefs[1]
  uSlope     = uCoefs[2]

  f = formatUncVec(slope,uSlope)
  D[selAA,'slope']    = f[,'y']
  D[selAA,'u_slope']  = f[,'uy']
  D[selAA,'CV_slope'] = signif(100 * abs(uSlope / slope),2)

  f = formatUncVec(intercept,uIntercept)
  D[selAA,'intercept']    = f[,'y']
  D[selAA,'u_intercept']  = f[,'uy']
  D[selAA,'CV_intercept'] = signif(100 * abs(uCoefs[1] / intercept),2)

  Sxy = sqrt(sum(residuals(reg)^2)/length(xo-2))
  lod  = 3 * Sxy / slope
  uLod = 3 * Sxy / slope^2 * uSlope
  f = formatUncVec(lod,uLod)
  D[selAA, 'LOD']    = f[,'y']
  D[selAA, 'u_LOD']  = f[,'uy']
  D[selAA, 'CV_LOD'] = D[selAA, 'CV_slope']

  loq  = 10 * Sxy / slope
  uLoq = 10 * Sxy / slope^2 * uSlope
  f = formatUncVec(loq,uLoq)
  D[selAA, 'LOQ']    = f[,'y']
  D[selAA, 'u_LOQ']  = f[,'uy']
  D[selAA, 'CV_LOQ'] = D[selAA, 'CV_slope']

  D[selAA, 'R2']  = signif(summary(reg)$r.squared,2)
  D[selAA, 'R20'] = signif(summary(reg0)$r.squared,2)

  if(makePlots) {
    #Data
    plot(
      xo, yo,
      pch = 16,
      col = cols[5],
      xlab = 'cAA', xlim = c(0,max(xo)),
      ylab = 'aire_AA / aire_IS', ylim = c(0,max(ratio)),
      main = paste(AA,'/',IS)
    )
    grid()
    segments(xo, yo - 2 * uyo,
             xo, yo + 2 * uyo,
             col = cols[5])
    # Fits
    x1 = c(0,xo)
    p = predict(reg,
                newdata = list(xo = x1),
                interval = 'conf')
    matlines(x1[!is.na(p[,1])],p,
             col = cols[4],
             lty=c(1,2,2))

    p0 = predict(reg0,
                 newdata = list(xo = x1),
                 interval = 'conf')
    matlines(x1[!is.na(p0[,1])],p0,
             col = cols[2],
             lty=c(1,2,2))
    # LOD
    Clod = D[selAA, 'LOD']
    abline(v=Clod , col=cols[3], lty = 2)
    mtext('LOD',side=3,col=cols[3],at=Clod, cex=0.75)

    box()
  }
}
if(makePlots)
  dev.off()

# Estimate means ####

props = c("m/z","CV","FWHM_m/z","FWHM_CV","Area",
          "ratio","slope0","slope","intercept","LOD","LOQ")
# signum = c(6   ,3   ,3         ,3        ,3     ,5,
#            5, 5, 5, 5, 5)
if(fit_dim == 0) {
  sel = c(1,3,5:11)
  props = props[sel]
  # signum = signum[sel]
}
# names(signum) = props

meanResTab = D[1:length(targets),]
meanResTab[,] = NA
meanResTab[,'Name'] = targets
rownames(meanResTab) = targets
for(targ in targets) {
  sel = D$Name == targ
  if(sum(sel)==0) next
  for(prop in props) {
    x  = D[sel,prop]
    ux = D[sel,paste0('u_',prop)]
    sel2 = !is.na(x+ux)
    if( sum(sel2) == 0 ) next
    if( sum(sel2) == 1 ) {
      wm  = x[sel2]
      uwm = NA
    } else {
      W   = fwm(x[sel2],ux[sel2])
      wm  = W$wm
      uwm = W$uwm
    }
    f = formatUncVec(wm,uwm)
    meanResTab[targ,prop]              = f[,'y']
    meanResTab[targ,paste0('u_',prop)] = f[,'uy']
    if(! prop %in% props[1:5])
      meanResTab[targ,paste0('CV_',prop)] = signif(100 * abs(uwm/wm),2)
  }
}
meanResTab[,'fit_dim'] = fit_dim
meanResTab[,'tag']     = 'Mean'

# Means for IS species
targIS = quant$IS
propsIS = props[-length(props)]
# signumIS = signum[-length(signum)]
# names(signumIS) = propsIS

meanResTabIS = D[1:length(targIS),]
meanResTabIS[,] = NA
meanResTabIS[,'Name'] = targIS
rownames(meanResTabIS) = targIS
for(targ in targIS) {
  sel = D$Name == targ
  if(sum(sel)==0) next
  for(prop in propsIS) {
    x  = D[sel,prop]
    ux = D[sel,paste0('u_',prop)]
    sel2 = !is.na(x+ux)
    if( sum(sel2) == 0 ) next
    if( sum(sel2) == 1 ) {
      wm  = x[sel2]
      uwm = NA
    } else {
      W   = fwm(x[sel2],ux[sel2])
      wm  = W$wm
      uwm = W$uwm
    }
    f = formatUncVec(wm,uwm)
    meanResTab[targ,prop]              = f[,'y']
    meanResTab[targ,paste0('u_',prop)] = f[,'uy']
    if(! prop %in% props[1:5])
      meanResTabIS[targ,paste0('CV_',prop)] = signif(100 * abs(uwm/wm), 2)
  }
}
meanResTabIS[,'fit_dim'] = fit_dim
meanResTabIS[,'tag']     = 'Mean'

# Save all ####

D = rbind(meanResTab,
          '',
          meanResTabIS,
          '',
          D)

write.csv(
  D,
  row.names = FALSE,
  file = paste0(
    tabRepo,
    dataTag,
    ifelse(userTag == '', '', '_'),
    userTag,
    '_compilation.csv'
  )
)

# END ####
