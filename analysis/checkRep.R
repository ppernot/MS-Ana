#
# Check repeatability of analysis
#
#===============================================
# 2020_07_17 [PP]
# - Adapt to new naming conventions in analysis.R
# - Estimate (weighted) means for properties
#
#===============================================
#
## Load packages and functions ####
source('functions.R')

# User configuration params ####
# taskTable = 'list_of_files_AA_2.csv'
# quantTable = 'listaacom_quantification.csv'

taskTable  = 'list_of_files_Francis_PP.csv'
quantTable = 'targets_paper_quantification.csv'

fit_dim = 2
userTag = paste0('fit_dim=',fit_dim)

const_fwhm = 0.7

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
# Add ratio columns
D = cbind(D, ratio = NA, u_ratio = NA)

# Get targets ####
quant = readTargetsFile(paste0(dataRepo, quantTable))
targets = quant$Name

# Manage graph labels according to fit_dim
ylab = ifelse(fit_dim==0,'m/z','CV')
rlab = ifelse(fit_dim==0,'m/z_exact','CV_ref')

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
  dratio = ratio * sqrt(
    daireAA^2/aireAA^2 +
      daireIS^2/aireIS^2
  )
  D[selAA,'ratio']   = signif(ratio,5)
  D[selAA,'u_ratio'] = signif(dratio,2)

  if(makePlots) {
    par(mfrow = c(2,3), mar = mar)

    # CV vs. dilution
    x = D[selAA, 'dilu']
    y = D[selAA, ylab]
    dy = 2 * D[selAA, paste0('u_',ylab)]
    y0 = D[selAA, rlab]

    plot(
      x,
      y,
      pch = 16,
      col = cols[5],
      log = '',
      xlab = 'Dilution',
      ylab = ylab,
      main = AA
    )
    segments(x, y - dy, x, y + dy, col = cols[5])
    lines(x, y0, lwd = 2, col = cols[3])
    text(min(x),y0[1],'CV_ref',adj=0,cex=0.75)
    cm = sort(unique(x))
    ym = c()
    for (i in 1:length(cm))
      ym[i] = mean(y[x == cm[i]], na.rm = TRUE)
    lines(cm, ym, lty = 2, col = cols[2])
    grid()

    # FWHM vs. dilution
    ftag = ifelse(fit_dim==0,'FWHM_m/z','FWHM_CV')
    x = D[selAA, 'dilu']
    y = D[selAA, ftag]
    dy = 2 * D[selAA, paste0('u_',ftag)]
    plot(
      x,
      y,
      pch = 16,
      col = cols[5],
      log = '',
      xlab = 'Dilution',
      ylab = ftag,
      main = AA
    )
    segments(x, y - dy, x, y + dy, col = cols[5])
    rect(0.8*min(x),0.8*const_fwhm,
         1.2*max(x),1.2*const_fwhm,
         col = cols_tr[4], border=NA)
    grid()

    # Area vs. dilution
    x = D[selAA, 'dilu']
    y = D[selAA, 'Area']
    dy = 2 * D[selAA, 'u_Area']
    plot(
      x,
      y,
      pch = 16,
      col = cols[5],
      log = '',
      xlab = 'Dilution',
      ylab = 'Area',
      main = AA
    )
    segments(x, y - dy, x, y + dy, col = cols[5])
    abline(h=area_min, lwd= 2, col=cols[2])
    io = order(x)
    xo = x[io]
    yo = y[io]
    reg = lm(yo~xo, weights = 1/(dy[io]/2)^2)
    p = predict(reg, interval = 'conf')
    matlines(xo[!is.na(yo)],p,
             col = cols[4],
             lty=c(1,2,2))
    grid()

    # CV vs. dilution
    x = D[selIS, 'dilu']
    y = D[selIS, ylab]
    dy = 2 * D[selIS, paste0('u_',ylab)]
    y0 = D[selIS, rlab]

    plot(
      x,
      y,
      pch = 16,
      col = cols[5],
      log = '',
      xlab = 'Dilution',
      ylab = ylab,
      main = IS
    )
    segments(x, y - dy, x, y + dy, col = cols[5])
    lines(x, y0, lwd = 2, col = cols[3])
    text(min(x),y0[1],'CV_ref',adj=0,cex=0.75)
    cm = sort(unique(x))
    ym = c()
    for (i in 1:length(cm))
      ym[i] = mean(y[x == cm[i]], na.rm = TRUE)
    lines(cm, ym, lty = 2, col = cols[2])
    grid()

    # FWHM vs. dilution
    x = D[selIS, 'dilu']
    y = D[selIS, ftag]
    dy = 2 * D[selIS, paste0('u_',ftag)]
    plot(
      x,
      y,
      pch = 16,
      col = cols[5],
      log = '',
      xlab = 'Dilution',
      ylab = 'FWHM',
      main = IS
    )
    segments(x, y - dy, x, y + dy, col = cols[5])
    rect(0.8*min(x),0.8*const_fwhm,
         1.2*max(x),1.2*const_fwhm,
         col = cols_tr[4], border=NA)
    grid()

    # Area vs. dilution
    x = D[selIS, 'dilu']
    y = D[selIS, 'Area']
    dy = 2 * D[selIS, 'u_Area']
    plot(
      x,
      y,
      pch = 16,
      col = cols[5],
      log = '',
      xlab = 'Dilution',
      ylab = 'Area',
      main = IS
    )
    segments(x, y - dy, x, y + dy, col = cols[5])
    abline(h=area_min, lwd= 2, col=cols[2])
    io = order(x)
    xo = x[io]
    yo = y[io]
    reg = lm(yo~xo, weights = 1/(dy[io]/2)^2)
    p = predict(reg, interval = 'conf')
    matlines(xo[!is.na(yo)],p,
             col = cols[4],
             lty=c(1,2,2))
    grid()
  }
}

# Estimate means ####

props = c("m/z","CV","FWHM_m/z","FWHM_CV","Area","ratio")
signum = c(6   ,3   ,3         ,3        ,3     ,5)
if(fit_dim == 0) {
  sel = c(1,3,5,6)
  props = props[sel]
  signum = signum[sel]
}
names(signum) = props

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
      uwn = NA
    } else {
      x   = x[sel2]
      ux  = ux[sel2]
      w   = 1/ux^2
      w   = w / sum(w)
      wm  = weighted.mean(x,w)
      s2  = var(x)
      uwm = sqrt( sum(w^2*(s2+ux^2)) )
    }
    meanResTab[targ,prop]              = signif(wm , signum[prop])
    meanResTab[targ,paste0('u_',prop)] = signif(uwm, 2)
  }
}
meanResTab[,'fit_dim'] = fit_dim
meanResTab[,'tag']     = 'Mean'

D = rbind(meanResTab,
          '',
          D)

dataTag = format(Sys.time(), "%Y_%m_%d_%H:%M:%S")
write.csv(
  D[, -1],
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