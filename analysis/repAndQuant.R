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
dataRepo  = "../data/TestQuant/"

taskTable  = 'Files_quanti_Gamme4.csv'
quantTable = 'targets_OA1_quantification.csv'

fit_dim = 2
userTag = paste0('fit_dim_',fit_dim)

dataTag = format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

const_fwhm = 0.7
area_min   = 10

makePlots = TRUE

serialize = TRUE # Serialize data on a day-by-day basis

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

# Data structure
if (serialize) {
  dayTags = sapply(Tasks[,"DMS_file"],getDayTag)
  names(dayTags)=NULL
}

# Add results columns
D = cbind(
  D,
  ratio = NA, CV_ratio = NA, u_ratio = NA,
  LOD = NA, CV_LOD = NA, u_LOD = NA,
  LOQ = NA, CV_LOQ = NA, u_LOQ = NA
  # slope0 = NA, CV_slope0 = NA, u_slope0 = NA,
  # slope = NA, CV_slope = NA, u_slope = NA,
  # intercept = NA, CV_intercept = NA, u_intercept = NA,
  # R2 = NA, R20 = NA
)

# Get targets ####
quant   = readTargetsFile(paste0(dataRepo, quantTable))
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

  # LOD / LOQ ####

  Data    = data.frame(x = cAA, y = ratio, t = dayTags )
  res     = lodCal(Data, serialize)

  tabRes  = cbind(AA,res$tab)
  colnames(tabRes) = c('AA',colnames(res$tab))
  file = paste0(
    tabRepo,
    dataTag,
    ifelse(userTag == '', '', '_'),
    userTag,
    '_LOD.csv'
  )
  write.table(
    tabRes,
    file = file,
    row.names = FALSE,
    sep = ',',
    col.names = it==1,
    append = TRUE
  )

  f = formatUncVec(res$LOD,res$LOD.u)
  D[selAA, 'LOD']    = f[,'y']
  D[selAA, 'u_LOD']  = f[,'uy']
  D[selAA, 'CV_LOD'] = signif(100 * abs(res$LOD.u / res$LOD),2)

  f = formatUncVec(res$LOQ,res$LOQ.u)
  D[selAA, 'LOQ']    = f[,'y']
  D[selAA, 'u_LOQ']  = f[,'uy']
  D[selAA, 'CV_LOQ'] = signif(100 * abs(res$LOQ.u / res$LOQ),2)

  if(makePlots) {
    for(iS in seq_along(res$series)) {
      tag = res$series[iS]

      x = res$fit[[tag]]$model$x
      y = res$fit[[tag]]$model$y

      if (iS == 1) {
        plot(
          x,
          y,
          type = 'p',
          pch = 16,
          col = 4,
          xlab = 'cAA',
          xlim = c(0, max(cAA)),
          ylab = 'aire_AA / aire_IS',
          ylim = c(0, max(ratio)),
          # xaxs = 'i',
          yaxs = 'i',
          main = AA
        )
        grid()
        abline(v = 0, h = 0)
      } else {
        points(x, y, pch = 16, col = 4)
      }
      abline(res$fit[[tag]], lty = 1, col = 'gray50')
      abline(v = res$tab[, 'LOD'], lty = 2, col = 2)
      abline(
        v = res$LOD,
        lty = 1,
        col = 2,
        lwd = 1.5
      )
      segments(
        res$LOD - 1.96 * res$LOD.u,
        0,
        res$LOD + 1.96 * res$LOD.u,
        0,
        col = 2,
        lwd = 4,
        lend = 2
      )
      segments(
        res$LOD - 1.96 * res$LOD.u,
        max(ratio),
        res$LOD + 1.96 * res$LOD.u,
        max(ratio),
        col = 2,
        lwd = 4,
        lend = 2
      )
      mtext('LOD', side = 3, at = res$LOD, col = 2, cex = 0.5)
      mtext(signif(res$LOD,3), side = 1, at = res$LOD, col = 2, cex = 0.5)

      # abline(
      #   v = res$LOQ,
      #   lty = 1,
      #   col = 3,
      #   lwd = 1
      # )
      # segments(
      #   res$LOQ - 1.96 * res$LOQ.u,
      #   0,
      #   res$LOQ + 1.96 * res$LOQ.u,
      #   0,
      #   col = 3,
      #   lwd = 4,
      #   lend = 2
      # )
    }

  }

}
if(makePlots)
  dev.off()

# Estimate means ####

# props = c("m/z","CV","FWHM_m/z","FWHM_CV","Area","ratio","LOD","LOQ")
# if(fit_dim == 0) {
#   sel = c(1,3,5:11)
#   props = props[sel]
# }
#
# meanResTab = D[1:length(targets),]
# meanResTab[,] = NA
# meanResTab[,'Name'] = targets
# rownames(meanResTab) = targets
# for(targ in targets) {
#   sel = D$Name == targ
#   if(sum(sel)==0) next
#   for(prop in props) {
#     x  = D[sel,prop]
#     ux = D[sel,paste0('u_',prop)]
#     sel2 = !is.na(x+ux)
#     if( sum(sel2) == 0 ) next
#     if( sum(sel2) == 1 ) {
#       wm  = x[sel2]
#       uwm = NA
#     } else {
#       W   = fwm(x[sel2],ux[sel2])
#       wm  = W$wm
#       uwm = W$uwm
#     }
#     f = formatUncVec(wm,uwm)
#     meanResTab[targ,prop]              = f[,'y']
#     meanResTab[targ,paste0('u_',prop)] = f[,'uy']
#     if(! (prop %in% props[1:5]) & ! (prop == 'LOD1'))
#       meanResTab[targ,paste0('CV_',prop)] = signif(100 * abs(uwm/wm),2)
#   }
# }
# meanResTab[,'fit_dim'] = fit_dim
# meanResTab[,'tag']     = 'Mean'
#
# # Means for IS species
# targIS = quant$IS
# propsIS = props[-length(props)]
#
# meanResTabIS = D[1:length(targIS),]
# meanResTabIS[,] = NA
# meanResTabIS[,'Name'] = targIS
# rownames(meanResTabIS) = targIS
# for(targ in targIS) {
#   sel = D$Name == targ
#   if(sum(sel)==0) next
#   for(prop in propsIS) {
#     x  = D[sel,prop]
#     ux = D[sel,paste0('u_',prop)]
#     sel2 = !is.na(x+ux)
#     if( sum(sel2) == 0 ) next
#     if( sum(sel2) == 1 ) {
#       wm  = x[sel2]
#       uwm = NA
#     } else {
#       W   = fwm(x[sel2],ux[sel2])
#       wm  = W$wm
#       uwm = W$uwm
#     }
#     f = formatUncVec(wm,uwm)
#     meanResTab[targ,prop]              = f[,'y']
#     meanResTab[targ,paste0('u_',prop)] = f[,'uy']
#     if(! (prop %in% props[1:5]) & ! (prop == 'LOD1'))
#       meanResTabIS[targ,paste0('CV_',prop)] = signif(100 * abs(uwm/wm), 2)
#   }
# }
# meanResTabIS[,'fit_dim'] = fit_dim
# meanResTabIS[,'tag']     = 'Mean'


# D = rbind(meanResTab,
#           '',
#           meanResTabIS,
#           '',
#           D)

# Save all ####

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
